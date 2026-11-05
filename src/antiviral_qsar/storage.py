from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import uuid
import fcntl

from .schema import Measurement


def now():
    return datetime.now(timezone.utc).isoformat()


class Store:
    def __init__(self, root: Path):
        self.root = root
        (root / "data").mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(root / "data/state.sqlite")
        self.db.row_factory = sqlite3.Row
        self.db.executescript("""
        PRAGMA journal_mode=WAL;
        PRAGMA foreign_keys=ON;
        CREATE TABLE IF NOT EXISTS runs (
          id TEXT PRIMARY KEY, created TEXT, updated TEXT, status TEXT, virus TEXT,
          options TEXT, inventory TEXT DEFAULT '{}', report TEXT DEFAULT '{}');
        CREATE TABLE IF NOT EXISTS jobs (
          run_id TEXT, source TEXT, kind TEXT, key TEXT, payload TEXT,
          status TEXT DEFAULT 'pending', attempts INTEGER DEFAULT 0, error TEXT DEFAULT '',
          updated TEXT, PRIMARY KEY(run_id,source,kind,key));
        CREATE TABLE IF NOT EXISTS observations (
          run_id TEXT, id TEXT, data TEXT, PRIMARY KEY(run_id,id));
        CREATE TABLE IF NOT EXISTS events (
          run_id TEXT, created TEXT, level TEXT, message TEXT);
        CREATE TABLE IF NOT EXISTS llm_cache (
          key TEXT PRIMARY KEY, created TEXT, value TEXT);
        CREATE TABLE IF NOT EXISTS quality_gates (
          signature TEXT PRIMARY KEY, data TEXT, updated TEXT);
        CREATE INDEX IF NOT EXISTS jobs_source_key_status ON jobs(source,key,status);
        """)

    def new_run(self, virus: dict, options: dict):
        rid = datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        self.db.execute("INSERT INTO runs(id,created,updated,status,virus,options) VALUES(?,?,?,?,?,?)",
                        (rid, now(), now(), "discovering", json.dumps(virus), json.dumps(options)))
        self.db.commit()
        (self.root / "data/runs" / rid).mkdir(parents=True, exist_ok=True)
        return rid

    def run(self, rid):
        row = self.db.execute("SELECT * FROM runs WHERE id=?", (rid,)).fetchone()
        if not row:
            raise ValueError(f"Run {rid} was not found")
        result = dict(row)
        for name in ("virus", "options", "inventory", "report"):
            result[name] = json.loads(result[name])
        return result

    def update_run(self, rid, **fields):
        allowed = {"status", "options", "inventory", "report"}
        if set(fields) - allowed:
            raise ValueError("Unsupported run fields")
        fields["updated"] = now()
        values = [json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v
                  for v in fields.values()]
        self.db.execute("UPDATE runs SET " + ",".join(f"{k}=?" for k in fields) + " WHERE id=?",
                        (*values, rid))
        self.db.commit()

    def add_job(self, rid, source, kind, key, payload):
        self.db.execute("INSERT OR IGNORE INTO jobs(run_id,source,kind,key,payload,updated) VALUES(?,?,?,?,?,?)",
                        (rid, source, kind, str(key), json.dumps(payload, ensure_ascii=False), now()))
        self.db.commit()

    def add_jobs(self, rid, source, kind, entries):
        timestamp = now()
        with self.db:
            self.db.executemany("INSERT OR IGNORE INTO jobs(run_id,source,kind,key,payload,updated) VALUES(?,?,?,?,?,?)",
                               ((rid, source, kind, str(key), json.dumps(payload, ensure_ascii=False), timestamp)
                                for key, payload in entries))

    def jobs(self, rid, kind=None, unfinished=False):
        sql, params = "SELECT * FROM jobs WHERE run_id=?", [rid]
        if kind:
            sql += " AND kind=?"
            params.append(kind)
        if unfinished:
            sql += " AND status NOT IN ('done','skipped')"
        rows = self.db.execute(sql + " ORDER BY rowid", params).fetchall()
        return [{**dict(r), "payload": json.loads(r["payload"])} for r in rows]

    def start_job(self, job):
        self.db.execute("UPDATE jobs SET status='running',attempts=attempts+1,error='',updated=? "
                        "WHERE run_id=? AND source=? AND kind=? AND key=?",
                        (now(), job["run_id"], job["source"], job["kind"], job["key"]))
        self.db.commit()

    def finish_job(self, job, observations=(), status="done", error=""):
        with self.db:
            for m in observations:
                self.db.execute("INSERT OR REPLACE INTO observations VALUES(?,?,?)",
                                (job["run_id"], m.observation_id, m.model_dump_json()))
            self.db.execute("UPDATE jobs SET status=?,error=?,updated=? "
                            "WHERE run_id=? AND source=? AND kind=? AND key=?",
                            (status, error[:2000], now(), job["run_id"], job["source"], job["kind"], job["key"]))

    def observations(self, rid):
        return [Measurement.model_validate_json(r[0]) for r in self.db.execute(
            "SELECT data FROM observations WHERE run_id=? ORDER BY id", (rid,))]

    def event(self, rid, message, level="info"):
        self.db.execute("INSERT INTO events VALUES(?,?,?,?)", (rid, now(), level, message))
        self.db.commit()

    def cache_get(self, key):
        r = self.db.execute("SELECT value FROM llm_cache WHERE key=?", (key,)).fetchone()
        return json.loads(r[0]) if r else None

    def cache_put(self, key, value):
        self.db.execute("INSERT OR REPLACE INTO llm_cache VALUES(?,?,?)", (key, now(), json.dumps(value)))
        self.db.commit()

    def previously_processed(self, source, key):
        return self.db.execute("SELECT COUNT(*) FROM jobs WHERE source=? AND key=? AND status='done'",
                               (source, str(key))).fetchone()[0] > 0

    def reuse_paper(self, rid, key, paper, signature):
        rows = self.db.execute("SELECT j.run_id,r.report FROM jobs j JOIN runs r ON r.id=j.run_id "
                               "WHERE j.source='europepmc' AND j.key=? AND j.status='done' "
                               "AND j.run_id<>? ORDER BY j.updated DESC", (key, rid)).fetchall()
        for row in rows:
            if json.loads(row['report']).get('model', {}).get('signature') != signature:
                continue
            previous = [m for m in self.observations(row['run_id'])
                        if m.source == 'europepmc' and m.document_id == paper.get('pmcid')]
            job = {'run_id':rid, 'source':'europepmc', 'kind':'paper', 'key':key}
            self.finish_job(job, previous)
            self.event(rid, f"reused_document:{key}:from:{row['run_id']}")
            return True
        return False

    def quality_passed(self, signature):
        r = self.db.execute("SELECT data FROM quality_gates WHERE signature=?", (signature,)).fetchone()
        return bool(r and json.loads(r[0]).get("passed"))

    @contextmanager
    def exclusive(self):
        with (self.root / "data/collector.lock").open("a") as f:
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise RuntimeError("Another collection is already running in this project") from None
            try:
                yield
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
