from __future__ import annotations

import csv
import json
import random
import shutil
import sys
from pathlib import Path
from collections import Counter

import typer
from rich.console import Console
from rich.table import Table

from .config import Config
from .model import LocalModel
from .network import Network
from .pipeline import Pipeline, export_run
from .schema import Virus
from .sources import taxonomy_candidates
from .storage import Store, now

app = typer.Typer(no_args_is_help=True, help="Collect evidence-backed antiviral measurements for QSAR")
model_app = typer.Typer(help="Local Qwen / MLX model")
audit_app = typer.Typer(help="Independent accuracy audit; an unlabelled sample is not a completed audit")
app.add_typer(model_app, name="model")
app.add_typer(audit_app, name="audit")
console = Console()


def components():
    cfg = Config()
    return cfg, Store(cfg.root)


def validate_budget(value):
    if value is None:
        return None
    if value.strip().lower() in ("all", "\u0432\u0441\u0435", "\u0432\u0441\u0451"):
        return "all"
    try:
        amount = int(value)
        if amount < 0:
            raise ValueError()
        return str(amount)
    except ValueError:
        raise typer.BadParameter("Enter a non-negative integer or all") from None


def show_inventory(inventory):
    table = Table(title="Discovery results")
    table.add_column("Source")
    table.add_column("Availability")
    table.add_column("Candidates")
    for source, entry in inventory.get("sources", {}).items():
        count = entry.get("assays", entry.get("publications", "—"))
        table.add_row(source, entry.get("status", "unknown"), str(count))
    console.print(table)
    literature = inventory.get("sources", {}).get("europepmc", {})
    console.print(f"Open full texts: {literature.get('open_fulltext', 0)}; "
                  f"previously processed: {literature.get('previously_processed', 0)}")
    console.print("These are candidates; relevance and measurement suitability are checked during processing.", style="dim")


def choose_budget(store, rid, requested):
    if requested is not None:
        return validate_budget(requested)
    eligible = [j for j in store.jobs(rid, "paper", unfinished=True)
                if j["payload"].get("pmcid") and j["payload"].get("isOpenAccess") == "Y"]
    console.print(f"Available papers remaining in this run: {len(eligible)}")
    if not eligible:
        return "0"
    if not sys.stdin.isatty():
        raise typer.BadParameter("No interactive terminal. Add --papers N or --papers all; the run has been saved.")
    while True:
        try:
            return validate_budget(typer.prompt("How many new papers should be processed? [number / all / 0]"))
        except typer.BadParameter as exc:
            console.print(str(exc), style="yellow")


@app.command()
def collect(virus: str = typer.Option(..., help="Virus name"),
            papers: str | None = typer.Option(None, help="N/all; prompts after discovery when omitted"),
            tax_id: str | None = typer.Option(None, help="Explicit NCBI Taxonomy ID"),
            assay_limit: int | None = typer.Option(None, min=1, help="Explicit assay limit per source for a pilot; unlimited by default"),
            no_model: bool = typer.Option(False, "--no-model", help="Use explicit assay descriptions only; uncertain cases go to review")):
    budget = validate_budget(papers)
    if no_model and budget not in (None, "0"):
        raise typer.BadParameter("--no-model only supports --papers 0")
    cfg, store = components()
    pipeline = Pipeline(cfg, store, use_model=not no_model)
    rid = None
    try:
        with store.exclusive():
            candidates = taxonomy_candidates(pipeline.net, virus, tax_id)
            if len(candidates) > 1:
                for i, candidate in enumerate(candidates, 1):
                    console.print(f"{i}. {candidate.name} (tax_id={candidate.tax_id})")
                if not sys.stdin.isatty():
                    raise typer.BadParameter("Ambiguous name. Retry with a --tax-id from the list.")
                selected = typer.prompt("Choose a number", type=int)
                if not 1 <= selected <= len(candidates):
                    raise typer.BadParameter("Invalid selection")
                resolved = candidates[selected-1]
            else:
                resolved = candidates[0]
            console.print(f"Resolved virus: {resolved.name} (NCBI:{resolved.tax_id})")
            rid = store.new_run(resolved.model_dump(), {"papers":budget, "assay_limit":assay_limit, "no_model":no_model})
            console.print(f"Run: {rid}")
            inventory = pipeline.discover(rid, resolved, assay_limit=assay_limit)
            show_inventory(inventory)
            budget = "0" if no_model else choose_budget(store, rid, budget)
            options = store.run(rid)["options"]
            options["papers"] = budget
            store.update_run(rid, options=options)
            folder, counts = pipeline.execute(rid, budget)
            console.print(f"Status: {store.run(rid)['status']}; output: {folder}\n{counts}")
    except KeyboardInterrupt:
        if rid:
            store.update_run(rid, status="interrupted")
            export_run(cfg, store, rid)
            console.print(f"Progress saved. Resume with: ./run.sh resume {rid}")
        raise typer.Exit(130) from None
    except typer.BadParameter:
        raise
    except Exception as exc:
        if rid:
            store.update_run(rid, status="failed")
            store.event(rid, str(exc), "error")
            export_run(cfg, store, rid)
        console.print(f"[red]{type(exc).__name__}: {exc}[/red]")
        raise typer.Exit(1) from None
    finally:
        pipeline.net.close()


@app.command()
def resume(run_id: str, papers: str | None = typer.Option(None, help="New budget for remaining papers: N/all"),
           retry_discovery: bool = typer.Option(True, help="Retry discovery for previously unavailable sources")):
    cfg, store = components()
    run = store.run(run_id)
    pipeline = Pipeline(cfg, store, use_model=not run["options"].get("no_model"))
    try:
        with store.exclusive():
            # Existing successful jobs remain untouched. Running/failed jobs are retried.
            if retry_discovery:
                pipeline.discover(run_id, Virus.model_validate(run["virus"]), assay_limit=run["options"].get("assay_limit"))
            show_inventory(store.run(run_id)["inventory"])
            if papers is None and 'selected_paper_keys' in run['options']:
                budget = None
            else:
                budget = choose_budget(store, run_id, papers if papers is not None else run['options'].get('papers'))
            if run["options"].get("no_model") and budget not in (None, "0"):
                raise typer.BadParameter("This run uses --no-model; resume with --papers 0")
            folder, counts = pipeline.execute(run_id, budget)
            console.print(f"Status: {store.run(run_id)['status']}; output: {folder}\n{counts}")
    except KeyboardInterrupt:
        store.update_run(run_id, status="interrupted")
        export_run(cfg, store, run_id)
        console.print(f"Progress saved: {run_id}")
        raise typer.Exit(130) from None
    except typer.BadParameter:
        raise
    except Exception as exc:
        store.update_run(run_id, status="failed")
        store.event(run_id, str(exc), "error")
        export_run(cfg, store, run_id)
        console.print(f"[red]{type(exc).__name__}: {exc}[/red]")
        raise typer.Exit(1) from None
    finally:
        pipeline.net.close()


@app.command()
def status(run_id: str | None = typer.Argument(None)):
    _, store = components()
    if run_id:
        run = store.run(run_id)
        console.print({"run_id":run_id, "status":run["status"], "virus":run["virus"],
                       "jobs":dict(Counter(j["kind"]+":"+j["status"] for j in store.jobs(run_id)))})
    else:
        for row in store.db.execute("SELECT id,status,virus FROM runs ORDER BY created DESC LIMIT 30"):
            console.print(row["id"], row["status"], json.loads(row["virus"])["name"])


@app.command()
def export(run_id: str):
    cfg, store = components()
    console.print(export_run(cfg, store, run_id))


@app.command()
def doctor(network: bool = typer.Option(False, help="Also check the public source APIs")):
    cfg, store = components()
    checks = {"python":sys.version.split()[0], "root":str(cfg.root),
              "free_disk_gb":round(shutil.disk_usage(cfg.root).free/1e9,1),
              "model_files":(cfg.model_dir/"config.json").exists(), "model_revision":cfg.model_revision,
              "model_server":LocalModel(cfg,store).health()}
    for name in ("rdkit", "mlx.core", "pymupdf", "openpyxl"):
        try:
            module = __import__(name, fromlist=["_"])
            checks[name] = "ok"
            if name == "mlx.core":
                checks["metal_available"] = module.metal.is_available()
        except Exception as exc:
            checks[name] = type(exc).__name__
    if network:
        net = Network(cfg)
        endpoints = {"chembl":"https://www.ebi.ac.uk/chembl/api/data/status.json",
                     "pubchem":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/SMILES/JSON",
                     "europepmc":"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=dengue&format=json&pageSize=1"}
        for source,url in endpoints.items():
            try:
                net.json(url, max_age=0)
                checks[source] = "available"
            except Exception as exc:
                checks[source] = str(exc)
        net.close()
    console.print_json(data=checks)


@model_app.command("download")
def model_download():
    cfg, store = components()
    console.print("Downloading model to", cfg.model_dir)
    console.print("Revision:", LocalModel(cfg, store).download())


@model_app.command("start")
def model_start():
    cfg, store = components()
    LocalModel(cfg, store).start()
    console.print("Model ready:", cfg.model_url)


@model_app.command("stop")
def model_stop():
    cfg, store = components()
    console.print("Stopped" if LocalModel(cfg, store).stop() else "No managed server found")


@audit_app.command("prepare")
def audit_prepare(run_id: str, size: int = typer.Option(100, min=1), seed: int = 42,
                  method: str = typer.Option("llm", help="llm / structured / all")):
    cfg, store = components()
    rows = store.observations(run_id)
    if method not in ("llm", "structured", "all"):
        raise typer.BadParameter("method must be llm, structured, or all")
    rows = [m for m in rows if method == "all" or m.extraction_method == method]
    random.Random(seed).shuffle(rows)
    chosen = rows[:size]
    signatures = {m.extractor_signature for m in chosen if m.extractor_signature}
    if len(signatures) > 1:
        raise typer.BadParameter("The sample contains multiple extractor versions; use a run with one version")
    folder = cfg.root / "data/runs" / run_id
    path = folder / f"audit-{method}-{seed}.csv"
    if path.exists():
        raise typer.BadParameter(f"File already exists: {path}; refusing to overwrite annotations")
    fields = ["observation_id", "correct", "reviewer", "note", "source_url", "source_snapshot",
              "evidence_locator", "compound", "smiles", "endpoint", "value_original", "units_original",
              "relation", "virus", "cell_line", "duration", "evidence", "context_evidence", "extraction_method"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for m in chosen:
            writer.writerow({**m.model_dump(), "correct":"", "reviewer":"", "note":""})
    manifest = {"run_id":run_id, "seed":seed, "method":method, "requested":size,
                "available":len(rows), "sample_ids":[m.observation_id for m in chosen],
                "signature":next(iter(signatures), cfg.extractor_signature), "created":now()}
    path.with_suffix(".json").write_text(json.dumps(manifest, indent=2))
    console.print(f"Audit sample: {path}; {len(chosen)} records. Manually label correct as 1/0.")


@audit_app.command("evaluate")
def audit_evaluate(path: Path, recall_file: Path | None = typer.Option(None, help="CSV: document_id,expected,matched,reviewer")):
    cfg, store = components()
    manifest = json.loads(path.with_suffix(".json").read_text())
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    ids = [r["observation_id"] for r in rows]
    if len(set(ids)) != len(ids) or set(ids) != set(manifest["sample_ids"]):
        raise typer.BadParameter("The sample was changed or contains duplicates; create a new audit")
    if not rows or any(r.get("correct") not in ("0", "1") or not r.get("reviewer", "").strip() for r in rows):
        raise typer.BadParameter("Independently review every row: correct=0/1 and reviewer are required")
    precision = sum(int(r["correct"]) for r in rows) / len(rows)
    methods = Counter(r["extraction_method"] for r in rows)
    gate = (manifest["method"] == "llm" and set(methods) == {"llm"} and
            len(rows) >= cfg.settings["quality"]["minimum_audit_records"] and
            precision >= cfg.settings["quality"]["minimum_precision"])
    report = {"sample_size":len(rows), "precision":precision, "methods":dict(methods),
              "passed":gate, "recall":None, "signature":manifest["signature"],
              "reviewers":sorted({r["reviewer"] for r in rows}), "audit_file":str(path), "evaluated_at":now()}
    if recall_file:
        with recall_file.open(encoding="utf-8-sig") as f:
            recalls = list(csv.DictReader(f))
        if not recalls or any(not r.get("reviewer") or int(r["expected"]) < int(r["matched"]) or int(r["matched"]) < 0 for r in recalls):
            raise typer.BadParameter("Invalid independent recall annotations")
        expected = sum(int(r["expected"]) for r in recalls)
        report["recall"] = sum(int(r["matched"]) for r in recalls) / expected if expected else None
        report["recall_documents"] = len(recalls)
    else:
        report["recall_note"] = "Not measured; requires complete independent annotation of selected documents"
    output = path.with_name(path.stem+"-result.json")
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2))
    if manifest["method"] == "llm":
        store.db.execute("INSERT OR REPLACE INTO quality_gates VALUES(?,?,?)",
                         (manifest["signature"], json.dumps(report), now()))
        store.db.commit()
    console.print_json(data=report)
