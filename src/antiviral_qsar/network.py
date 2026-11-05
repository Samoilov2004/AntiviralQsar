from __future__ import annotations

import json
import hashlib
import time
from urllib.parse import urlparse
import httpx

from .config import Config, digest
from .storage import now


class Network:
    def __init__(self, config: Config):
        self.config = config
        self.root = config.root
        self.client = httpx.Client(timeout=config.settings["network"]["timeout_seconds"],
                                   follow_redirects=True, headers={"User-Agent": "AntiviralQSAR/0.1 research"})
        self.last_request: dict[str, float] = {}
        self.requests = 0
        self.cache_hits = 0
        self.snapshots: dict[str, str] = {}

    def get(self, url, params=None, *, max_age=None, max_bytes=100_000_000):
        if urlparse(url).scheme != "https":
            raise ValueError("Only HTTPS source URLs are supported")
        request_url = str(httpx.URL(url, params=params)) if params else url
        key = digest(request_url)
        blob = self.root / "data/raw" / f"{key}.bin"
        meta = blob.with_suffix(".json")
        if blob.exists() and meta.exists() and (max_age is None or time.time() - meta.stat().st_mtime < max_age):
            self.cache_hits += 1
            metadata = json.loads(meta.read_text())
            snapshot = metadata.get("snapshot", str(blob.relative_to(self.root)))
            self.snapshots[request_url] = snapshot
            return (self.root / snapshot).read_bytes(), snapshot
        host = urlparse(url).netloc
        retries = int(self.config.settings["network"]["retries"])
        error = None
        for attempt in range(retries):
            delay = max(0, 1 / self.config.settings["network"]["requests_per_second"] -
                        (time.monotonic() - self.last_request.get(host, 0)))
            time.sleep(delay)
            self.last_request[host] = time.monotonic()
            try:
                self.requests += 1
                with self.client.stream("GET", url, params=params) as response:
                    if response.status_code in (429, 500, 502, 503, 504):
                        pause = response.headers.get("Retry-After", "")
                        time.sleep(min(float(pause) if pause.isdigit() else 2 ** attempt, 30))
                    response.raise_for_status()
                    data = bytearray()
                    for block in response.iter_bytes():
                        data.extend(block)
                        if len(data) > max_bytes:
                            raise ValueError(f"Source exceeds {max_bytes} byte limit: {url}")
                    content = bytes(data)
                    # Do not overwrite evidence from older runs when a URL's contents change.
                    content_key = hashlib.sha256(content).hexdigest()
                    version = self.root / "data/raw" / f"{key}-{content_key[:16]}.bin"
                    version.write_bytes(content)
                    blob.write_bytes(content)
                    metadata = {"url": request_url, "final_url": str(response.url), "retrieved_at": now(),
                                "content_type": response.headers.get("Content-Type", ""),
                                "sha256": content_key, "snapshot": str(version.relative_to(self.root))}
                    meta.write_text(json.dumps(metadata, indent=2))
                    version.with_suffix(".json").write_text(json.dumps(metadata, indent=2))
                    self.snapshots[request_url] = str(version.relative_to(self.root))
                    return content, str(version.relative_to(self.root))
            except (httpx.HTTPError, OSError) as exc:
                error = exc
                if isinstance(exc, httpx.HTTPStatusError) and exc.response.status_code not in (429,500,502,503,504):
                    break
                if attempt + 1 < retries:
                    time.sleep(2 ** attempt)
        status = f":HTTP_{error.response.status_code}" if isinstance(error, httpx.HTTPStatusError) else ""
        raise RuntimeError(f"Source request failed: {request_url}: {type(error).__name__}{status}") from error

    def json(self, url, params=None, **kwargs):
        content, path = self.get(url, params, **kwargs)
        return json.loads(content), path

    def close(self):
        self.client.close()
