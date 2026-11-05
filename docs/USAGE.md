# AntiviralQSAR usage guide

A local CLI for collecting small-molecule antiviral measurements from ChEMBL, PubChem BioAssay, and Europe PMC. It turns a virus name into evidence-backed CSV files for QSAR research. Qwen runs locally through MLX; QSAR model training is outside this project's scope.

**Experimental extraction:** literature records remain in the review table until the extractor passes an independent audit. The 95% precision target has not been validated. See [the pilot report](../PILOT.md).

## Quick start

Requires macOS on Apple Silicon, Python 3.11, and an internet connection for source APIs and the initial model download.

```bash
git clone https://github.com/Samoilov2004/AntiviralQsar.git
cd AntiviralQsar
./setup.sh
./run.sh collect --virus "Vaccinia virus"
```

Setup creates a local `.venv`, installs the pinned dependencies, and downloads the pinned Qwen model. If Python is not on your PATH, set `AQ_PYTHON=/path/to/python3.11` when running setup.

After discovery, the CLI shows the candidate counts and asks how many papers to process: a number, `all`, or `0`. The budget counts candidate documents, which may contain no suitable measurements. `0` processes databases only. Database measurements are unlimited by default.

After collecting database results, the agent compares article DOI and PMID with the references of the collected ChEMBL and PubChem assays. DOI URLs, case and PMID prefixes are normalized; an exact match on either identifier is sufficient. Articles without a match are processed first. Within each group, the existing relevance ranking is used. `--papers 10` is a total budget, not ten papers per group.

Matched articles are processed only when imported records lack a structure, cell line, viral strain, duration, method or usable concentration, or when a compound lacks CC50 or an antiviral endpoint. IC50 and EC50 are not required simultaneously. An absent endpoint is a search opportunity, not a claim that the paper reports it; CC50 is still not paired with an antiviral assay automatically. An assay with no imported measurements or incomplete retrieval also leaves a gap. Articles with no detected gaps are skipped without spending the paper budget or calling Qwen.

Publication matching covers this run's database assays, including their companion cytotoxicity assays. It does not prove that an unmatched article is absent everywhere in either database. Missing publication identifiers, incomplete reference lookups and multi-publication assays are reported explicitly; uncertain coverage does not cause an article to be skipped. Interrupted runs retain their saved selection and order.

```bash
./run.sh collect --virus "Vaccinia virus" --papers 10
./run.sh collect --virus "Vaccinia virus" --papers 0
./run.sh status
./run.sh resume RUN_ID
./run.sh export RUN_ID
```

A compact terminal progress display replaces per-experiment output:

```text
Assays ━━━━━╸────────────────  21%  669/3194  2,525 left  ETA 0:42:10
```

The numbers above are illustrative. ETA is estimated from the current phase's observed processing speed; it becomes available after initial completions and can vary with assay or article size. Saved work is included in the count when resuming. Papers display a second bar for their fragments. Detailed processing events go to `logs/RUN_ID.log`; failures and coverage are also saved in the run report. Redirected output contains phase summaries instead of animation or thousands of per-item lines.

Ctrl+C saves progress. `resume RUN_ID` continues the saved paper selection without prompting again. `resume RUN_ID --papers N` explicitly selects a new budget from the remaining papers. Completed jobs are not processed again. Ambiguous virus names require a selection or an explicit `--tax-id`; noninteractive collection requires `--papers`.

For a small technical trial, use `--assay-limit 20` to select at most 20 primary assays per database. This explicit restriction is recorded in the report. `--no-model --papers 0` skips inference and sends uncertain assay contexts to review.

## Where the model lives

GitHub contains **the code, settings, tests, dependency lock, and `model.lock.json`**. Model weights, environments, caches, collected articles, and results are local files excluded by `.gitignore`.

`model.lock.json` identifies the Hugging Face repository and a full commit SHA. On another machine, `./setup.sh` or `./run.sh model download` fetches that exact version into `models/qwen3-8b-4bit/`. The installed revision is recorded in `models/manifest.json`. An existing installation with a different model revision is not silently overwritten; use a separate checkout when changing model versions.

The default model is [mlx-community/Qwen3-8B-4bit](https://huggingface.co/mlx-community/Qwen3-8B-4bit). Its weights retain their upstream license; the repository's code license does not replace the model's license.

```bash
./run.sh model download
./run.sh model start
./run.sh model stop
./run.sh doctor --network
```

The model server starts automatically when needed, listens on `127.0.0.1:8080`, and handles one extraction at a time. Articles are split into fragments. Source texts are processed locally. The server remains available after collection; `model stop` releases its memory. Model logs are in `logs/model.log`. Adjust timeouts, port, and token budgets in `config.toml`.

All runtime files are contained in the checkout. When moving an existing installation, stop the collector first and recreate `.venv` with `./setup.sh`; Python virtual environments should not be copied between paths. Existing `models/`, `data/`, and caches can be retained locally without committing them.

## Output and provenance

Each run writes to `data/runs/RUN_ID/`:

| File | Contents |
|---|---|
| `measurements.csv` | Measurements passing automatic checks; this is not manual validation |
| `review.csv` | Unresolved records and explicit `review_reasons` |
| `all_observations.csv` | Every observation, including identified copies and provenance |
| `report.json` | Coverage, errors, source queries, sampling limits, timing, and model usage |
| `manifest.json` | Settings, model and rule versions, and source snapshot references |
| `code/` | Code snapshots for reproducible execution |

One row represents one measurement. Fields include original and standardized SMILES, virus/strain, cell line, endpoint, original value/units, censoring relation, normalized concentration, conditions, source, and supporting evidence.

IC50, EC50, and CC50 remain separate. Different experiments are not averaged, and missing values are not guessed. Concentrations are normalized to `value_nm`; `p_activity` is computed only for exact values. Inequalities are preserved. Unsupported units, including mass concentrations without a verified conversion, remain under review. CC50 is not automatically paired with antiviral activity to infer a selectivity index.

`source_snapshot` points to an immutable cached response; `evidence_locator` identifies a record, table row, or page. Raw documents and API responses live in `data/raw/`, state in `data/state.sqlite`, and inference traces in `data/llm/`. Identical values from different experiments are retained. Cross-source copies are merged only with a matching original assay identity, structure, value, and conditions.

For matched papers, extraction adds candidate measurements or missing context without overwriting database values. `supplement_of` links proposed supplements to database observation IDs; ambiguous links and conflicting values remain under review. `literature_role=already_covered` means the extracted fields add no information to an imported record. These rows remain in `all_observations.csv` and raw snapshots but are omitted from `measurements.csv` and `review.csv`; this is a coverage decision, not proof of identical experiments. Cached full extractions are reassessed against the current database results. `report.json` includes per-paper matches, gaps, skip reasons, reference errors and supplementary-row counts under `paper_coverage`.

## Independent audit

All LLM-extracted measurements remain in `review.csv` until an independent sample of at least 100 records reaches 95% complete-record precision for that extractor version. Assay contexts determined only by the model have a separate `assay_classifier_requires_review` flag.

```bash
./run.sh audit prepare RUN_ID --size 100 --method llm --seed 42
# Independently compare every sampled record with the original source.
# Fill in correct=1/0, reviewer, and note in the generated CSV.
./run.sh audit evaluate data/runs/RUN_ID/audit-llm-42.csv
./run.sh export RUN_ID
```

Check the entire relationship: structure, experiment, virus/cell context, endpoint, value, units, and inequality. An incorrect or unsupported component makes the record incorrect. Preparing an empty annotation template does not pass the audit. Fewer than 100 available records does not lower the threshold. Audit results are tied to the recorded extractor version.

Recall requires complete independent annotation of selected documents. Supply a CSV with `document_id,expected,matched,reviewer` using `audit evaluate ... --recall-file recall.csv`. Without it, recall is reported as `null`. Tests verify software behavior, not scientific extraction precision.

## Development

```bash
python3.11 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check src tests
```

The tests run without a model download or live source APIs; CI uses this CPU-only path. Local inference additionally requires the `mlx` extra and Apple Silicon. `requirements.lock.txt` records the tested macOS environment used by `setup.sh`.

## Coverage and limitations

The sources are ChEMBL, PubChem BioAssay, and Europe PMC. Candidate counts are not counts of usable QSAR observations. Full texts are retrieved from open XML, with a fallback to explicitly declared open PDF links. Text PDFs and CSV/TSV, XLSX, TXT, and XML attachments are supported. Scans, chemical drawings, old XLS files, and DOCX are not parsed; access restrictions are not bypassed.

Compound names are resolved only through explicit, unambiguous matches. Document-local labels such as `7b` require an associated structure record or remain under review. Open-access metadata does not guarantee that a full text can actually be downloaded. Failed sources, unavailable documents, partial extraction, and unprocessed papers are reported explicitly.

The code is distributed under the existing [GNU GPL v3 license](../LICENSE). Source datasets, publications, and model weights retain their respective licenses.
