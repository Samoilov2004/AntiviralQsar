# AntiviralQSAR

Collect SMILES, cell lines and IC50/EC50/CC50 measurements for QSAR from ChEMBL, PubChem BioAssay and Europe PMC. Qwen reads papers locally; every result keeps a link to its source.

Requires an Apple Silicon Mac and Python 3.11. From the project folder, install once and collect data for **VACV (Vaccinia virus)**:

```bash
./setup.sh
./run.sh collect --virus "Vaccinia virus" --papers 10
```

This processes database results and up to 10 papers: new papers first; papers matched to database records by DOI/PMID only to fill gaps. Use `--papers 0` for databases only, or omit it to choose interactively.

Results go to `data/runs/<run_id>/`: `measurements.csv`, `review.csv` and `report.json`. Paper extractions still need review before use in a model.

Setup downloads Qwen automatically. Weights, caches and results stay local and are excluded from Git.

[More commands and audit details](docs/USAGE.md) · [License](LICENSE)
