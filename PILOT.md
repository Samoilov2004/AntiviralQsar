# Initial pilot — November 5, 2026

The development pilot used Dengue virus (NCBI Taxonomy 12637), 20 candidate assays from each database, and 10 candidate papers. This was a technical sample, not a complete collection. Production collection has no implicit assay cap.

| Outcome | Count |
|---|---:|
| Observations | 188 |
| Passed automatic checks | 14 |
| Require review | 174 |
| Structured database observations | 65 |
| LLM-extracted observations | 123 |
| Papers fully processed | 5 |
| Papers partially processed | 1 |
| Full texts unavailable | 4 |

The run finished with `partial` status. Structured results were revalidated without changing their identifiers, original numerical values, units, or structures. Previous exports were retained locally. These metrics describe the initial development pilot, not a fresh benchmark of subsequent code changes.

**The 95% precision target has not been independently validated.** All 123 literature observations remained in review; none passed every structural, contextual, and evidence check in that pilot. This is an automatic eligibility result, not a measured precision score. A reproducible sample of 100 records was prepared for independent annotation; recall was not measured.

Typical issues included mixed-virus tables, local compound labels, image-only structures, unresolved names, and evidence quotes that did not match the source exactly. Truncated responses caused partial extraction. Later changes split such fragments immediately instead of repeating the same output budget.

Processing after discovery took 6,866.73 seconds (114.4 minutes), including 100 model requests. Inference accounted for 6,423.41 seconds, with 276,454 input and 148,141 output tokens. The model was Qwen3-8B-4bit on an Apple M3 Pro with 36 GB unified memory. The server's recorded peak physical footprint was 17.7 GB across its lifetime, including additional checks; sampled RSS did not capture all Metal memory. These are observations for this workload, not guaranteed runtime or memory limits.

Raw sources, result files, annotations, and model weights are intentionally excluded from Git. Run the collector locally to produce a dataset and inspect its report before using it for QSAR.
