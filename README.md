# Codon Optimization Species Agent

### [Open the Live Application →](https://abusuraihsakhri.github.io/codon-optimization-species-agent/)

A dependency-free Python toolkit for codon-usage analysis and deterministic synonymous-codon optimization using bundled reference RSCU tables for *Escherichia coli*, human, and *Saccharomyces cerevisiae*.

The same Python analysis engine is used by the command-line interface and the browser application.

## Features

- Codon Adaptation Index (CAI) calculation from a supplied RSCU reference table.
- Overall GC content plus GC1, GC2, and GC3 positional content.
- Rare-codon screening using a configurable RSCU threshold.
- Deterministic synonymous-codon selection using the highest bundled RSCU value for each amino acid.
- Coding DNA/RNA validation, including frame alignment and protein/CDS translation consistency checks.
- A simple inverted-repeat hairpin screen with a heuristic pseudo-energy score.
- CSV batch processing with explicit per-row validation status and error messages.
- Browser execution through Pyodide; no application backend is required.

## Scientific scope

CAI is a codon-usage bias metric based on relative adaptiveness to a reference set, originally described by Sharp and Li (1987). It should not be interpreted as a direct measurement or guarantee of protein expression.

The RSCU tables in this repository are simplified bundled reference tables. They are not cell-line-specific datasets and should not be presented as HEK293, CHO, tissue-specific, strain-specific, or current experimental codon-usage measurements without separate validation.

The hairpin routine is a sequence heuristic based on complementary stems. Its returned score is **not** a thermodynamic folding free energy and does not replace tools such as RNAfold or experimental structure measurements.

## Browser application

GitHub Pages deployment is configured from the repository root. The page loads a pinned Pyodide runtime and executes `codon_optimization/engine.py` directly in the browser.

Entered sequences stay in the browser; the application has no sequence-processing backend. Initial page use requires network access to download Pyodide from jsDelivr.

## CLI usage

The project requires Python 3.10 or newer and has no runtime Python dependencies.

```bash
python cli.py cai --sequence ATGGCTAAGGATGAAGAG --organism e_coli
python cli.py gc --sequence ATGCGATCGGCT --json
python cli.py rare --sequence AGAAGAAGA --organism e_coli --threshold 0.30
python cli.py optimize --protein MAKDE --dna ATGGCTAAGGATGAA --organism e_coli
python cli.py hairpin --sequence GCGCAAAAGCGC --min-stem 4 --max-loop 10
```

RNA input containing `U` is normalized to DNA `T`. Ambiguous nucleotide symbols are rejected rather than silently omitted.

### Batch CSV

```bash
python cli.py batch -i sample.csv -o results.csv
```

Expected columns:

| Column | Required | Description |
| --- | --- | --- |
| `sequence_id` | No | Row identifier; generated when absent |
| `gene_name` | No | Descriptive label |
| `organism` | No | `e_coli`, `human`, or `yeast`; CLI default is used when absent |
| `dna_sequence` | One of DNA/protein | Codon-aligned coding DNA or RNA |
| `protein_sequence` | One of DNA/protein | Standard one-letter amino-acid sequence |

Output includes `status` and `message` columns so malformed rows are preserved and reported instead of being silently skipped.

## Local development

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install -e ".[dev]"
python -m pytest -p no:zarr -q
```

Additional checks used by CI:

```bash
python -m compileall -q codon_optimization cli.py codon_optimization_app.py tests
node --check assets/app.js
python cli.py batch -i sample.csv -o out_smoke.csv
```

## Project structure

```text
.
├── index.html
├── assets/
│   ├── app.js
│   └── styles.css
├── codon_optimization/
│   ├── __init__.py
│   ├── cli.py
│   └── engine.py
├── tests/
│   └── test_codon_optimization.py
├── .github/workflows/
│   ├── ci.yml
│   └── pages.yml
├── cli.py
├── codon_optimization_app.py
├── sample.csv
└── pyproject.toml
```

## Browser compatibility

The web application requires a modern browser with WebAssembly support. Python execution is provided by Pyodide 314.0.7. The layout supports light and dark themes, keyboard focus states, desktop view, and a compact mobile input/results switcher.

## References

- Sharp PM, Li WH. The codon adaptation index—a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*. 1987;15(3):1281-1295. doi:10.1093/nar/15.3.1281.
- Pyodide documentation: browser-based Python runtime and deployment guidance.

## License

MIT. See [LICENSE](LICENSE).
