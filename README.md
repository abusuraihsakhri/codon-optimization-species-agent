# Codon Optimization Species Agent

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![Pytest](https://img.shields.io/badge/Pytest-Passing-brightgreen.svg?logo=pytest&logoColor=white)
![CI](https://img.shields.io/badge/CI-Active-brightgreen.svg)

A high-performance computational biology toolkit for **Codon Adaptation Index (CAI)** calculation, **GC content optimization**, **species-specific Relative Synonymous Codon Usage (RSCU)** adaptation, **mRNA hairpin structure detection**, and **rare codon minimization** across synthetic biology and recombinant biotherapeutics hosts (*Escherichia coli*, *Homo sapiens*, and *Saccharomyces cerevisiae*).

---

## 🔬 Scientific & Clinical Fundamentals

Heterologous protein expression often encounters rate-limiting translation bottlenecks caused by:
1. **Codon bias discordance**: Host organism tRNA pools vary widely in isoacceptor abundance. Depleted cognate tRNAs stall ribosomal elongation and promote ribosome collisions or mistranslation.
2. **Pathological GC composition**: Incompatible GC ratios trigger transcriptional termination, premature polyadenylation, or poor mRNA stability.
3. **Unfavorable secondary mRNA folding**: Stable secondary hairpin loops ($\Delta G < -5.0 \text{ kcal/mol}$) near the 5' untranslated region (5' UTR) or ribosomal binding site (RBS) severely impair ribosomal initiation and yield.

### 1. Codon Adaptation Index (CAI) Formulation

Codon Adaptation Index (Sharp & Li, 1987) measures the degree of codon preference bias towards highly expressed genes:

$$w_i = \frac{\text{RSCU}_i}{\max_{j \in \text{syn}(i)} \text{RSCU}_j}$$

$$\text{CAI} = \exp \left( \frac{1}{L} \sum_{k=1}^L \ln(w_k) \right) = \left( \prod_{k=1}^L w_k \right)^{1/L}$$

Where:
- $\text{RSCU}_i$: Relative Synonymous Codon Usage of codon $i$.
- $\max_{j \in \text{syn}(i)} \text{RSCU}_j$: Maximum RSCU value among all synonymous codons coding for the same amino acid.
- $w_i$: Relative adaptiveness value of codon $i$ ($0 < w_i \le 1.0$).
- $L$: Total number of sense codons in the coding sequence (excluding start/stop codons).

### 2. Relative Synonymous Codon Usage (RSCU) Formulation

$$\text{RSCU}_{ij} = \frac{X_{ij}}{\frac{1}{n_i} \sum_{j=1}^{n_i} X_{ij}}$$

Where $X_{ij}$ is the observed frequency of the $j$-th codon for the $i$-th amino acid, and $n_i$ represents the degree of synonymous degeneracy ($1 \le n_i \le 6$).

---

## 📊 Species-Specific Codon Preference & tRNA Abundance

Optimal synonymous codon choices for key amino acids across reference expression hosts:

| Amino Acid | Degeneracy | *E. coli* Preferred (RSCU) | Human / HEK293 Preferred (RSCU) | *S. cerevisiae* Preferred (RSCU) |
|:---|:---:|:---|:---|:---|
| **Leu (L)** | 6 | `CTG` (5.25) | `CTG` (2.18) | `TTA` (1.33) / `TTG` (1.27) |
| **Arg (R)** | 6 | `CGT` (3.44) | `AGA` (1.15) / `AGG` (1.10) | `AGA` (2.67) |
| **Gly (G)** | 4 | `GGT` (2.86) | `GGC` (1.35) | `GGT` (1.60) |
| **Ala (A)** | 4 | `GCT` (1.84) | `GCC` (1.52) | `GCT` (1.87) |
| **Val (V)** | 4 | `GTT` (1.53) | `GTG` (1.35) | `GTT` (1.40) / `GTG` (1.40) |
| **Thr (T)** | 4 | `ACC` (2.50) | `ACC` (1.40) | `ACT` (1.67) |
| **Pro (P)** | 4 | `CCG` (2.83) | `CCC` (1.23) | `CCA` (1.87) |
| **Ile (I)** | 3 | `ATC` (1.82) | `ATC` (1.33) | `ATT` (1.27) |
| **Glu (E)** | 2 | `GAA` (1.58) | `GAG` (1.15) | `GAA` (1.60) |
| **Asp (D)** | 2 | `GAT` (1.14) | `GAT` (1.08) | `GAT` (1.20) |
| **Lys (K)** | 2 | `AAA` (1.65) | `AAG` (1.13) | `AAA` (1.47) |
| **Asn (N)** | 2 | `AAC` (1.63) | `AAC` (1.08) | `AAT` (1.07) |
| **Gln (Q)** | 2 | `CAG` (1.66) | `CAG` (1.48) | `CAA` (1.67) |
| **His (H)** | 2 | `CAC` (1.19) | `CAC` (1.17) | `CAC` (1.13) |
| **Tyr (Y)** | 2 | `TAC` (1.28) | `TAC` (1.17) | `TAC` (1.07) |
| **Cys (C)** | 2 | `TGC` (1.54) | `TGC` (1.28) | `TGC` (1.33) |
| **Phe (F)** | 2 | `TTC` (1.26) | `TTT` (1.07) | `TTT` (1.13) |
| **Ser (S)** | 6 | `TCT` (1.30) | `AGC` (1.28) / `TCC` (1.22) | `TCT` (1.87) |
| **Met (M)** | 1 | `ATG` (1.00) | `ATG` (1.00) | `ATG` (1.00) |
| **Trp (W)** | 1 | `TGG` (1.00) | `TGG` (1.00) | `TGG` (1.00) |

---

## 🚀 Quickstart & CLI Usage

The command-line interface provides both granular analysis subcommands and high-throughput batch execution.

### 1. Batch Processing via CSV (Primary Pipeline)

Process multi-sequence libraries with automatic CAI calculation, GC optimization, and host adaptation:

```bash
# Run batch optimization with short flags
python cli.py batch -i sample.csv -o out_smoke.csv

# Or use standard long flags and specify target organism
python cli.py batch --input sample.csv --output out_smoke.csv --organism e_coli
```

#### Expected Input Format (`sample.csv`)

| Column Name | Description | Example |
|:---|:---|:---|
| `sequence_id` | Unique sequence identifier | `SEQ_001` |
| `gene_name` | Gene or therapeutic target | `GFP_reporter` |
| `organism` | Host organism (`e_coli`, `human`, `yeast`) | `e_coli` |
| `dna_sequence` | Input coding DNA (length multiple of 3) | `ATGGCTAGCAAAGGAGAA...` |
| `protein_sequence` | Target translated amino acid sequence | `MASKGEELFTGVVPILV...` |

### 2. Codon Adaptation Index (CAI)

```bash
# Calculate CAI for an E. coli sequence
python cli.py cai --sequence ATGGCTAAGGATGAAGAG --organism e_coli

# Include per-codon w_i relative adaptiveness in JSON
python cli.py cai --sequence ATGGCTAAGGATGAAGAG --organism e_coli --json --verbose
```

### 3. GC Content & Positional Bias Analysis

Analyze overall GC% and positional bias ($GC_1$, $GC_2$, $GC_3$):

```bash
python cli.py gc --sequence ATGCGATCGGCT --json
```

### 4. Rare Codon Detection

Identify inhibitory rare codons below an RSCU threshold (default: $0.3$):

```bash
python cli.py rare --sequence AGAAGGAGGCTA --organism e_coli --threshold 0.30
```

### 5. Single Sequence Codon Optimization

Optimize an amino acid sequence directly against host tRNA abundance:

```bash
python cli.py optimize --protein MAKDE --dna ATGGCTAAGGATGAA --organism e_coli
```

### 6. mRNA Hairpin Secondary Structure Scan

Detect potential secondary hairpin structures with thermodynamic energy estimation:

```bash
python cli.py hairpin --sequence GCGCAAAAGCGC --min-stem 4 --max-loop 10
```

---

## 🧪 Testing & Verification

Run the comprehensive unit and integration test suite:

```bash
python -m pytest -p no:zarr -v
```

Execute the batch CLI smoke test:

```bash
python cli.py batch -i sample.csv -o out_smoke.csv
python -c "import csv; assert len(list(csv.DictReader(open('out_smoke.csv')))) >= 4; print('Batch smoke test passed!')"
rm out_smoke.csv
```

---

## 📦 Project Architecture

```
codon-optimization-species-agent/
├── .github/workflows/ci.yml   # Multi-version CI matrix (Python 3.10, 3.11, 3.12)
├── codon_optimization/
│   ├── __init__.py            # Package exports
│   ├── cli.py                 # CLI implementation (cai, gc, rare, optimize, hairpin, batch)
│   └── engine.py              # RSCU tables, CAI formulas, GC analyzers, hairpin prediction
├── tests/
│   ├── __init__.py            # Test package root
│   ├── test_codon_optimization.py  # CAI, GC, rare, hairpin, CLI & batch unit tests
│   └── test_enrichment.py     # Advanced optimization engine tests
├── cli.py                     # Root CLI entrypoint
├── codon_optimization_app.py  # Standalone launcher
├── enrichment.py              # Extended multi-species design & translation kinetics suite
├── sample.csv                 # Realistic biological test dataset (GFP, Insulin, Nanobody, Pol)
├── pyproject.toml             # Packaging configuration
└── README.md                  # Domain documentation & mathematical specifications
```

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
