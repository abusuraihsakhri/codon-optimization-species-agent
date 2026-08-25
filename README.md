# Codon Optimization

Codon Adaptation Index (CAI), GC content analysis, RSCU tables for common organisms, mRNA hairpin detection, rare codon identification, and optimal codon selection.

## Features

- **CAI Calculation**: CAI = (Π w_i)^(1/L) where w_i = RSCU_i / max(RSCU)
- **RSCU Tables**: E. coli, Human, Yeast (S. cerevisiae)
- **GC Content**: Overall and per-codon-position analysis
- **Hairpin Detection**: mRNA secondary structure avoidance
- **Rare Codon Identification**: Based on RSCU threshold
- **Codon Optimization**: Select optimal codons for target organism

## Quick Start

```bash
# Calculate CAI
python cli.py cai --sequence ATGGCTAAGGATGAAGAG --organism e_coli

# GC content analysis
python cli.py gc --sequence ATGCGATCGATCG

# Detect hairpins
python cli.py hairpin --sequence GCGCAAAAGCGC

# Identify rare codons
python cli.py rare --sequence AGAAGAAGA --organism e_coli --threshold 0.3

# Optimize codons
python cli.py optimize --protein MAKDEFGHIL --dna ATGGCTAAGGATGAAGAGTTTATTCAT --organism e_coli
```

## Python API

```python
from codon_optimization import (
    calculate_cai, calculate_gc_content, detect_hairpins,
    identify_rare_codons, optimize_codons, full_optimization,
    E_COLI_RSCU, HUMAN_RSCU, YEAST_RSCU,
)

# CAI
result = calculate_cai('ATGGCTAAGGATGAAGAG', E_COLI_RSCU)
print(f"CAI: {result.cai:.4f}")

# Optimize for E. coli
optimized = optimize_codons('MAKDEFGHIL', E_COLI_RSCU)
print(f"Optimized: {optimized}")

# Full optimization with comparison
result = full_optimization('MAKDEFGHIL', 'ATGGCTAAGGATGAAGAG', 'e_coli')
print(f"CAI: {result.original_cai:.4f} -> {result.optimized_cai:.4f}")
```

## Supported Organisms

| Organism | Key Features |
|----------|-------------|
| E. coli | CTG (Leu), CGT (Arg), GGT (Gly) preferred |
| Human | CTG (Leu), GCC (Ala), AGA/AGG (Arg) preferred |
| Yeast | TTA (Leu), AGA (Arg), GGT (Gly) preferred |

## CAI Interpretation

| CAI Value | Interpretation |
|-----------|---------------|
| 1.0 | All optimal codons |
| > 0.8 | Highly adapted |
| 0.5-0.8 | Moderately adapted |
| < 0.5 | Poorly adapted |

## License

MIT License.
