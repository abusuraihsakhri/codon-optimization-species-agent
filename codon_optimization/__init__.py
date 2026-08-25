"""
Codon Optimization: CAI, GC content, RSCU, hairpin detection, rare codons.
"""
__version__ = "3.0.0"

from .engine import (
    CAIResult, GCContentResult, HairpinResult, RareCodon, OptimizationResult,
    CODON_TABLE, AA_TO_CODONS,
    ORGANISM_RSCU, E_COLI_RSCU, HUMAN_RSCU, YEAST_RSCU,
    calculate_rscu, calculate_cai,
    calculate_gc_content, gc_content_by_position, analyze_gc_content,
    detect_hairpins,
    identify_rare_codons,
    select_optimal_codon, optimize_codons,
    full_optimization,
)
