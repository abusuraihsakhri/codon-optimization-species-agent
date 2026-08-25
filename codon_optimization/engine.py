"""
Codon Optimization Engine
Codon Adaptation Index (CAI), GC content optimization, RSCU tables,
mRNA hairpin detection, rare codon identification, and optimal codon selection.
"""
import math
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field


# --- Standard Genetic Code ---
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Build reverse mapping: amino acid -> list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)


# --- RSCU Tables for Common Organisms ---
# RSCU = Relative Synonymous Codon Usage
# Values represent observed frequency / expected frequency (if uniform)

# E. coli RSCU (simplified from Kazusa codon usage database)
E_COLI_RSCU = {
    'TTT': 0.74, 'TTC': 1.26, 'TTA': 0.41, 'TTG': 0.15,
    'CTT': 0.33, 'CTC': 0.14, 'CTA': 0.07, 'CTG': 5.25,
    'ATT': 1.03, 'ATC': 1.82, 'ATA': 0.07, 'ATG': 1.00,
    'GTT': 1.53, 'GTC': 0.57, 'GTA': 0.49, 'GTG': 1.37,
    'TCT': 1.30, 'TCC': 0.85, 'TCA': 0.35, 'TCG': 0.25,
    'CCT': 0.47, 'CCC': 0.12, 'CCA': 0.58, 'CCG': 2.83,
    'ACT': 1.54, 'ACC': 2.50, 'ACA': 0.20, 'ACG': 0.72,
    'GCT': 1.84, 'GCC': 1.44, 'GCA': 0.58, 'GCG': 1.14,
    'TAT': 0.72, 'TAC': 1.28, 'TAA': 1.86, 'TAG': 0.10,
    'CAT': 0.81, 'CAC': 1.19, 'CAA': 0.34, 'CAG': 1.66,
    'AAT': 0.37, 'AAC': 1.63, 'AAA': 1.65, 'AAG': 0.35,
    'GAT': 1.14, 'GAC': 0.86, 'GAA': 1.58, 'GAG': 0.42,
    'TGT': 0.46, 'TGC': 1.54, 'TGA': 0.04, 'TGG': 1.00,
    'CGT': 3.44, 'CGC': 1.94, 'CGA': 0.10, 'CGG': 0.12,
    'AGT': 0.35, 'AGC': 1.15, 'AGA': 0.07, 'AGG': 0.04,
    'GGT': 2.86, 'GGC': 2.42, 'GGA': 0.13, 'GGG': 0.59,
}

# Human RSCU (simplified)
HUMAN_RSCU = {
    'TTT': 1.07, 'TTC': 0.93, 'TTA': 0.45, 'TTG': 0.62,
    'CTT': 0.65, 'CTC': 0.80, 'CTA': 0.55, 'CTG': 2.18,
    'ATT': 0.97, 'ATC': 1.33, 'ATA': 0.70, 'ATG': 1.00,
    'GTT': 0.55, 'GTC': 0.72, 'GTA': 0.38, 'GTG': 1.35,
    'TCT': 0.98, 'TCC': 1.22, 'TCA': 0.72, 'TCG': 0.38,
    'CCT': 0.88, 'CCC': 1.23, 'CCA': 0.93, 'CCG': 0.96,
    'ACT': 0.90, 'ACC': 1.40, 'ACA': 0.92, 'ACG': 0.78,
    'GCT': 0.98, 'GCC': 1.52, 'GCA': 0.82, 'GCG': 0.68,
    'TAT': 0.83, 'TAC': 1.17, 'TAA': 0.65, 'TAG': 0.35,
    'CAT': 0.83, 'CAC': 1.17, 'CAA': 0.52, 'CAG': 1.48,
    'AAT': 0.92, 'AAC': 1.08, 'AAA': 0.87, 'AAG': 1.13,
    'GAT': 1.08, 'GAC': 0.92, 'GAA': 0.85, 'GAG': 1.15,
    'TGT': 0.72, 'TGC': 1.28, 'TGA': 0.50, 'TGG': 1.00,
    'CGT': 0.35, 'CGC': 0.92, 'CGA': 0.48, 'CGG': 1.05,
    'AGT': 0.72, 'AGC': 1.28, 'AGA': 1.15, 'AGG': 1.10,
    'GGT': 0.50, 'GGC': 1.35, 'GGA': 0.92, 'GGG': 1.23,
}

# Yeast (S. cerevisiae) RSCU (simplified)
YEAST_RSCU = {
    'TTT': 1.13, 'TTC': 0.87, 'TTA': 1.33, 'TTG': 1.27,
    'CTT': 0.53, 'CTC': 0.20, 'CTA': 0.87, 'CTG': 0.80,
    'ATT': 1.27, 'ATC': 0.93, 'ATA': 0.80, 'ATG': 1.00,
    'GTT': 1.40, 'GTC': 0.67, 'GTA': 0.53, 'GTG': 1.40,
    'TCT': 1.87, 'TCC': 1.07, 'TCA': 1.07, 'TCG': 0.20,
    'CCT': 1.27, 'CCC': 0.53, 'CCA': 1.87, 'CCG': 0.33,
    'ACT': 1.67, 'ACC': 1.00, 'ACA': 1.00, 'ACG': 0.33,
    'GCT': 1.87, 'GCC': 0.87, 'GCA': 0.93, 'GCG': 0.33,
    'TAT': 0.93, 'TAC': 1.07, 'TAA': 1.00, 'TAG': 0.50,
    'CAT': 0.87, 'CAC': 1.13, 'CAA': 1.67, 'CAG': 0.33,
    'AAT': 1.07, 'AAC': 0.93, 'AAA': 1.47, 'AAG': 0.53,
    'GAT': 1.20, 'GAC': 0.80, 'GAA': 1.60, 'GAG': 0.40,
    'TGT': 0.67, 'TGC': 1.33, 'TGA': 0.50, 'TGG': 1.00,
    'CGT': 0.40, 'CGC': 0.13, 'CGA': 0.13, 'CGG': 0.07,
    'AGT': 1.00, 'AGC': 1.00, 'AGA': 2.67, 'AGG': 1.33,
    'GGT': 1.60, 'GGC': 0.60, 'GGA': 1.20, 'GGG': 0.60,
}

ORGANISM_RSCU = {
    'e_coli': E_COLI_RSCU,
    'human': HUMAN_RSCU,
    'yeast': YEAST_RSCU,
}


# --- Data Models ---

@dataclass
class CAIResult:
    """Codon Adaptation Index result."""
    cai: float
    log_cai: float
    codon_scores: List[Tuple[str, str, float]]  # (codon, aa, w_i)
    num_codons: int


@dataclass
class GCContentResult:
    """GC content analysis."""
    gc_content: float
    gc_at_positions: Dict[int, float]  # GC% at codon positions 1, 2, 3
    total_gc: int
    total_bases: int


@dataclass
class HairpinResult:
    """mRNA hairpin detection result."""
    has_hairpin: bool
    position: int
    stem_length: int
    free_energy: float  # simplified
    sequence: str


@dataclass
class RareCodon:
    """A rare codon in the sequence."""
    position: int
    codon: str
    amino_acid: str
    rscu: float
    frequency: float


@dataclass
class OptimizationResult:
    """Complete codon optimization result."""
    original_sequence: str
    optimized_sequence: str
    original_cai: float
    optimized_cai: float
    original_gc: float
    optimized_gc: float
    num_changes: int
    changes: List[Tuple[int, str, str]]  # (position, original, optimized)


# --- CAI Calculation ---

def calculate_rscu(codon_counts: Dict[str, int]) -> Dict[str, float]:
    """Calculate Relative Synonymous Codon Usage from codon counts.
    
    RSCU_ij = X_ij / (1/n_i × Σ X_ij)
    where X_ij is the count of codon j for amino acid i, n_i is the number
    of synonymous codons for amino acid i.
    """
    # Group codons by amino acid
    aa_codon_counts = {}
    for codon, count in codon_counts.items():
        aa = CODON_TABLE.get(codon)
        if aa and aa != '*':
            if aa not in aa_codon_counts:
                aa_codon_counts[aa] = {}
            aa_codon_counts[aa][codon] = count
    
    rscu = {}
    for aa, codons in aa_codon_counts.items():
        total = sum(codons.values())
        n = len(AA_TO_CODONS.get(aa, []))
        if total == 0 or n == 0:
            continue
        for codon in codons:
            rscu[codon] = codons[codon] / (total / n)
    
    return rscu


def calculate_cai(dna_sequence: str, rscu_table: Dict[str, float]) -> CAIResult:
    """Calculate Codon Adaptation Index.
    
    CAI = (Π w_i)^(1/L)
    where w_i = RSCU_i / max(RSCU for that amino acid)
    L = number of codons (excluding start/stop)
    
    Args:
        dna_sequence: DNA coding sequence (must be multiple of 3)
        rscu_table: RSCU values for the target organism
    Returns:
        CAIResult with CAI value and per-codon scores
    """
    sequence = dna_sequence.upper().replace(' ', '').replace('\n', '')
    
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) must be a multiple of 3")
    
    # Calculate max RSCU per amino acid
    max_rscu = {}
    for codon, rscu_val in rscu_table.items():
        aa = CODON_TABLE.get(codon)
        if aa and aa != '*':
            if aa not in max_rscu or rscu_val > max_rscu[aa]:
                max_rscu[aa] = rscu_val
    
    codon_scores = []
    log_sum = 0.0
    count = 0
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        
        if aa == '*' or aa == 'X':
            codon_scores.append((codon, aa, 0.0))
            continue
        
        rscu_val = rscu_table.get(codon, 0.0)
        max_val = max_rscu.get(aa, 1.0)
        
        if max_val > 0:
            w = rscu_val / max_val
        else:
            w = 0.0
        
        # Clamp to avoid log(0)
        w = max(w, 1e-10)
        
        codon_scores.append((codon, aa, round(w, 4)))
        log_sum += math.log(w)
        count += 1
    
    if count == 0:
        return CAIResult(cai=0.0, log_cai=0.0, codon_scores=codon_scores, num_codons=0)
    
    cai = math.exp(log_sum / count)
    
    return CAIResult(
        cai=round(cai, 6),
        log_cai=round(log_sum / count, 6),
        codon_scores=codon_scores,
        num_codons=count,
    )


# --- GC Content ---

def calculate_gc_content(dna_sequence: str) -> float:
    """Calculate overall GC content as a fraction."""
    seq = dna_sequence.upper().replace(' ', '').replace('\n', '')
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for b in seq if b in 'GC')
    return gc / len(seq)


def gc_content_by_position(dna_sequence: str) -> Dict[int, float]:
    """Calculate GC content at each codon position (1, 2, 3)."""
    seq = dna_sequence.upper().replace(' ', '').replace('\n', '')
    counts = {1: [0, 0], 2: [0, 0], 3: [0, 0]}
    
    for i, base in enumerate(seq):
        pos = (i % 3) + 1
        counts[pos][1] += 1  # total
        if base in 'GC':
            counts[pos][0] += 1  # GC count
    
    return {
        pos: counts[pos][0] / counts[pos][1] if counts[pos][1] > 0 else 0.0
        for pos in [1, 2, 3]
    }


def analyze_gc_content(dna_sequence: str) -> GCContentResult:
    """Full GC content analysis."""
    seq = dna_sequence.upper().replace(' ', '').replace('\n', '')
    gc = sum(1 for b in seq if b in 'GC')
    return GCContentResult(
        gc_content=gc / len(seq) if seq else 0.0,
        gc_at_positions=gc_content_by_position(seq),
        total_gc=gc,
        total_bases=len(seq),
    )


# --- mRNA Hairpin Detection ---

def detect_hairpins(dna_sequence: str, min_stem: int = 4,
                     max_loop: int = 10) -> List[HairpinResult]:
    """Detect potential hairpin structures in mRNA.
    
    A hairpin consists of:
    - A stem (complementary base pairing)
    - A loop (unpaired region)
    
    Uses simple Watson-Crick complementarity (A-U/T, G-C).
    Minimum stem length: 4 bp.
    """
    seq = dna_sequence.upper().replace(' ', '').replace('\n', '')
    # Convert T to U for RNA
    rna = seq.replace('T', 'U')
    
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    hairpins = []
    
    for i in range(len(rna)):
        for loop_size in range(3, min(max_loop + 1, (len(rna) - i) // 2)):
            for stem_len in range(min_stem, min(20, (len(rna) - i - loop_size) // 2 + 1)):
                # Check if we have enough sequence
                left_start = i
                left_end = i + stem_len
                right_start = i + stem_len + loop_size
                right_end = right_start + stem_len
                
                if right_end > len(rna):
                    continue
                
                # Check complementarity
                matches = 0
                for j in range(stem_len):
                    left_base = rna[left_start + j]
                    right_base = rna[right_end - 1 - j]
                    if complement.get(left_base) == right_base:
                        matches += 1
                
                # Require at least 70% complementarity
                if matches >= stem_len * 0.7:
                    # Simplified free energy: -1.5 kcal/mol per GC pair, -1.0 per AU pair
                    energy = 0.0
                    for j in range(stem_len):
                        left_base = rna[left_start + j]
                        right_base = rna[right_end - 1 - j]
                        if complement.get(left_base) == right_base:
                            if left_base in 'GC':
                                energy -= 1.5
                            else:
                                energy -= 1.0
                    
                    hairpins.append(HairpinResult(
                        has_hairpin=True,
                        position=i,
                        stem_length=stem_len,
                        free_energy=round(energy, 2),
                        sequence=rna[left_start:right_end],
                    ))
    
    # Return only the strongest hairpin (most negative energy)
    if hairpins:
        best = min(hairpins, key=lambda h: h.free_energy)
        return [best]
    return []


# --- Rare Codon Identification ---

def identify_rare_codons(dna_sequence: str, rscu_table: Dict[str, float],
                          threshold: float = 0.3) -> List[RareCodon]:
    """Identify rare codons in a sequence.
    
    A codon is "rare" if its RSCU value is below the threshold.
    
    Args:
        dna_sequence: DNA coding sequence
        rscu_table: RSCU values for the target organism
        threshold: RSCU threshold for "rare" (default 0.3)
    Returns:
        List of RareCodon objects
    """
    seq = dna_sequence.upper().replace(' ', '').replace('\n', '')
    rare = []
    
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        
        if aa == '*' or aa == 'X':
            continue
        
        rscu = rscu_table.get(codon, 0.0)
        
        if rscu < threshold:
            # Calculate frequency from RSCU
            synonymous = AA_TO_CODONS.get(aa, [])
            total_rscu = sum(rscu_table.get(c, 0) for c in synonymous)
            freq = rscu / total_rscu if total_rscu > 0 else 0
            
            rare.append(RareCodon(
                position=i // 3,
                codon=codon,
                amino_acid=aa,
                rscu=round(rscu, 4),
                frequency=round(freq, 4),
            ))
    
    return rare


# --- Optimal Codon Selection ---

def select_optimal_codon(amino_acid: str, rscu_table: Dict[str, float]) -> str:
    """Select the optimal codon for an amino acid based on RSCU.
    
    Returns the codon with the highest RSCU value.
    """
    codons = AA_TO_CODONS.get(amino_acid, [])
    if not codons:
        raise ValueError(f"Unknown amino acid: {amino_acid}")
    
    best_codon = codons[0]
    best_rscu = rscu_table.get(codons[0], 0.0)
    
    for codon in codons[1:]:
        rscu = rscu_table.get(codon, 0.0)
        if rscu > best_rscu:
            best_rscu = rscu
            best_codon = codon
    
    return best_codon


def optimize_codons(protein_sequence: str, rscu_table: Dict[str, float],
                     gc_target: Optional[float] = None) -> str:
    """Optimize a protein sequence's codons for a target organism.
    
    For each amino acid, selects the codon with the highest RSCU value.
    Optionally adjusts for GC content target.
    
    Args:
        protein_sequence: Amino acid sequence
        rscu_table: RSCU values for target organism
        gc_target: Target GC content (0-1), or None for pure CAI optimization
    Returns:
        Optimized DNA sequence
    """
    optimized = []
    
    for aa in protein_sequence.upper():
        if aa == '*':
            # Stop codon - pick the most common one
            stop_codons = ['TAA', 'TAG', 'TGA']
            best_stop = max(stop_codons, key=lambda c: rscu_table.get(c, 0))
            optimized.append(best_stop)
        elif aa in AA_TO_CODONS:
            codon = select_optimal_codon(aa, rscu_table)
            optimized.append(codon)
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    
    return ''.join(optimized)


# --- Full Optimization ---

def full_optimization(protein_sequence: str, original_dna: str,
                       target_organism: str = 'e_coli') -> OptimizationResult:
    """Run complete codon optimization analysis.
    
    1. Calculate CAI of original sequence
    2. Optimize codons for target organism
    3. Calculate CAI of optimized sequence
    4. Compare GC content
    """
    rscu = ORGANISM_RSCU.get(target_organism, E_COLI_RSCU)
    
    original_cai = calculate_cai(original_dna, rscu)
    optimized_dna = optimize_codons(protein_sequence, rscu)
    optimized_cai = calculate_cai(optimized_dna, rscu)
    
    original_gc = calculate_gc_content(original_dna)
    optimized_gc = calculate_gc_content(optimized_dna)
    
    # Track changes
    changes = []
    orig_upper = original_dna.upper().replace(' ', '').replace('\n', '')
    for i in range(0, min(len(orig_upper), len(optimized_dna)) - 2, 3):
        orig_codon = orig_upper[i:i+3]
        opt_codon = optimized_dna[i:i+3]
        if orig_codon != opt_codon:
            changes.append((i // 3, orig_codon, opt_codon))
    
    return OptimizationResult(
        original_sequence=original_dna,
        optimized_sequence=optimized_dna,
        original_cai=original_cai.cai,
        optimized_cai=optimized_cai.cai,
        original_gc=original_gc,
        optimized_gc=optimized_gc,
        num_changes=len(changes),
        changes=changes,
    )
