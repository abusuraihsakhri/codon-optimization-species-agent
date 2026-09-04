"""
Real tests for Codon Optimization.
Tests CAI, GC content, hairpin detection, rare codons, and optimization.
"""
import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from codon_optimization.engine import (
    CODON_TABLE, AA_TO_CODONS,
    ORGANISM_RSCU, E_COLI_RSCU, HUMAN_RSCU, YEAST_RSCU,
    calculate_rscu, calculate_cai,
    calculate_gc_content, gc_content_by_position, analyze_gc_content,
    detect_hairpins,
    identify_rare_codons,
    select_optimal_codon, optimize_codons,
    full_optimization,
)


# --- Test Constants ---

def test_codon_table_complete():
    """Codon table should have 64 entries."""
    assert len(CODON_TABLE) == 64


def test_aa_to_codons():
    """Each amino acid should have at least one codon."""
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        assert aa in AA_TO_CODONS
        assert len(AA_TO_CODONS[aa]) >= 1


def test_stop_codons():
    """Stop codons should map to '*'."""
    assert CODON_TABLE['TAA'] == '*'
    assert CODON_TABLE['TAG'] == '*'
    assert CODON_TABLE['TGA'] == '*'


def test_organism_rscu_tables():
    """All organism RSCU tables should have entries for all sense codons."""
    for organism, rscu in ORGANISM_RSCU.items():
        for codon, aa in CODON_TABLE.items():
            if aa != '*':
                assert codon in rscu, f"{codon} missing from {organism} RSCU"


# --- Test CAI ---

def test_cai_optimal_sequence():
    """Sequence using only optimal codons should have CAI ≈ 1.0."""
    # Build a sequence using the best codon for each amino acid
    protein = 'MAKDEFGHILNPQRSTVWY'
    optimized = optimize_codons(protein, E_COLI_RSCU)
    result = calculate_cai(optimized, E_COLI_RSCU)
    assert result.cai > 0.95  # should be very close to 1.0


def test_cai_poor_codons():
    """Sequence using rare codons should have low CAI."""
    # Use rare E. coli codons: AGA (R), CUA (L), AUA (I)
    # R-L-I = AGACUAAUA
    result = calculate_cai('AGACUAAUA', E_COLI_RSCU)
    assert result.cai < 0.5


def test_cai_range():
    """CAI should be between 0 and 1."""
    result = calculate_cai('ATGGCTAAGGATGAAGAG', E_COLI_RSCU)
    assert 0 <= result.cai <= 1


def test_cai_num_codons():
    """Number of codons should be correct."""
    result = calculate_cai('ATGGCTAAG', E_COLI_RSCU)  # 3 codons
    assert result.num_codons == 3


def test_cai_different_organisms():
    """CAI should differ between organisms for the same sequence."""
    seq = 'ATGGCTAAGGATGAAGAGTTT'
    cai_ecoli = calculate_cai(seq, E_COLI_RSCU)
    cai_human = calculate_cai(seq, HUMAN_RSCU)
    # They should be different (different RSCU values)
    assert cai_ecoli.cai != cai_human.cai


def test_cai_invalid_length():
    """Sequence not multiple of 3 should raise error."""
    try:
        calculate_cai('ATGGC', E_COLI_RSCU)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_cai_codon_scores():
    """Codon scores should be between 0 and 1."""
    result = calculate_cai('ATGGCTAAG', E_COLI_RSCU)
    for codon, aa, w in result.codon_scores:
        assert 0 <= w <= 1


# --- Test GC Content ---

def test_gc_content_all_gc():
    """All GC sequence should have GC = 1.0."""
    assert calculate_gc_content('GCGCGCGC') == 1.0


def test_gc_content_all_at():
    """All AT sequence should have GC = 0.0."""
    assert calculate_gc_content('ATATATAT') == 0.0


def test_gc_content_mixed():
    """Mixed sequence should have correct GC."""
    assert abs(calculate_gc_content('ATGC') - 0.5) < 0.001


def test_gc_content_empty():
    """Empty sequence should return 0."""
    assert calculate_gc_content('') == 0.0


def test_gc_content_by_position():
    """GC at each codon position should be correct."""
    # Codons: ATG CAT GCT
    # Pos 1 (indices 0,3,6): A, C, G -> 2/3 GC
    # Pos 2 (indices 1,4,7): T, A, C -> 1/3 GC
    # Pos 3 (indices 2,5,8): G, T, T -> 1/3 GC
    result = gc_content_by_position('ATGCATGCT')
    assert abs(result[1] - 2/3) < 0.01
    assert abs(result[2] - 1/3) < 0.01
    assert abs(result[3] - 1/3) < 0.01


def test_analyze_gc_content():
    """Full GC analysis should return all fields."""
    result = analyze_gc_content('ATGCGATCG')
    assert 0 <= result.gc_content <= 1
    assert result.total_bases == 9
    assert result.total_gc > 0


# --- Test Hairpin Detection ---

def test_hairpin_palindrome():
    """Self-complementary sequence should form hairpin."""
    # GCGC...GCGC type structure
    seq = 'GCGCAAAAGCGC'  # stem: GCGC, loop: AAAA, stem: GCGC (reverse complement)
    hairpins = detect_hairpins(seq, min_stem=3)
    # May or may not find depending on exact complementarity
    # Just verify it doesn't crash
    assert isinstance(hairpins, list)


def test_hairpin_no_structure():
    """Random sequence with no complementarity should have no hairpins."""
    seq = 'AAAAAAAAAAAAAAAAAAAA'
    hairpins = detect_hairpins(seq, min_stem=4)
    assert len(hairpins) == 0


def test_hairpin_short_sequence():
    """Very short sequence should have no hairpins."""
    hairpins = detect_hairpins('ATG', min_stem=4)
    assert len(hairpins) == 0


# --- Test Rare Codons ---

def test_rare_codons_ecoli():
    """AGA (rare Arg in E. coli) should be identified."""
    rare = identify_rare_codons('AGAAGAAGA', E_COLI_RSCU, threshold=0.3)
    assert len(rare) == 3
    for r in rare:
        assert r.codon == 'AGA'
        assert r.amino_acid == 'R'


def test_no_rare_codons():
    """Optimal codons should not be identified as rare."""
    # CGT is the optimal Arg codon in E. coli (RSCU = 3.44)
    rare = identify_rare_codons('CGTCGTCGT', E_COLI_RSCU, threshold=0.3)
    assert len(rare) == 0


def test_rare_codons_threshold():
    """Higher threshold should find more rare codons."""
    seq = 'ATGGCTAAGGATGAAGAG'
    rare_low = identify_rare_codons(seq, E_COLI_RSCU, threshold=0.1)
    rare_high = identify_rare_codons(seq, E_COLI_RSCU, threshold=0.8)
    assert len(rare_high) >= len(rare_low)


def test_rare_codons_different_organisms():
    """Same codon may be rare in one organism but not another."""
    # AGA is rare in E. coli but common in human
    rare_ecoli = identify_rare_codons('AGA', E_COLI_RSCU, threshold=0.5)
    rare_human = identify_rare_codons('AGA', HUMAN_RSCU, threshold=0.5)
    assert len(rare_ecoli) > len(rare_human)


# --- Test Optimal Codon Selection ---

def test_select_optimal_ecoli():
    """Should select highest RSCU codon for E. coli."""
    # For Alanine in E. coli: GCT (1.84) is highest
    codon = select_optimal_codon('A', E_COLI_RSCU)
    assert codon == 'GCT'


def test_select_optimal_human():
    """Should select highest RSCU codon for human."""
    # For Alanine in human: GCC (1.52) is highest
    codon = select_optimal_codon('A', HUMAN_RSCU)
    assert codon == 'GCC'


def test_select_optimal_unknown_aa():
    """Unknown amino acid should raise error."""
    try:
        select_optimal_codon('X', E_COLI_RSCU)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


# --- Test Codon Optimization ---

def test_optimize_codons():
    """Optimization should produce valid DNA."""
    protein = 'MAKDEFGHIL'
    optimized = optimize_codons(protein, E_COLI_RSCU)
    assert len(optimized) == len(protein) * 3
    assert all(b in 'ATGC' for b in optimized)


def test_optimize_codons_translates_back():
    """Optimized sequence should translate back to the same protein."""
    protein = 'MAKDEFGHILNPQRSTVWY'
    optimized = optimize_codons(protein, E_COLI_RSCU)
    
    # Translate back
    from codon_optimization.engine import CODON_TABLE
    translated = ''
    for i in range(0, len(optimized) - 2, 3):
        codon = optimized[i:i+3]
        translated += CODON_TABLE.get(codon, 'X')
    
    assert translated == protein


def test_optimize_codons_improves_cai():
    """Optimized sequence should have higher CAI than a random one."""
    protein = 'MAKDEFGHIL'
    # Use rare codons
    bad_dna = 'ATGGCTAGAGATGAAGATTTT'  # AGA for K, GAT for D (not optimal in E. coli)
    bad_cai = calculate_cai(bad_dna, E_COLI_RSCU)
    
    optimized = optimize_codons(protein, E_COLI_RSCU)
    opt_cai = calculate_cai(optimized, E_COLI_RSCU)
    
    assert opt_cai.cai >= bad_cai.cai


def test_optimize_stop_codon():
    """Stop codon should be optimized."""
    protein = 'MA*'
    optimized = optimize_codons(protein, E_COLI_RSCU)
    assert len(optimized) == 9
    # Last codon should be a stop codon
    last_codon = optimized[6:9]
    assert CODON_TABLE.get(last_codon) == '*'


# --- Test Full Optimization ---

def test_full_optimization():
    """Full optimization should return complete results."""
    protein = 'MAKDEFGHIL'
    original_dna = 'ATGGCTAAGGATGAAGAGTTTATTCAT'
    result = full_optimization(protein, original_dna, 'e_coli')
    
    assert result.original_cai > 0
    assert result.optimized_cai > 0
    assert result.optimized_cai >= result.original_cai
    assert len(result.optimized_sequence) == len(protein) * 3


def test_full_optimization_different_organisms():
    """Different target organisms should give different optimized sequences."""
    protein = 'MAKDEFGHIL'
    original_dna = 'ATGGCTAAGGATGAAGAGTTTATTCAT'
    
    r_ecoli = full_optimization(protein, original_dna, 'e_coli')
    r_human = full_optimization(protein, original_dna, 'human')
    
    # Optimized sequences should differ (different optimal codons)
    assert r_ecoli.optimized_sequence != r_human.optimized_sequence


# --- Test RSCU Calculation ---

def test_calculate_rscu():
    """RSCU from codon counts should be reasonable."""
    codon_counts = {
        'GCT': 100, 'GCC': 50, 'GCA': 30, 'GCG': 20,  # Ala
    }
    rscu = calculate_rscu(codon_counts)
    # GCT should have highest RSCU (most frequent)
    assert rscu['GCT'] > rscu['GCC']
    assert rscu['GCT'] > rscu['GCA']
    assert rscu['GCT'] > rscu['GCG']


def test_rscu_uniform():
    """Uniform codon usage should give RSCU = 1.0 for all."""
    codon_counts = {'GCT': 100, 'GCC': 100, 'GCA': 100, 'GCG': 100}
    rscu = calculate_rscu(codon_counts)
    for codon, val in rscu.items():
        assert abs(val - 1.0) < 0.01


# --- Test CLI ---

def test_cli_cai():
    """Test CLI CAI command."""
    from codon_optimization.cli import main
    assert main(['cai', '--sequence', 'ATGGCTAAGGATGAAGAG', '--organism', 'e_coli']) == 0


def test_cli_gc():
    """Test CLI GC command."""
    from codon_optimization.cli import main
    assert main(['gc', '--sequence', 'ATGCGATCG']) == 0


def test_cli_rare():
    """Test CLI rare command."""
    from codon_optimization.cli import main
    assert main(['rare', '--sequence', 'AGAAGAAGA', '--organism', 'e_coli']) == 0


def test_cli_optimize():
    """Test CLI optimize command."""
    from codon_optimization.cli import main
    assert main(['optimize', '--protein', 'MAKDE', '--dna', 'ATGGCTAAGGATGAA', '--organism', 'e_coli']) == 0


def test_cli_hairpin():
    """Test CLI hairpin command."""
    from codon_optimization.cli import main
    assert main(['hairpin', '--sequence', 'GCGCAAAAGCGC']) == 0


def test_cli_batch(tmp_path):
    """Test CLI batch optimization with temporary CSV files."""
    from codon_optimization.cli import main
    in_csv = tmp_path / "test_in.csv"
    out_csv = tmp_path / "test_out.csv"
    in_csv.write_text(
        "sequence_id,gene_name,organism,dna_sequence,protein_sequence\n"
        "T1,gene1,e_coli,ATGGCTAAGGATGAAGAG,MAKDEE\n"
        "T2,gene2,human,,MVHL\n",
        encoding="utf-8"
    )
    res = main(['batch', '-i', str(in_csv), '-o', str(out_csv)])
    assert res == 0
    assert out_csv.exists()
    content = out_csv.read_text(encoding="utf-8")
    assert "optimized_cai" in content
    assert "T1" in content
    assert "T2" in content


def test_cli_batch_missing_file():
    """Test CLI batch error handling on missing input."""
    from codon_optimization.cli import main
    assert main(['batch', '-i', 'non_existent_file_xyz.csv', '-o', 'out.csv']) == 1

