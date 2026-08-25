"""
CLI module for Codon Optimization.
"""
import argparse
import json
import sys

from .engine import (
    ORGANISM_RSCU,
    calculate_cai, calculate_gc_content, gc_content_by_position, analyze_gc_content,
    detect_hairpins, identify_rare_codons,
    select_optimal_codon, optimize_codons, full_optimization,
)


def cmd_cai(args):
    """Calculate Codon Adaptation Index."""
    rscu = ORGANISM_RSCU.get(args.organism, ORGANISM_RSCU['e_coli'])
    result = calculate_cai(args.sequence, rscu)
    
    output = {
        'cai': result.cai,
        'log_cai': result.log_cai,
        'num_codons': result.num_codons,
        'organism': args.organism,
    }
    
    if args.json:
        output['codon_scores'] = [
            {'codon': c, 'aa': aa, 'w': w}
            for c, aa, w in result.codon_scores
        ]
        print(json.dumps(output, indent=2))
    else:
        print(f"CAI Analysis ({args.organism})")
        print("=" * 40)
        print(f"  CAI:      {result.cai:.6f}")
        print(f"  ln(CAI):  {result.log_cai:.6f}")
        print(f"  Codons:   {result.num_codons}")
        if args.verbose:
            print(f"\n  {'Codon':>6} {'AA':>3} {'w_i':>8}")
            print("  " + "-" * 20)
            for c, aa, w in result.codon_scores:
                print(f"  {c:>6} {aa:>3} {w:>8.4f}")
    
    return 0


def cmd_gc(args):
    """Analyze GC content."""
    result = analyze_gc_content(args.sequence)
    
    output = {
        'gc_content': round(result.gc_content, 4),
        'gc_position_1': round(result.gc_at_positions.get(1, 0), 4),
        'gc_position_2': round(result.gc_at_positions.get(2, 0), 4),
        'gc_position_3': round(result.gc_at_positions.get(3, 0), 4),
        'total_gc': result.total_gc,
        'total_bases': result.total_bases,
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("GC Content Analysis")
        print("=" * 40)
        print(f"  Overall:    {result.gc_content:.4f} ({result.gc_content*100:.1f}%)")
        print(f"  Position 1: {result.gc_at_positions.get(1, 0):.4f}")
        print(f"  Position 2: {result.gc_at_positions.get(2, 0):.4f}")
        print(f"  Position 3: {result.gc_at_positions.get(3, 0):.4f}")
        print(f"  GC bases:   {result.total_gc}/{result.total_bases}")
    
    return 0


def cmd_hairpin(args):
    """Detect mRNA hairpins."""
    hairpins = detect_hairpins(args.sequence, args.min_stem, args.max_loop)
    
    if args.json:
        output = [
            {
                'position': h.position,
                'stem_length': h.stem_length,
                'free_energy': h.free_energy,
                'sequence': h.sequence,
            }
            for h in hairpins
        ]
        print(json.dumps(output, indent=2))
    else:
        if hairpins:
            print(f"Found {len(hairpins)} hairpin(s):")
            for h in hairpins:
                print(f"  Position: {h.position}")
                print(f"  Stem:     {h.stem_length} bp")
                print(f"  Energy:   {h.free_energy} kcal/mol")
                print(f"  Sequence: {h.sequence}")
        else:
            print("No hairpins detected.")
    
    return 0


def cmd_rare(args):
    """Identify rare codons."""
    rscu = ORGANISM_RSCU.get(args.organism, ORGANISM_RSCU['e_coli'])
    rare = identify_rare_codons(args.sequence, rscu, args.threshold)
    
    if args.json:
        output = [
            {
                'position': r.position, 'codon': r.codon,
                'amino_acid': r.amino_acid, 'rscu': r.rscu,
            }
            for r in rare
        ]
        print(json.dumps(output, indent=2))
    else:
        print(f"Rare Codons ({args.organism}, threshold={args.threshold})")
        print("=" * 50)
        if rare:
            print(f"{'Pos':>5} {'Codon':>6} {'AA':>3} {'RSCU':>8}")
            print("-" * 25)
            for r in rare:
                print(f"{r.position:>5} {r.codon:>6} {r.amino_acid:>3} {r.rscu:>8.4f}")
            print(f"\nTotal: {len(rare)} rare codons")
        else:
            print("No rare codons found.")
    
    return 0


def cmd_optimize(args):
    """Optimize codons for target organism."""
    result = full_optimization(args.protein, args.dna, args.organism)
    
    output = {
        'original_cai': result.original_cai,
        'optimized_cai': result.optimized_cai,
        'original_gc': round(result.original_gc, 4),
        'optimized_gc': round(result.optimized_gc, 4),
        'num_changes': result.num_changes,
        'optimized_sequence': result.optimized_sequence,
    }
    
    if args.json:
        output['changes'] = [
            {'position': p, 'original': o, 'optimized': n}
            for p, o, n in result.changes
        ]
        print(json.dumps(output, indent=2))
    else:
        print(f"Codon Optimization ({args.organism})")
        print("=" * 50)
        print(f"  Original CAI:  {result.original_cai:.6f}")
        print(f"  Optimized CAI: {result.optimized_cai:.6f}")
        print(f"  Original GC:   {result.original_gc:.4f}")
        print(f"  Optimized GC:  {result.optimized_gc:.4f}")
        print(f"  Changes:       {result.num_changes}")
        print(f"\n  Optimized DNA:")
        print(f"  {result.optimized_sequence}")
    
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='codon-optimization-species-agent',
        description='Codon Optimization: CAI, GC content, hairpin detection, rare codons',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # CAI
    p_cai = subparsers.add_parser('cai', help='Calculate Codon Adaptation Index')
    p_cai.add_argument('--sequence', required=True, help='DNA coding sequence')
    p_cai.add_argument('--organism', default='e_coli', choices=list(ORGANISM_RSCU.keys()))
    p_cai.add_argument('--verbose', action='store_true', help='Show per-codon scores')
    p_cai.add_argument('--json', action='store_true', help='JSON output')
    
    # GC
    p_gc = subparsers.add_parser('gc', help='GC content analysis')
    p_gc.add_argument('--sequence', required=True, help='DNA sequence')
    p_gc.add_argument('--json', action='store_true', help='JSON output')
    
    # Hairpin
    p_hp = subparsers.add_parser('hairpin', help='Detect mRNA hairpins')
    p_hp.add_argument('--sequence', required=True, help='DNA/RNA sequence')
    p_hp.add_argument('--min-stem', type=int, default=4, help='Minimum stem length')
    p_hp.add_argument('--max-loop', type=int, default=10, help='Maximum loop size')
    p_hp.add_argument('--json', action='store_true', help='JSON output')
    
    # Rare
    p_rare = subparsers.add_parser('rare', help='Identify rare codons')
    p_rare.add_argument('--sequence', required=True, help='DNA coding sequence')
    p_rare.add_argument('--organism', default='e_coli', choices=list(ORGANISM_RSCU.keys()))
    p_rare.add_argument('--threshold', type=float, default=0.3, help='RSCU threshold')
    p_rare.add_argument('--json', action='store_true', help='JSON output')
    
    # Optimize
    p_opt = subparsers.add_parser('optimize', help='Optimize codons')
    p_opt.add_argument('--protein', required=True, help='Protein sequence')
    p_opt.add_argument('--dna', required=True, help='Original DNA sequence')
    p_opt.add_argument('--organism', default='e_coli', choices=list(ORGANISM_RSCU.keys()))
    p_opt.add_argument('--json', action='store_true', help='JSON output')
    
    args = parser.parse_args(argv)
    
    commands = {
        'cai': cmd_cai,
        'gc': cmd_gc,
        'hairpin': cmd_hairpin,
        'rare': cmd_rare,
        'optimize': cmd_optimize,
    }
    
    return commands[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
