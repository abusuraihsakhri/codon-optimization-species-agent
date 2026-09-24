"""
CLI module for Codon Optimization.
"""
import argparse
import json
import sys

from .engine import (
    CODON_TABLE, ORGANISM_RSCU,
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


import csv
from pathlib import Path


def cmd_batch(args):
    """Run batch codon optimization or CAI evaluation on a CSV file."""
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file {args.input} not found.", file=sys.stderr)
        return 1

    rows_out = []
    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            seq_id = (row.get("sequence_id") or f"SEQ_{idx:03d}").strip()
            gene = (row.get("gene_name") or "unnamed").strip()
            org = (row.get("organism") or args.organism).strip()
            dna = (row.get("dna_sequence") or "").strip()
            prot = (row.get("protein_sequence") or "").strip()

            original_cai = 0.0
            original_gc = 0.0
            optimized_dna = ""
            optimized_cai = 0.0
            optimized_gc = 0.0
            num_changes = 0
            status = "ok"
            message = ""

            try:
                if org not in ORGANISM_RSCU:
                    valid = ", ".join(sorted(ORGANISM_RSCU))
                    raise ValueError(f"Unknown organism '{org}'. Choose one of: {valid}")
                if not dna and not prot:
                    raise ValueError("Provide dna_sequence, protein_sequence, or both")

                rscu = ORGANISM_RSCU[org]

                if dna:
                    cai_res = calculate_cai(dna, rscu)
                    original_cai = round(cai_res.cai, 4)
                    original_gc = round(calculate_gc_content(dna), 4)

                if prot:
                    optimized_dna = optimize_codons(prot, rscu)
                    optimized_cai = round(calculate_cai(optimized_dna, rscu).cai, 4)
                    optimized_gc = round(calculate_gc_content(optimized_dna), 4)

                    if dna:
                        full_res = full_optimization(prot, dna, target_organism=org)
                        num_changes = full_res.num_changes
                else:
                    normalized_dna = "".join(dna.split()).upper().replace("U", "T")
                    prot = "".join(
                        CODON_TABLE[normalized_dna[i:i + 3]]
                        for i in range(0, len(normalized_dna), 3)
                    )
                    optimized_dna = optimize_codons(prot, rscu)
                    optimized_cai = round(calculate_cai(optimized_dna, rscu).cai, 4)
                    optimized_gc = round(calculate_gc_content(optimized_dna), 4)
                    full_res = full_optimization(prot, normalized_dna, target_organism=org)
                    num_changes = full_res.num_changes

            except (ValueError, TypeError, NotImplementedError) as exc:
                status = "error"
                message = str(exc)

            rows_out.append({
                "sequence_id": seq_id,
                "gene_name": gene,
                "organism": org,
                "status": status,
                "message": message,
                "original_cai": original_cai,
                "optimized_cai": optimized_cai,
                "original_gc": original_gc,
                "optimized_gc": optimized_gc,
                "num_changes": num_changes,
                "optimized_dna": optimized_dna,
            })

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sequence_id", "gene_name", "organism", "status", "message",
        "original_cai", "optimized_cai",
        "original_gc", "optimized_gc",
        "num_changes", "optimized_dna",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    errors = sum(row["status"] == "error" for row in rows_out)
    print(
        f"Batch processing completed: {len(rows_out)} sequences written to {args.output}"
        + (f" ({errors} row error(s))" if errors else "")
    )
    return 0 if errors == 0 else 2

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

    # Batch
    p_batch = subparsers.add_parser('batch', help='Batch optimize sequences from CSV')
    p_batch.add_argument('-i', '--input', required=True, help='Input CSV file path')
    p_batch.add_argument('-o', '--output', required=True, help='Output CSV file path')
    p_batch.add_argument('--organism', default='e_coli', choices=list(ORGANISM_RSCU.keys()), help='Target host organism')
    
    args = parser.parse_args(argv)
    
    commands = {
        'cai': cmd_cai,
        'gc': cmd_gc,
        'hairpin': cmd_hairpin,
        'rare': cmd_rare,
        'optimize': cmd_optimize,
        'batch': cmd_batch,
    }
    
    try:
        return commands[args.command](args)
    except (ValueError, TypeError, NotImplementedError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())

