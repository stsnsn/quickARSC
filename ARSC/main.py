# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created to compute the N/C/S-ARSC metrics as described in the following publication:
    Mende et al., Nature Microbiology, 2017 https://doi.org/10.1038/s41564-017-0008-3

Original citations for calculation metrics in Mende et al. 2017:
    Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic composition. Science 293 297–300.
    Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.
"""

import sys
import argparse
import os
from multiprocessing import Pool
from statistics import mean, stdev
from ARSC import __version__
from ARSC.utils import collect_faa_files, process_faa_auto, collect_fna_files, process_fna_pipeline, detect_nucleotide_file, dispatch_process
from ARSC.core import aa_dictionary

quickARSC_LOGO = """
             _      _              _____   _____  _____ 
            (_)    | |       /\   |  __ \ / ____|/ ____|
  __ _ _   _ _  ___| | __   /  \  | |__) | (___ | |     
 / _` | | | | |/ __| |/ /  / /\ \ |  _  / \___ \| |     
| (_| | |_| | | (__|   <  / ____ \| | \ \ ____) | |____ 
 \__, |\__,_|_|\___|_|\_\/_/    \_\_|  \_\_____/ \_____|
    | |                                                 
    |_|                                                 
          Compute ARSC (N/C/S) from FASTA files         
    Usage: arsc [options] <input.faa or input_directory>
"""

class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    parser = argparse.ArgumentParser(description=f"{quickARSC_LOGO}\n\n", formatter_class=CustomFormatter)
    parser.add_argument("input", help="Positional input: fasta file or directory")
    parser.add_argument("-p", "--per-sequence", action="store_true", help="Process each sequence individually")
    parser.add_argument("-a", "--aa-composition", action="store_true", help="Include amino acid composition ratios")
    parser.add_argument("-o", "--output", help="Output TSV file (default: stdout)")
    parser.add_argument("-s", "--stats", action="store_true", help="Output summary statistics to stderr")
    parser.add_argument("--no-auto-detection", action="store_true", help="Disable automatic nucleotide sequence detection (skip nucleotide-detected files)")
    parser.add_argument("--no-header", action="store_true", help="Suppress header line")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads")
    parser.add_argument("-d", "--decimal-places", default=6, type=int, help="Decimal places")
    parser.add_argument("--min-length", type=int, help="Minimum sequence length")
    parser.add_argument("--max-length", type=int, help="Maximum sequence length")
    parser.add_argument("-n", "--nucleotide", action="store_true", help="Nucleotide mode for calculate GC contents and ARSCs from fna/fna.gz (please install Prodigal)")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    # NOTE: -n / --nucleotide option was commented out
    if not hasattr(args, 'nucleotide'):
        args.nucleotide = False

    # 入力パス
    target_input = args.input
    if not target_input:
        parser.error("missing input: provide a .faa/.faa.gz file or directory")

    if args.stats and not os.path.isdir(target_input):
        parser.error("--stats can only be used with directory input")

    if args.nucleotide:
        items = list(collect_fna_files(args.input))
        initial_mode = 'fna'
    else:
        items = list(collect_faa_files(args.input))
        initial_mode = 'faa'

    print(f"quickARSC Version: {__version__}", file=sys.stderr)
    print(f"Found {len(items)} files to process.", file=sys.stderr)
    print(f"Using {args.threads} threads.", file=sys.stderr)

    # --- 各ファイルの処理モードを決め、並列で dispatch ---
    task_args = []
    if initial_mode == 'fna':
        for item in items:
            task_args.append((item, 'fna', args.per_sequence))
    else:
        for item in items:
            handle = item.get('handle')
            is_nuc = detect_nucleotide_file(handle)
            if is_nuc:
                if args.no_auto_detection:
                    # User requested to disable auto-detection: treat everything as amino-acid sequences
                    print(f"Note: {item.get('name')} looks like nucleotide sequences but --no-auto-detection set; treating as protein.", file=sys.stderr)
                    task_args.append((item, 'faa', args.per_sequence))
                else:
                    print(f"Warning: {item.get('name')} looks like nucleotide sequences — switching to nucleotide processing", file=sys.stderr)
                    task_args.append((item, 'fna', args.per_sequence))
            else:
                task_args.append((item, 'faa', args.per_sequence))

    with Pool(args.threads) as pool:
        results = pool.starmap(dispatch_process, task_args)

    # フィルタリング
    filtered_results = []
    for r in results:
        if 'error' in r:
            print(f"Skipping genome due to error: {r['error']}", file=sys.stderr)
            continue  # エラーがある場合は結果に含めない
        if args.per_sequence:
            r['sequences'] = [seq for seq in r['sequences'] if (args.min_length is None or seq['length'] >= args.min_length) and (args.max_length is None or seq['length'] <= args.max_length)]
            if r['sequences']:
                filtered_results.append(r)
        else:
            length = r.get('total_aa_length', 0)
            if (args.min_length is not None and length < args.min_length) or (args.max_length is not None and length > args.max_length):
                continue
            filtered_results.append(r)

    results = filtered_results
    print(f"After filtering: {len(results)} results.", file=sys.stderr)

    # 統計
    if args.stats:
        n_vals, c_vals, s_vals, mw_vals = [], [], [], []
        for r in results:
            if 'error' in r: continue
            target_list = r['sequences'] if args.per_sequence else [r]
            for data in target_list:
                n_vals.append(data.get('N_ARSC', 0))
                c_vals.append(data.get('C_ARSC', 0))
                s_vals.append(data.get('S_ARSC', 0))
                mw_vals.append(data.get('MW_ARSC', 0))
        
        if n_vals:
            print("\n" + "="*70, file=sys.stderr)
            label = "Per-Sequence" if args.per_sequence else "Per-File"
            print(f"SUMMARY STATISTICS ({label})".center(70), file=sys.stderr)
            print("="*70, file=sys.stderr)
            print(f"{'Metric':<12} {'Mean':<16} {'Stdev':<16} {'Min':<16} {'Max':<16}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            for name, vals in [('N_ARSC', n_vals), ('C_ARSC', c_vals), ('S_ARSC', s_vals), ('AvgResMW', mw_vals)]:
                sd = stdev(vals) if len(vals) > 1 else 0
                print(f"{name:<12} {mean(vals):<16.{args.decimal_places}f} {sd:<16.{args.decimal_places}f} {min(vals):<16.{args.decimal_places}f} {max(vals):<16.{args.decimal_places}f}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            print(f"{'Count':<12} {len(n_vals):<16}", file=sys.stderr)
            print("="*70 + "\n", file=sys.stderr)

    # 出力
    decimal_fmt = f"{{:.{args.decimal_places}f}}"
    out_handle = open(args.output, "w") if args.output else sys.stdout
    aa_keys = sorted(aa_dictionary.keys())

    try:
        if not args.no_header:
            h = ["query"]
            # -n が指定された場合に GC, base_ATGC を追加
            if args.nucleotide:
                h.extend(["genomic_GC", "base_A", "base_T", "base_G", "base_C"])
            # -p が指定された場合に sequence_id を追加
            if args.per_sequence:
                h.append("sequence_id")
            h.extend(["N_ARSC", "C_ARSC", "S_ARSC", "AvgResMW", "length" if args.per_sequence else "TotalLength"])
            if args.aa_composition:
                h.extend(aa_keys)
            out_handle.write("\t".join(h) + "\n")

        # データ行
        for r in results:
            if 'error' in r: continue

            # -p
            if args.per_sequence:
                for seq in r.get('sequences', []):
                    row = [r['genome']]

                    if args.nucleotide:
                        row.extend([
                            decimal_fmt.format(r.get('GC', 0)),
                            decimal_fmt.format(r.get('base_A', 0)),
                            decimal_fmt.format(r.get('base_T', 0)),
                            decimal_fmt.format(r.get('base_G', 0)),
                            decimal_fmt.format(r.get('base_C', 0))
                        ])

                    row.append(seq['sequence_id'])

                    row.extend([
                        decimal_fmt.format(seq.get('N_ARSC') or 0),
                        decimal_fmt.format(seq.get('C_ARSC') or 0),
                        decimal_fmt.format(seq.get('S_ARSC') or 0),
                        decimal_fmt.format(seq.get('MW_ARSC') or 0),
                        str(seq.get('length', 0))
                    ])

                    if args.aa_composition:
                        comp = seq.get('aa_composition', {})
                        row.extend([decimal_fmt.format(comp.get(aa, 0)) for aa in aa_keys])

                    out_handle.write("\t".join(row) + "\n")

            # 通常
            else:
                row = [r['genome']]
                # -n
                if args.nucleotide:
                    row.extend([
                        decimal_fmt.format(r.get('GC', 0)),
                        decimal_fmt.format(r.get('base_A', 0)),
                        decimal_fmt.format(r.get('base_T', 0)),
                        decimal_fmt.format(r.get('base_G', 0)),
                        decimal_fmt.format(r.get('base_C', 0))
                    ])
                row.extend([
                    decimal_fmt.format(r.get('N_ARSC') or 0),
                    decimal_fmt.format(r.get('C_ARSC') or 0),
                    decimal_fmt.format(r.get('S_ARSC') or 0),
                    decimal_fmt.format(r.get('MW_ARSC') or 0),
                    str(r.get('total_aa_length', 0))
                ])
                if args.aa_composition:
                    comp = r.get('aa_composition', {})
                    row.extend([decimal_fmt.format(comp.get(aa, 0)) for aa in aa_keys])
                out_handle.write("\t".join(row) + "\n")

    finally:
        if args.output:
            out_handle.close()
            print(f"Output written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()