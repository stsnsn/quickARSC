# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created to compute the N/C-ARSC metrics as described in the following publication:
    Mende et al., Nature Microbiology, 2017 https://doi.org/10.1038/s41564-017-0008-3

Original citations for calculation metrics in Mende et al. 2017:
    Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic composition. Science 293 297–300.
    Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.
"""

import sys
import argparse
from multiprocessing import Pool
from ASTUR import __version__
from ASTUR.utils import collect_faa_files, process_faa_auto


ASTUR_LOGO = """
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;+;;xx+X;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xX+xxxxxX++;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xxxXXXXxXX+;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXxxX+;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;XXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXX+;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXx;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXX+;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXX+;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXX+;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXXx;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;+++;+XXXXXXXXXXXx;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXXXXXXx;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXX+;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;+xXXXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;+XXXXXXXXXXXXXx+XXXXXXXXx+;;;;;;;;
;;;;;;;;;;;;;;+XXXXXXXXXXXXXXx;;xXXXXXXx+;;;;;;;;;
;;;;;;;;;;;;+XXXXXXXXXXXXXXX+;;;;XXXxX;;;;;;;;;;;;
;;;;;;;;;+xXXXXXXXXXXXXXXXx;;;;;;;+;;;;;;;;;;;;;;;
;;;;;;;+XXXXXXXXXXXXXXXXx;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;++XXXXXXXXXXXXx+;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;xxxxXXXXXXXx+;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;xxXxXxX++;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;xx+X++;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 .d888888  .d88888b  d888888P dP     dP  888888ba 
d8'    88  88.    "'    88    88     88  88    `8b
88aaaaa88a `Y88888b.    88    88     88 a88aaaa8P'
88     88        `8b    88    88     88  88   `8b.
88     88  d8'   .8P    88    Y8.   .8P  88     88
88     88   Y88888P     dP    `Y88888P'  dP     dP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
"""

class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    parser = argparse.ArgumentParser(description=f"{ASTUR_LOGO}\n\nCompute ARSC from .faa/.faa.gz files", formatter_class=CustomFormatter)
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-i", "--input_dir", required=False, help="A faa or faa.gz file, or directory")
    parser.add_argument("input", nargs="?", help="Positional input: .faa/.faa.gz file or directory")
    parser.add_argument("-o", "--output", help="Output TSV file w/ header (optional). If omitted, print to stdout w/o header.")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads")


    args = parser.parse_args()

    # Allow either flag (-i/--input_dir) or positional input (not both)
    if args.input_dir is not None and args.input is not None:
        parser.error("cannot specify both positional input and -i/--input_dir; use one or the other")
    if args.input_dir is None and args.input is not None:
        args.input_dir = args.input
    if args.input_dir is None:
        parser.error("missing input: provide a .faa/.faa.gz file or directory")


    try:
        items = list(collect_faa_files(args.input_dir))
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"ASTUR Version: {__version__}", file=sys.stderr)
    print(f"Found {len(items)} files to process.", file=sys.stderr)
    print(f"Using {args.threads} threads.", file=sys.stderr)

    # Multiprocessing
    with Pool(args.threads) as pool:
        results = pool.map(process_faa_auto, items)


    # Output
    if args.output:
        with open(args.output, "w") as out:
            out.write("File\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW\n")
            for r in results:
                out.write(f"{r['genome']}\t{r['N_ARSC']}\t{r['C_ARSC']}\t{r['S_ARSC']}\t{r['MW_ARSC']}\n")
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        for r in results:
            print(f"{r['genome']}\t{r['N_ARSC']}\t{r['C_ARSC']}\t{r['S_ARSC']}\t{r['MW_ARSC']}")

if __name__ == "__main__":
    main()