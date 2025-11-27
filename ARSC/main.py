# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created by removing unnecessary parts from the script provided in the following paper
    and by adding components such as S-ARSC or AvgResMW, as well as incorporating the multiprocessing module.
    Mende et al., Nature Microbiology, 2017 https://doi.org/10.1038/s41564-017-0008-3

Original citations for calculation metrics:
    Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic composition. Science 293 297–300.
    Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.
"""

# USAGE: python ARSC.py -i <input directory of protein files> -o <outputfile> -t <num_threads>
# e.g., python ARSC.py -i protein_faa/ -o ARSC.tsv -t 4

import os
import argparse
from ARSC.core import process_faa
from ARSC import __version__
from multiprocessing import Pool

# Main
# -----------------------
def main():
    parser = argparse.ArgumentParser(description="Compute ARSC from .faa files")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory of .faa files")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    faa_files = [os.path.join(args.input_dir, f)
                 for f in os.listdir(args.input_dir)
                 if f.endswith(".faa")]

    print(f"Found {len(faa_files)} faa files.")
    print(f"Using {args.threads} threads.")

    with Pool(args.threads) as pool:
        results = pool.map(process_faa, faa_files)

    # write output
    with open(args.output, "w") as out:
        out.write("Genome\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW\n")
        for r in results:
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n")

    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()