# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created to build a file processing utility for ARSC computations.
"""


import os
import re
import gzip
from ARSC.core import process_faa

# Remove extensions
def get_genome_name(path):
    base = os.path.basename(path)
    root, ext = os.path.splitext(base)  # ext = ".faa" or ".gz"
    if ext == ".gz":
        root, _ = os.path.splitext(root)  # remove ".faa"
    return root


# Return items: {"handle": <file_path>, "name": <genome_name>}
def collect_faa_files(input_path):
    # --- directory ---
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if f.endswith(".faa") or f.endswith(".faa.gz"):
                fpath = os.path.join(input_path, f)
                genome = get_genome_name(f)
                yield {"handle": fpath, "name": genome}
        return

    # --- single file ---
    if os.path.isfile(input_path):
        if input_path.endswith(".faa") or input_path.endswith(".faa.gz"):
            genome =  get_genome_name(input_path)
            yield {"handle": input_path, "name": genome}
            return
        else:
            raise ValueError(f"input must be .faa or .faa.gz: {input_path}")

    raise ValueError(f"input path does not exist: {input_path}")


def process_faa_auto(item, per_sequence=False):
    """
    item: {"handle": path_str, "name": genome_name}
    per_sequence: bool, whether to process sequences individually
    """
    handle = item["handle"]
    name = item["name"]

    if handle.endswith(".gz"):
        with gzip.open(handle, "rt") as f:
            return process_faa(f, name=name, per_sequence=per_sequence)
    else:
        return process_faa(handle, name=name, per_sequence=per_sequence)
