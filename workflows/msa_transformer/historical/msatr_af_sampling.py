#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run MSA Transformer on a single MSA file (Stockholm or FASTA format).
Generates contact predictions and saves them to the current directory.

Author: schaferjw (modified)
"""

import argparse
import os
from pathlib import Path
from Bio import AlignIO

from modules.MSA_TR_af_sampling import MSA_TR

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MSA Transformer on a single MSA file.")
    parser.add_argument("--msa", type=str, default=None, help="MSA file in Stockholm format")
    parser.add_argument("--fasta", type=str, default=None, help="Optional MSA file in FASTA format")
    parser.add_argument("--b_factor", type=str, default='n', help="(y/n) Color contact map based on pLDDT scores")
    parser.add_argument("--Beg", type=int, default=0, help="Trim residues from beginning of MSA")
    parser.add_argument("--End", type=int, default=0, help="Trim residues from end of MSA")
    parser.add_argument("--nseq", type=int, default=256, help="Number of sequences to subsample (diverse selection)")
    args = parser.parse_args()

    # Determine input format and read alignment
    msa_path = args.fasta if args.fasta else args.msa
    msa_format = "fasta" if args.fasta else "stockholm"
    align = AlignIO.read(msa_path, msa_format)

    # Extract query length
    seqs = [str(record.seq) for record in align]
    query_length = len(seqs[0].replace('-', ''))

    # Derive output name from input FASTA/MSA filename
    name = Path(msa_path).stem
    files_tr = [Path(msa_path).name]
    msa_dir = str(Path(msa_path).parent)
    output = [2.0]

    print(f"Running MSA Transformer on: {msa_path}")
    print(f"Query length: {query_length}")
    print(f"Subsampling {args.nseq} sequences (diverse selection)")
    print(f"Saving results as: full_{name}.csv")
    MSA_TR(msa_dir, files_tr, output, query_length, args.Beg, args.End, name, nseq=args.nseq)
