#!/usr/bin/env python3

from Bio.AlignIO import read
import gzip as gz
from tqdm import tqdm
from argparse import ArgumentParser
import sys

def trim_gaps(alignment, max_frac=0.9):
    num_seqs = len(alignment)
    align_len = alignment.get_alignment_length()
    max_num = round(max_frac * num_seqs)
    sys.stderr.write(f"Trimming columns with {max_frac}*{num_seqs} = {max_num} or more gaps\n")
    cols = []
    # find columns to keep
    for col in tqdm(range(align_len-1), desc="analyzing columns", unit=" columns", leave=False):
        gaps = alignment[:, col].count('-')
        if gaps >= max_num:
            continue
        cols.append(col)
    sys.stderr.write(f"Found {len(cols)} columns to keep\n")
    records = {}
    for i, record in enumerate(tqdm(alignment, desc="trimming gaps", unit=" sequences", leave=False)):
        records[record.id] = "".join(alignment[i, col] for col in cols)
    return records

def main(input, output, frac):
    sys.stderr.write(f"Reading alignment from {input}\n")
    if input.endswith(".gz"):
        with gz.open('data/chesters_tree/chesters_expanded_aligned.fasta.gz', 'rt') as f:
            alignment = read(f, 'fasta')
    else:
        with open(input, 'r') as f:
            alignment = read(f, 'fasta')
    sys.stderr.write(f"Read alignment with {len(alignment)} sequences and {alignment.get_alignment_length()} positions\n")
    records = trim_gaps(alignment)
    sys.stderr.write(f"Writing trimmed alignment to {output}\n")
    with open(output, 'w') as f:
        for key, value in records.items():
            f.write(f">{key}\n{value}\n")


if __name__ == "__main__":
    parser = ArgumentParser(description="Trim gaps from an alignment")
    parser.add_argument("-i", "--input", help="Input alignment file")
    parser.add_argument("-o", "--output", help="Output alignment file")
    parser.add_argument("-f", "--frac", type=float, default=0.9, help="Maximum fraction of gaps to keep column")
    args = parser.parse_args()
    main(args.input, args.output, args.frac)
    