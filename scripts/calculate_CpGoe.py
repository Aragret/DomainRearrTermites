#!/usr/bin/python3

import argparse
from Bio import SeqIO

parser=argparse.ArgumentParser()

parser.add_argument("--input", help="Input fasta file")
parser.add_argument("--output", help="Output table with CpGoe values for each input entry")

args=parser.parse_args()

with open(args.output, 'w') as outfile:
    outfile.write('Name\tLength\tCG\tCpGoe\n')

    for rec in SeqIO.parse(args.input, "fasta"):
        CG = rec.seq.upper().count('CG')
        seq_len = len(rec.seq)
        frCG = CG/seq_len
        frC = rec.seq.upper().count('C')/seq_len
        frG = rec.seq.upper().count('G')/seq_len
        
        if frC == 0 or frG == 0:
            continue

        CpGoe = frCG/(frC*frG)
        
        outfile.write('{0}\t{1}\t{2}\t{3}\n'.format(rec.description, seq_len, frC+frG, CpGoe))
