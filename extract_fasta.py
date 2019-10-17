# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys

from Bio import SeqIO

# Input: FASTA.
# Output: The specified sequence.

parser = argparse.ArgumentParser("Output the specified sequence in a FASTA file.")
parser.add_argument('--fasta', required = True, type = argparse.FileType('rU'))
parser.add_argument('--seq-id', required = True)
parser.add_argument('--reverse-complement', action = "store_true")
parser.add_argument('--remove-gaps', action = "store_true")
args = parser.parse_args()

print("Reading FASTA…", file = sys.stderr)
seq = None
for rec in SeqIO.parse(args.fasta, format = "fasta"):
	print("Found record: %s" % rec.id, file = sys.stderr)
	if rec.id == args.seq_id:
		print(">%s" % rec.id)
		if args.reverse_complement:
			table = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-', 'N': 'N'}
			for i in reversed(range(len(rec.seq))):
				if args.remove_gaps and '-' == rec.seq[i]:
					continue
				sys.stdout.write(table[rec.seq[i]])
		else:
			# FIXME: this is really inefficient.
			for i in range(len(rec.seq)):
				if args.remove_gaps and '-' == rec.seq[i]:
					continue
				sys.stdout.write(rec.seq[i])
		sys.stdout.write("\n")
		sys.exit(0)

print("Sequence identifier not found in FASTA.", file = sys.stderr)
sys.exit(os.EX_DATAERR)
