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
args = parser.parse_args()

print("Reading FASTA…", file = sys.stderr)
seq = None
for rec in SeqIO.parse(args.fasta, format = "fasta"):
	print("Found record: %s" % rec.id, file = sys.stderr)
	if rec.id == args.seq_id:
		print(">%s" % rec.id)
		# FIXME: this is really inefficient.
		for i in range(len(rec.seq)):
			sys.stdout.write(rec.seq[i])
		sys.stdout.write("\n")
		sys.exit(0)

print("Sequence identifier not found in FASTA.", file = sys.stderr)
sys.exit(os.EX_DATAERR)
