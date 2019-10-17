# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys

from Bio import SeqIO

# Input: FASTA.
# Output: The specified subsequence in the given sequence.

parser = argparse.ArgumentParser("Output the subsequence at the specified position in a sequence.")
parser.add_argument('--fasta', required = True, type = argparse.FileType('rU'))
parser.add_argument('--seq-id', required = True)
parser.add_argument('--pos', type = int)
parser.add_argument('--pos1', type = int)
parser.add_argument('--length', type = int)
args = parser.parse_args()

pos = 0
if args.pos1 is not None:
	pos = args.pos1 - 1
else:
	pos = args.pos

print("Reading FASTA…", file = sys.stderr)
seq = None
for rec in SeqIO.parse(args.fasta, format = "fasta"):
	print("Found record: %s" % rec.id, file = sys.stderr)
	if rec.id == args.seq_id:
		print("%d: %s" % (pos, rec.seq[pos : pos + args.length]))
		sys.exit(0)

print("Sequence identifier not found in FASTA.", file = sys.stderr)
sys.exit(os.EX_DATAERR)
