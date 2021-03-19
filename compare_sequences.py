
# 
# Copyright (c) 2018-2020 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser("Compare the first sequences in two FASTA files.")
parser.add_argument('--lhs', required = True, type = str)
parser.add_argument('--rhs', required = True, type = str)
parser.add_argument('--lhs-is-txt', action = "store_true")
parser.add_argument('--continue-after-mismatch', action = "store_true")
parser.add_argument('--case-insensitive', action = "store_true")
args = parser.parse_args()


class Sequence(object):
	def __init__(self, seq):
		self.id = "<Unnamed>"
		self.seq = seq

def read_sequence(path, is_txt):
	if is_txt:
		with open(path, 'r') as f:
		    return Sequence(f.read())
	else:
		it = SeqIO.parse(path, format = "fasta")
		rec = next(it)
		return rec

rec1 = read_sequence(args.lhs, args.lhs_is_txt)
rec2 = read_sequence(args.rhs, False)

print(f"Comparing “{rec1.id}” to “{rec2.id}”…", file = sys.stderr)
if len(rec1.seq) != len(rec2.seq):
	print("Sequence lengths differ.", file = sys.stderr)
	sys.exit(1)

found_mismatch = False
for i, (c1, c2) in enumerate(zip(rec1.seq, rec2.seq)):
	cc1, cc2 = c1, c2
	if args.case_insensitive:
		cc1 = cc1.lower()
		cc2 = cc2.lower()
	if cc1 != cc2:
		print(f"Sequences differ at offset {i} ({c1} v. {c2})", file = sys.stderr)
		found_mismatch = True
		if not args.continue_after_mismatch:
			sys.exit(1)

if not found_mismatch:
	print("No differences found.", file = sys.stderr)
