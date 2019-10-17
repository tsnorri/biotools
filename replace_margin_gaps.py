# 
# Copyright (c) 2018-2019 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys

from Bio import SeqIO

# Input: FASTA.
# Output: The sequences with gaps in the beginning and end replaced with N.

parser = argparse.ArgumentParser("Replace gap characters in the beginning and end of each sequence with N characters.")
parser.add_argument('--fasta', required = True, type = argparse.FileType('rU'))
args = parser.parse_args()

print("Reading FASTA…", file = sys.stderr)
seq = None
for rec in SeqIO.parse(args.fasta, format = "fasta"):
	print(">%s" % rec.id)
	
	seq_len = len(rec.seq)
	gap_part_length = seq_len
	for i in range(seq_len):
		if rec.seq[i] != '-':
			gap_part_length = i
			break
	
	# Check whether the end of the sequence was reached.
	if gap_part_length == seq_len:
		continue
	
	end_gap_part_length = seq_len
	for i in range(seq_len):
		if rec.seq[seq_len - 1 - i] != '-':
			end_gap_part_length = i
			break
	
	for i in range(gap_part_length):
		sys.stdout.write('N')
	for i in range(gap_part_length, seq_len - end_gap_part_length):
		sys.stdout.write(rec.seq[i])
	for i in range(end_gap_part_length):
		sys.stdout.write('N')
	
	sys.stdout.write("\n")
