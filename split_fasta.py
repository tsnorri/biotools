# 
# Copyright (c) 2019 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import gzip
import os
import sys

from Bio import SeqIO

# Input: FASTA.
# Output: FASTA, one sequence per file.


def open_file(fname, mode, use_gzip):
	if use_gzip:
		return gzip.open(fname, mode = mode)
	else:
		return open(fname, mode = mode)
		

def safe_filename(string):
	# Used the idea from https://stackoverflow.com/a/7406369/856976
	return "".join([c for c in string if c.isalpha() or c.isdigit()])

parser = argparse.ArgumentParser("Split the input FASTA into multiple files.")
parser.add_argument('--fasta', required = True, type = str, help = "Input file")
parser.add_argument('--gzip-input', action = 'store_true', help = "Read compressed input")
parser.add_argument('--gzip-output', action = 'store_true', help = "Compress the output")
args = parser.parse_args()

print("Reading FASTA…", file = sys.stderr)
with open_file(args.fasta, "rt", args.gzip_input) as src_file:
	for rec in SeqIO.parse(src_file, format = "fasta"):
		identifier = rec.id
		print("Handling record %s…" % identifier, file = sys.stderr)
		output_fname_fmt = "%s.fa.gz" if args.gzip_output else "%s.fa"
		output_fname = output_fname_fmt % safe_filename(identifier) 
		with open_file(output_fname, "xt", args.gzip_output) as dst_file:
			SeqIO.write([rec], dst_file, "fasta")
