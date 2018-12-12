# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
from Bio import SeqIO
import os
import sys

# Input: FASTA.
# Output: The sequence lengths.

parser = argparse.ArgumentParser("Output the lengths of the sequences in a FASTA file.")
parser.add_argument('--fasta', required = True, type = argparse.FileType('rU'))
args = parser.parse_args()

for rec in SeqIO.parse(args.fasta, format = "fasta"):
	print("%s\t%d" % (rec.id, len(rec.seq)))
