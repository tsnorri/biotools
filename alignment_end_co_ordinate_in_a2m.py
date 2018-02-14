# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
from Bio import SeqIO
import sys

# Parse the command line arguments.
parser = argparse.ArgumentParser(description = "Find the last co-ordinate of a sequence of given length in an A2M multiple alignment.")
parser.add_argument("--start", required = True, type = int)		# Start co-ordinate, zero-based.
parser.add_argument("--length", required = True, type = int)	# Sequence length.

args = parser.parse_args()
start = args.start
length = args.length

# Get the first record.
record = next(SeqIO.parse(sys.stdin, "fasta"))
seq = record.seq
valid_characters = ['A', 'C', 'G', 'T', 'N']

end = start
i = 0
while i < length:
	c = seq[end]
	end += 1
	
	if '-' == c:
		continue
	
	assert c in valid_characters
	i += 1
	
print("%d" % end)
