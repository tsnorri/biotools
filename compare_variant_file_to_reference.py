#
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
from Bio import SeqIO
import sys
import vcf


if __name__ == "__main__":
	# Parse the command line arguments.
	parser = argparse.ArgumentParser(description = "Check the REF column of a variant file.")
	parser.add_argument('--vcf', type = argparse.FileType('r'), help = "Variant file")
	parser.add_argument('--reference', type = argparse.FileType('r'), help = "Reference FASTA")
	args = parser.parse_args()
	
	# Take the first sequence.
	fasta_record = next(SeqIO.parse(args.reference, format = "fasta"))
	reference_seq = fasta_record.seq
	
	vcf_reader = vcf.Reader(args.vcf)
	for row in vcf_reader:
		pos = row.POS
		
		actual = row.REF
		expected = reference_seq[pos - 1 : pos - 1 + len(actual)]
		if actual != expected:
			print("Mismatch at VCF position %d: expected '%s', got '%s'" % (pos, expected, actual), file = sys.stderr)
