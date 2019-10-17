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
	parser.add_argument('--chr', type = str, required = False, default = None, help = "Filter by chromosome")
	args = parser.parse_args()
	
	# Take the first sequence.
	fasta_record = next(SeqIO.parse(args.reference, format = "fasta"))
	print("Using sequence with ID '%s' as the reference." % fasta_record.id, file = sys.stderr)
	reference_seq = fasta_record.seq
	
	vcf_reader = vcf.Reader(args.vcf)
	for row in vcf_reader:
		if not (args.chr is None or row.CHROM == args.chr):
			continue

		pos = row.POS
		actual = str(row.REF).upper()
		expected = str(reference_seq[pos - 1 : pos - 1 + len(actual)]).upper()
		if actual != expected:
			print("Mismatch at VCF position %d: expected '%s', got '%s'" % (pos, expected, actual), file = sys.stderr)
