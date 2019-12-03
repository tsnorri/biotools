#
# Copyright (c) 2019 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
from Bio import SeqIO
import itertools
import sys
import vcf
import vcf.utils

# Input: two or more VCF files
# Output: single VCF file with the variant at each position taken from the first input that has a variant at that position.

def combine_vcf(input_files, output_file, chr_id = None, reference = None, reference_seq = None):
	first_reader = vcf.Reader(input_files[0])
	writer = vcf.Writer(output_file, first_reader)

	ref_rec = None
	ref_seq = reference_seq
	if ref_seq is None and reference is not None:
		ref_rec = next(SeqIO.parse(reference, format = "fasta"))
		print("Using record “%s” as the reference." % ref_rec.id, file = sys.stderr)
		ref_seq = ref_rec.seq

	for records in vcf.utils.walk_together(*itertools.chain((first_reader,), map(lambda x: vcf.Reader(x), input_files[1:])), vcf_record_sort_key = lambda x: (x.POS,)):
		for rec in records:
			if rec is not None:
				# Mangle REF and check.
				if ref_seq is not None:
					pos0 = rec.POS - 1
					current_ref_part = rec.REF
					ref_part = ref_seq[pos0:(pos0 + len(current_ref_part))]
					# FIXME compare ref_part to each ALT and fix the GT values.
					assert '-' not in ref_part
					rec.REF = ref_part

				# Mangle CHROM.
				if chr_id is not None:
					rec.CHROM = chr_id
				writer.write_record(rec)
				break

if __name__ == "__main__":
	parser = argparse.ArgumentParser("Combine multiple VCF files, use the first one as a template.")
	parser.add_argument('input', nargs = '+', type = argparse.FileType('rU'), help = "Input VCF files.")
	parser.add_argument('--chr', type = str, default = None, help = "Replace the chromosome identifier with the given value.")
	parser.add_argument('--reference', type = argparse.FileType('rU'), default = None, help = "Replace REF column values using the given reference.")
	args = parser.parse_args()

	combine_vcf(args.input, sys.stdout, chr_id = args.chr, reference = args.reference)
