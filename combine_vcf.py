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

parser = argparse.ArgumentParser("Combine multiple VCF files, use the first one as a template.")
parser.add_argument('input', nargs = '+', type = argparse.FileType('rU'), help = "Input VCF files.")
parser.add_argument('--chr', type = str, default = None, help = "Replace the chromosome identifier with the given value.")
parser.add_argument('--reference', type = argparse.FileType('rU'), default = None, help = "Replace REF column values using the given reference.")
args = parser.parse_args()

first_reader = vcf.Reader(args.input[0])
writer = vcf.Writer(sys.stdout, first_reader)

ref_rec = None
if args.reference is not None:
	ref_rec = next(SeqIO.parse(args.reference, format = "fasta"))
	print("Using record “%s” as the reference." % ref_rec.id, file = sys.stderr)

for records in vcf.utils.walk_together(*itertools.chain((first_reader,), map(lambda x: vcf.Reader(x), args.input[1:])), vcf_record_sort_key = lambda x: (x.POS,)):
	for rec in records:
		if rec is not None:
			# Mangle REF and check.
			if ref_rec is not None:
				pos0 = rec.POS - 1
				current_ref_part = rec.REF
				ref_part = ref_rec.seq[pos0:(pos0 + len(current_ref_part))]
				# FIXME compare ref_part to each ALT and fix the GT values.
				rec.REF = ref_part

			# Mangle CHROM.
			if args.chr is not None:
				rec.CHROM = args.chr
			writer.write_record(rec)
			break
