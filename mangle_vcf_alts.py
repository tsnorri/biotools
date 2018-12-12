#
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
import random
import sys
import vcf


def generate_sequence(length):
	return ''.join(random.choice("ACGT") for _ in range(length))


def handle_alts(rec):
	for allele in rec.alleles[1:]:
		if type(allele == vcf.model._Substitution):
			length = len(allele.sequence)
			allele.sequence = generate_sequence(length)
			# Check that we didn’t get the original sequence.
			while allele.sequence == rec.alleles[0]:
				allele.sequence = generate_sequence(length)
		elif type(allele == vcf.model._AltRecord and allele.type == 'DEL'):
			pass
		else:
			return False
	return True


if __name__ == "__main__":
	# Parse the command line arguments.
	parser = argparse.ArgumentParser(description = "Mangle insertions and MNVs in a variant file.")
	parser.add_argument('--input-vcf', type = argparse.FileType('r'), required = True, help = "Input VCF")
	parser.add_argument('--output-vcf', type = argparse.FileType('w'), required = True, help = "Output VCF")
	parser.add_argument('--chr', type = str, required = False, default = None, help = "Filter by chromosome")
	parser.add_argument('--random-seed', type = int, required = False, default = 0, help = "Random seed")
	args = parser.parse_args()
	
	# Set up the RNG.
	random.seed(args.random_seed)
	
	# Iterate the records. Output only snips and indels and generate random strings for those.
	vcf_reader = vcf.Reader(args.input_vcf)
	vcf_writer = vcf.Writer(args.output_vcf, template = vcf_reader)
	for rec in vcf_reader:
		if not (args.chr is None or rec.CHROM == args.chr):
			continue
		
		if handle_alts(rec):
			vcf_writer.write_record(rec)
		else:
			print("Skipping variant at position %d…" % rec.POS, file = sys.stderr)
