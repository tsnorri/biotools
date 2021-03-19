# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
import os
import sys
from Bio import SeqIO

# Input: hla.dat from IMGT.
# Ouput: 
#	1.	Tab-separated file that contains the following:
#		Gene name (one of HLA-A, -B, -C, -DQA1, -DQB1, -DRB1)
#		Allele identifier
#		HLA type identifier
#		Exon co-ordinates from gene start as zero-based half-open intervals formatted as follows:
#			exon_number_1:start_pos-end_pos,exon_number_2:start_pos-end_pos
#		Exon sequence
#	2.	FASTA


parser = argparse.ArgumentParser("Output HLA-A, -B, -C, -DQA1, -DQB1, -DRB1 gene sequences from IMGT data.")
parser.add_argument('--input', required = True, type = argparse.FileType('rU'))
parser.add_argument('--genes', type = argparse.FileType('rU'))
parser.add_argument("--fasta", action = 'store_true', default = False)
parser.add_argument("--separate-files", action = 'store_true', default = False)	# Store each gene into a separate file.
parser.add_argument("--exons-only", action = 'store_true', default = False)
parser.add_argument("--separate-exons", action = 'store_true', default = False)	# Store each exon into a separate file (requires separate-files).
args = parser.parse_args()


# Default output.
gene_names = {
	"HLA-A":	(2, 3),
	"HLA-B":	(2, 3),
	"HLA-C":	(2, 3),
	"HLA-DQA1":	(2,),
	"HLA-DQB1":	(2,),
	"HLA-DRB1":	(2,)
}


def handle_record(rec, gene_name, exon_numbers, output_fasta, exons_only, separate_exons, dst_files):
	features = record.features
	seq = record.seq
	exon_descs = []

	# Extract exon co-ordinates from gene start.
	for feature in features:
		if feature.type == 'exon':
			qualifiers = feature.qualifiers
			numbers = qualifiers["number"]
			assert 1 == len(numbers)
			number = int(numbers[0])

			if number in exon_numbers:
				loc = feature.location
				desc = (number, loc.start, loc.end)
				exon_descs.append(desc)

	# If exons only were requested, check their co-ordinates and read substrings from seq.
	exon_seqs = {}
	if exons_only:
		for desc in exon_descs:
			(num, start, end) = desc
			end = 1 + end
			subseq = seq[start:end]
			exon_seqs[num] = str(subseq)

	exon_descs_formatted = ["%d:%d-%d" % (num, start, end) for (num, start, end) in exon_descs]

	# Attempt to parse the HLA type_name.
	desc_parts = rec.description.split(",")
	hla_type_parts = desc_parts[0].split("-")
	hla_type = hla_type_parts[1]
	
	if separate_exons:
		for (exon_number, exon_seq) in exon_seqs.items():
			if 0 < len(exon_seq):
				file_identifier = "%s-%s" % (gene_name, exon_number)
				dst_file = sys.stdout
				if file_identifier in dst_files:
					dst_file = dst_files[file_identifier]

				if output_fasta:
					print(">%s\t%s" % (rec.id, len(exon_seq)), file = dst_file)
					print(exon_seq, file = dst_file)
				else:
					print("%s\t%s\t%s\t%s" % (rec.id, hla_type, len(exon_seq), exon_seq), file = dst_file)
	else:
		if exons_only:
			seq = "".join(exon_seqs.values())
		if 0 < len(seq):
			dst_file = sys.stdout
			if gene_name in dst_files:
				dst_file = dst_files[gene_name]
			if output_fasta:
				print(">%s\t%s\t%s\t%s" % (gene_name, rec.id, ",".join(exon_descs_formatted), len(seq)), file = dst_file)
				print(seq, file = dst_file)
			else:
				print("%s\t%s\t%s\t%s\t%s\t%s" % (gene_name, rec.id, hla_type, ",".join(exon_descs_formatted), len(seq), seq), file = dst_file)


# Read the gene names if needed.
if args.genes:
	gene_names = {}
	for line in args.genes:
		if line.startswith("#"):
			continue

		fields = line.strip().split("\t")
		gene_name = fields[0]
		exon_numbers = fields[1:]
		exon_numbers = [int(x) for x in exon_numbers]
		gene_names[gene_name] = exon_numbers


# If separate files were requested, open files only if they do not exist.
dst_files = {}
if args.separate_files:
	flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
	file_handles = {}

	for gene_name in gene_names:
		fname = None

		identifiers = None
		if args.separate_exons:
			identifiers = ["%s-%s" % (gene_name, x) for x in gene_names[gene_name]]
		else:
			identifiers = [gene_name]

		for identifier in identifiers:
			if args.fasta:
				fname = "%s.fa" % identifier
			else:
				fname = "%s.tsv" % identifier
		
			fd = os.open(fname, flags)
			file_handles[identifier] = fd

	# Should use with statement instead.
	dst_files = {key: os.fdopen(fd, 'w') for (key, fd) in file_handles.items()}


# Parse the input.
for record in SeqIO.parse(args.input, "imgt"):
	# Check the gene name.
	for feature in record.features:
		if feature.type == 'CDS':
			genes = feature.qualifiers["gene"]
			assert 1 == len(genes)
			gene_name = genes[0]
			if gene_name in gene_names:
				#print("Record: %s gene: %s" % (record.id, gene_name), file = sys.stderr)
				handle_record(record, gene_name, gene_names[gene_name], args.fasta, args.exons_only, args.separate_exons, dst_files)
