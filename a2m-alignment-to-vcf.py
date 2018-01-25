# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
from Bio import SeqIO
import sys


def output_vcf_record(vcf, chrom, pos, rec_idx, ref, alt_ids, alt_seqs_by_id_):
	assert 0 < len(ref)
	
	# Map ALT sequences to sequence identifiers that have an ALT.
	# Also process the sequences by removing '-'.
	alt_ids_by_seq = {}
	alt_seqs_by_id = {}
	for seq_id, seq in alt_seqs_by_id_.items():
		new_seq = "".join(filter(lambda c: '-' != c, seq))
		new_seq = new_seq if 0 < len(new_seq) else '<DEL>'
		
		alt_ids_by_seq.setdefault(new_seq, set()).add(seq_id)
		alt_seqs_by_id[seq_id] = new_seq
	
	## Debug
	#for k, v in alt_seqs_by_id_.items():
	#	sys.stdout.write("".join(v))
	#print("")
	#for k, v in alt_ids_by_seq.items():
	#	print("k: %s" % k)
	#	print("v: %d" % len(v))
	#	#for val in v:
	#	#	print("v: %s" % v)
	#sys.exit(0)
	
	# Map ALT sequences to indices.
	alt_seqs = [seq for seq in alt_ids_by_seq.keys()]
	alt_idxs_by_seq = { seq: i for i, seq in enumerate(alt_seqs) }

	# Output a record.
	vcf.write(str(chrom))				# CHROM
	vcf.write("\t")
	vcf.write(str(pos))					# POS
	vcf.write("\t")
	vcf.write("p%d" % rec_idx)			# ID
	vcf.write("\t")
	vcf.write("".join(ref))				# REF
	vcf.write("\t")
	vcf.write(",".join(alt_seqs))		# ALT
	vcf.write("\t")
	vcf.write(".")						# QUAL
	vcf.write("\t")
	vcf.write(".")						# FILTER
	vcf.write("\t")
	vcf.write(".")						# INFO
	vcf.write("\t")
	vcf.write("GT")						# FORMAT
	
	for alt_id in alt_ids:
		vcf.write("\t")
		alt_idx = 0
		if alt_id in alt_seqs_by_id:
			alt_seq = alt_seqs_by_id[alt_id]
			alt_idx = alt_idxs_by_seq[alt_seq]
			
		vcf.write(str(alt_idx))

	vcf.write("\n")


# Parse the command line arguments.
parser = argparse.ArgumentParser(description = "Convert A2M multialignment to VCF and FASTA.")
parser.add_argument("--chrom", required = True, type = int)
parser.add_argument('--fasta', required = True, type = argparse.FileType('w'))
parser.add_argument('--vcf', required = True, type = argparse.FileType('w'))

args = parser.parse_args()
chrom = args.chrom
fasta = args.fasta
vcf = args.vcf

# Read the sequences.
sequences = []
rec_ids = []
for record in SeqIO.parse(sys.stdin, "fasta"):
	sequences.append(record)
	rec_ids.append(record.id)

# Output VCF header.
vcf.write("##fileformat=VCFv4.2\n")
vcf.write("##source=a2m-alignment-to-vcf.py\n")

# Output the format header.
vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
for rec_id in rec_ids:
	vcf.write("\t")
	vcf.write(rec_id)
vcf.write("\n")

# Output the first sequence without the gaps.
gap_indices = []
first_record = sequences[0]
ref = first_record.seq
fasta.write(">")
fasta.write(first_record.id)
fasta.write("\n")
for i, c in enumerate(ref):
	if '-' == c:
		gap_indices.append(i)
	else:
		fasta.write(c)
fasta.write("\n")
gap_indices = frozenset(gap_indices)

# Output the VCF records.
insertions = 0
records = 0
current_rec_pos = 0
current_ref = []
current_alts = {}
for i, rc in enumerate(ref):
	
	# VCF recuires us to have a reference sequence for each insertion. Hence, do the following steps:
	# - Suppose that we're at index i. Maintain the reference sequence and ALTs at index i - 1.
	# - If ref[i] != '-', output the previous ALTs (in current_alts) and reset the variables.
	# - Otherwise, append the non-matching characters in to the ALT sequences.
	# - The case where position 1 contains a polymorphism is not handled.

	if '-' == rc:
		insertions += 1
	else:
		if 0 < len(current_alts):
			records += 1
			output_vcf_record(vcf, chrom, 1 + current_rec_pos, records, current_ref, rec_ids, current_alts)
			
		current_rec_pos = i - insertions
		current_ref = [rc]
		current_alts = {}
	
	# Check if the reference differs from any of the sequences
	# and update current_alts accordingly.
	have_diff = False
	it = iter(sequences)
	next(it) # Skip the reference.
	for record in it:
		rec_id = record.id
		ac = record.seq[i]
		seq = None
		if ac != rc:
			have_diff = True

			# Append ac to the ALTs.
			if rec_id in current_alts:
				seq = current_alts[rec_id]
			else:
				seq = current_ref[:-1]	# Slice up to but not including the last element.
				current_alts[rec_id] = seq
			seq.append(ac)

# Output the final record if needed.
if 0 < len(current_alts):
	records += 1
	output_vcf_record(vcf, chrom, 1 + current_rec_pos, records, current_ref, rec_ids, current_alts)
