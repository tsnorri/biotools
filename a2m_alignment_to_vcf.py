# 
# Copyright (c) 2018-2019 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
from Bio import SeqIO
import sys


def get_sequences(f):
	seq_objs = [seq for seq in SeqIO.parse(f, "fasta")]
	seq_ids = [seq.id for seq in seq_objs]
	sequences = [seq.seq for seq in seq_objs]
	return seq_ids, sequences


def format_ref(seq):
	"""Remove gaps."""
	return str(seq).replace('-', '')


def format_nonempty_alt(ref, seq):
	"""Remove gaps and replace N with REF values."""
	assert len(ref) == len(seq)
	retval = list(seq)
	for i, c in enumerate(ref):
		if 'N' == retval[i]:
			retval[i] = c
		# Check for an added gap character.
		if '-' == retval[i]:
			retval[i] = ''
	retval = "".join(retval)
	return retval


def format_alt(ref, seq):
	"""Return <DEL> for empty ALTs, otherwise call format_nonempty_alt."""
	for c in seq:
		if '-' != c:
			return format_nonempty_alt(ref, seq)
	return "<DEL>"


def format_sample_name(n):
	"""Remove unusual characters."""
	return n.translate(str.maketrans("*:/", "_--"))


def handle_range(base_pos, chrom, sequences, seq_ids, rs, re, join_gt_by):
	ref, *rest = sequences
	ref_part = ref[rs:re]
	ref_part_formatted = format_ref(ref_part)
	assert not 'N' in ref_part_formatted
	
	# Format each ALT.
	alts = [format_alt(ref_part, seq[rs:re]) for seq in rest]

	# Get the unique values.
	unique_alts_without_ref = frozenset(filter(lambda x: x != ref_part_formatted, alts))
	if 0 == len(unique_alts_without_ref):
		return

	# Assign indices to them.
	alt_idxs = {alt: 1 + i for i, alt in enumerate(unique_alts_without_ref)}

	# Get the unique values in the order just determined.
	alts_in_order = [x[0] for x in sorted(alt_idxs.items(), key = lambda kv: kv[1])]

	# Add REF to alt_idxs.
	alt_idxs[ref_part_formatted] = 0

	# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
	print("%d\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s" % (chrom, base_pos + rs, ref_part_formatted, ",".join(alts_in_order), join_gt_by.join([str(alt_idxs[alt]) for alt in alts])))


parser = argparse.ArgumentParser(description = "Transform a multiple alignment in A2M format into variants. The first sequence needs to be the reference.")
parser.add_argument('--input', type = argparse.FileType('r'), required = True, help = "Input MSA")
parser.add_argument('--chr', type = int, required = True, help = "Chromosome number")
parser.add_argument('--base-position', type = int, default = 0, help = "Base position to be added to the co-ordinates")
parser.add_argument('--mangle-sample-names', action = 'store_true', help = "Replace unusual characters in sample names")
parser.add_argument('--specific-type', nargs = '*', type = str, default = [], action = "store", help = "Instead of writing one haploid sample for each HLA type, output one diploid donor with the given HLA types.")
args = parser.parse_args()

if len(args.specific_type) not in (0, 2):
	print(len(args.specific_type))
	print("Exactly zero or two --specific-type arguments needed.", file = sys.stderr)
	sys.exit(1)

seq_ids, sequences = get_sequences(args.input)

# Check for gap characters in the first column.
for s in sequences:
	if '-' == s[0]:
		print("Found a gap character in the first column.", file = sys.stderr)
		sys.exit(1)
	
# Make a bit vector of positions that are followed by a gap.
gap_follows = len(sequences[0]) * [False]
for i, c in enumerate(sequences[0][1:]):
	gap_follows[i] = ('-' == c)

# Filter by --specific-type.
join_gt_by = "\t"
if 0 != len(args.specific_type):
	join_gt_by = "|"
	lhs_seq_idx = None
	rhs_seq_idx = None
	try:
		lhs_seq_idx = 1 + seq_ids[1:].index(args.specific_type[0])
	except ValueError:
		print("Sequence with id %s not found." % args.specific_type[0], file = sys.stderr)
		sys.exit(1)
	try:
		rhs_seq_idx = 1 + seq_ids[1:].index(args.specific_type[1])
	except ValueError:
		print("Sequence with id %s not found." % args.specific_type[1], file = sys.stderr)
		sys.exit(1)
	seq_ids = [seq_ids[0], seq_ids[lhs_seq_idx], seq_ids[rhs_seq_idx]]
	sequences = [sequences[0], sequences[lhs_seq_idx], sequences[rhs_seq_idx]]

# Output the VCF header.
formatted_sample_names = None
if args.mangle_sample_names:
	formatted_sample_names = [format_sample_name(x) for x in seq_ids[1:]]
	assert len(formatted_sample_names) == len(frozenset(formatted_sample_names))
else:
	formatted_sample_names = seq_ids[1:]

print("##fileformat=VCFv4.2")
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
print("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % ("\t".join(formatted_sample_names) if 0 == len(args.specific_type) else "SAMPLE1"))

# Range start, end.
rs = 0
re = 0
for i, next_is_gap in enumerate(gap_follows):
	re = i + 1
	if next_is_gap:
		continue
	
	handle_range(args.base_position, args.chr, sequences, seq_ids, rs, re, join_gt_by)
	rs = re
