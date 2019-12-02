# 
# Copyright (c) 2018-2019 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

import argparse
from Bio import SeqIO
from collections import defaultdict
from functools import reduce
import itertools
import random
import sys


def get_sequences(f):
	seq_objs = [seq for seq in SeqIO.parse(f, "fasta")]
	seq_ids = [seq.id for seq in seq_objs]
	sequences = [seq.seq for seq in seq_objs]
	return seq_ids, sequences


def format_ref(seq):
	"""Remove gaps."""
	return str(seq).replace('-', '')


def format_nonempty_alt(seq, wildcard_handler):
	"""Remove gaps and replace N with REF values."""
	l = list(seq)
	wildcard_handler(l)
	retval = "".join(l)
	return retval


def format_alt(seq, wildcard_handler):
	"""Return <DEL> for empty ALTs, otherwise call format_nonempty_alt."""
	for c in seq:
		if '-' != c:
			return format_nonempty_alt(seq, wildcard_handler)
	return "<DEL>"


def format_sample_name(n):
	"""Remove unusual characters."""
	return n.translate(str.maketrans("*:/", "_--"))


def wildcard_handler_ref(seq, ref):
	assert len(ref) == len(seq)
	for i, c in enumerate(ref):
		if 'N' == seq[i]:
			seq[i] = c
		# Check for an added gap character.
		if '-' == seq[i]:
			seq[i] = ''


def wildcard_handler_random(seq):
	for i, c in enumerate(seq):
		if '-' == c:
			seq[i] = ''
		elif 'N' == c:
			seq[i] = random.choice(["A", "C", "G", "T"])


def handle_range(base_pos, gap_csum, chrom, sequences, seq_ids, wildcard_handling, rs, re, join_gt_by):
	ref, *rest = sequences
	ref_part = ref[rs:re]
	ref_part_formatted = format_ref(ref_part)
	assert not 'N' in ref_part_formatted
	
	# Format each ALT.
	alts = None
	if "reference" == wildcard_handling:
		wildcard_handler = lambda x: wildcard_handler_ref(x, ref_part)
		alts = [format_alt(seq[rs:re], wildcard_handler) for seq in rest]
	elif "random" == wildcard_handling:
		# Determine the unique ALTs for handling.
		alt_occs = reduce(lambda dst, kv: dst[str(kv[1])].append(kv[0]) or dst, enumerate(map(lambda x: x[rs:re], rest)), defaultdict(list))
		# Transform each key, combine unique values again.
		alt_occs = reduce(lambda dst, kv: dst[format_alt(kv[0], wildcard_handler_random)].extend(kv[1]) or dst, alt_occs.items(), defaultdict(list))

		# Calculate the resulting list size.
		res_count = reduce(lambda dst, x: dst + len(x), alt_occs.values(), 0)
		# Create the ALT list.
		alts = [None] * res_count
		for alt, idxs in alt_occs.items():
			for i in idxs:
				alts[i] = alt

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
	print("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s" % (chrom, base_pos + rs - gap_csum[rs], ref_part_formatted, ",".join(alts_in_order), join_gt_by.join([str(alt_idxs[alt]) for alt in alts])))


def find_seq_idxs(specific_sequences, seq_ids, include_first):
	if include_first:
		# Return REF
		yield 0

	# Handle the given ids.
	for seq_id in specific_sequences:
		try:
			yield 1 + seq_ids[1:].index(seq_id)
		except ValueError:
			print("Sequence with id %s not found." % seq_id, file = sys.stderr)
			sys.exit(1)


parser = argparse.ArgumentParser(description = "Transform a multiple alignment in A2M format into variants. The first sequence needs to be the reference.")
parser.add_argument('--input', type = argparse.FileType('r'), required = True, help = "Input MSA")
parser.add_argument('--chr', type = str, required = True, help = "Chromosome identifier")
parser.add_argument('--base-position', type = int, default = 1, help = "Base position to be added to the co-ordinates (default = 1)")
parser.add_argument('--mangle-sample-names', action = 'store_true', help = "Replace unusual characters in sample names")
parser.add_argument('--no-dels', action = 'store_true', help = "Do not use the <DEL> structural variant")
parser.add_argument('--specific-sequences', nargs = '*', type = str, default = [], action = "store", help = "Instead of writing one haploid sample for each input sequence, output one haploid or diploid donor using the given sequence identifier")
parser.add_argument('--omit-sequences', nargs = '*', type = str, default = [], action = "store", help = "Omit the given sequences from the output")
parser.add_argument('--wildcard-handling', choices = ['reference', 'random'], default = 'reference', help = "Specify how N values are to be handled")
parser.add_argument('--random-seed', type = int, default = 0, help = "Random seed (for use with --wildcard-handling=random)")
args = parser.parse_args()

if 2 < len(args.specific_sequences):
	print("At most two --specific-types arguments needed.", file = sys.stderr)
	sys.exit(1)

random.seed(args.random_seed)

seq_ids, sequences = get_sequences(args.input)

# Check for gap characters in the first column.
for s in sequences:
	if '-' == s[0]:
		print("Found a gap character in the first column.", file = sys.stderr)
		sys.exit(1)
	
# Make a boolean vector of positions that are followed by a gap.
gap_follows = len(sequences[0]) * [False]
for i, c in enumerate(sequences[0][1:]):
	gap_follows[i] = ('-' == c)

# Count the gaps up to each aligned position.
gap_csum = list(itertools.accumulate(itertools.chain([0], gap_follows[:-1])))

# Filter by --specific-types.
join_gt_by = "\t"
if 0 != len(args.specific_sequences):
	join_gt_by = "|"
	seq_idxs = list(find_seq_idxs(args.specific_sequences, seq_ids, True))
	seq_ids = [seq_ids[idx] for idx in seq_idxs]
	sequences = [sequences[idx] for idx in seq_idxs]
elif 0 != len(args.omit_sequences):
	seq_idxs = frozenset(find_seq_idxs(args.omit_sequences, seq_ids, False))
	seq_ids = [seq_ids[i] for i, _ in enumerate(seq_ids) if not(i in seq_idxs)]
	sequences = [sequences[i] for i, _ in enumerate(sequences) if not(i in seq_idxs)]

# Make a boolean vector of positions that are followed by a gap in any sequence.
gap_follows_in_any = None
if args.no_dels:
	gap_follows_in_any = len(sequences[0]) * [False]
	for i in range(len(sequences[0]) - 1):
		for s in sequences:
			if '-' == s[1 + i]:
				gap_follows_in_any[i] = True
				break

# Output the VCF header.
formatted_sample_names = None
if args.mangle_sample_names:
	formatted_sample_names = [format_sample_name(x) for x in seq_ids[1:]]
	assert len(formatted_sample_names) == len(frozenset(formatted_sample_names))
else:
	formatted_sample_names = seq_ids[1:]

print("##fileformat=VCFv4.2")
print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
if not args.no_dels:
	print('##ALT=<ID=DEL,Description="Deletion">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % ("\t".join(formatted_sample_names) if 0 == len(args.specific_sequences) else "SAMPLE1"))

# Range start, end.
rs = 0
re = 0
for i, next_is_gap in enumerate(gap_follows_in_any if args.no_dels else gap_follows):
	re = i + 1
	if next_is_gap:
		continue
	
	handle_range(args.base_position, gap_csum, args.chr, sequences, seq_ids, args.wildcard_handling, rs, re, join_gt_by)
	rs = re
