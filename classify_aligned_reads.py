#
# Copyright (c) 2020 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
import pysam
import sys


class Result(object):
	def __init__(self):
		self.unmapped_count = 0
		self.mapped_primary_count = 0


parser = argparse.ArgumentParser(description = 'Gather statistics from aligned reads')
parser.add_argument("--alignments", type = str, required = True, help = "Aligned reads")
args = parser.parse_args()


results_by_chr = {}


with pysam.AlignmentFile(args.alignments, 'r') as samfile:
	for i, seg in enumerate(samfile.fetch(until_eof = True), start = 1):
		if 0 == (i % 1000000):
			print("Segment %d…" % i, file = sys.stderr)

		rname = samfile.get_reference_name(seg.reference_id)
		res = results_by_chr.setdefault(rname, Result())

		if seg.is_unmapped:
			res.unmapped_count += 1
			continue

		# Is mapped.
		if not (seg.is_secondary or seg.is_supplementary):
			res.mapped_primary_count += 1


print("CHROM\tMAPPED_PRIMARY_COUNT\tUNMAPPED_COUNT", file = sys.stdout)
for chrom, res in results_by_chr.items():
	print("%s\t%d\t%d" % (chrom, res.mapped_primary_count, res.unmapped_count), file = sys.stdout)
