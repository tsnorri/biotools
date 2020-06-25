#
# Copyright (c) 2020 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
import pysam
import sys


parser = argparse.ArgumentParser(description = 'List the edit distances of the segments in a BAM file')
parser.add_argument("--alignments", type = str, required = True, help = "Aligned reads")
args = parser.parse_args()

with pysam.AlignmentFile(args.alignments, 'r') as samfile:
	for i, seg in enumerate(samfile.fetch(until_eof = True), start = 1):
		if 0 == (i % 1000000):
			print("Segment %d…" % i, file = sys.stderr)

		if not (seg.is_unmapped):
			rname = samfile.get_reference_name(seg.reference_id)
			print("%s\t%s\t%d" % (rname, seg.query_name, seg.get_tag("NM")), file = sys.stdout)
