# 
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
# 

# Add genotype information to a VCF file.

import sys

for line in sys.stdin:
	if line.startswith("##"):
		sys.stdout.write(line)
	elif line.startswith("#"):
		sys.stdout.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
		line = line.rstrip("\n")
		sys.stdout.write("%s\tFORMAT\tSAMPLE1\n" % line)
	else:
		line = line.rstrip("\n")
		sys.stdout.write("%s\tGT\t1\n" % line)
