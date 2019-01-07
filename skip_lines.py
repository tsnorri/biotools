#
# Copyright (c) 2018 Tuukka Norri
# This code is licensed under MIT license (see LICENSE for details).
#

import argparse
import sys


if __name__ == "__main__":
	# Parse the command line arguments.
	parser = argparse.ArgumentParser(description = "Output all the input lines except the ones listed.")
	parser.add_argument('--lines', metavar = 'line', type = int, nargs='+', help = 'skipped line numbers')
	args = parser.parse_args()

	skipped_lines = args.lines
	
	# Read from stdin line by line. If the line number matches one in skipped_lines, skip.
	for lineno, line in enumerate(sys.stdin, 1):
		if 0 == lineno % 1000000:
			print("Line %d…" % lineno, file = sys.stderr)

		if not(lineno in skipped_lines):
			sys.stdout.write(line)
