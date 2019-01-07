//
// Copyright (c) 2019 Tuukka Norri
// This code is licensed under MIT license (see LICENSE for details).
//

#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include "cmdline.h"

typedef std::uint64_t line_number;


namespace {

	template <typename t_value>
	inline auto to_unsigned(t_value val) -> std::make_unsigned_t <t_value>
	{
		return val;
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);

	std::ios_base::sync_with_stdio(false);  // Don't use C style IO after calling cmdline_parser.

	// Store the lines to be skipped.
	std::set <line_number> skipped_lines;
	for (std::size_t i(0); i < args_info.lines_given; ++i)
	{
		auto const lineno(args_info.lines_arg[i]);
		if (lineno <= 0)
		{
			std::cerr << "Line numbers must be positive; got " << lineno << ".\n";
			continue;
		}

		auto const lineno_u(to_unsigned(lineno));
		if (std::numeric_limits <line_number>::max() < lineno_u)
		{
			std::cerr << "Unable to handle line number " << lineno << ".\n";
			continue;
		}

		skipped_lines.insert(lineno_u);
	}

	// Free the allocated memory.
	cmdline_parser_free(&args_info);

	// List the skipped line numbers as they were parsed.
	{
		std::cerr << "Skipping lines: ";
		bool first(true);
		for (auto const line : skipped_lines)
		{
			if (!first)
				std::cerr << ", ";
			std::cerr << line;
			first = false;
		}
		std::cerr << '\n';
	}

	// Handle the input.
	line_number lineno(0);
	std::string line;
	while (std::getline(std::cin, line))
	{
		++lineno;

		// Output status.
        if (0 == lineno % 10000000)
			std::cerr << "Line " << lineno << "…\n";

		if (0 == skipped_lines.count(lineno))
			std::cout << line << '\n';
	}

	return EXIT_SUCCESS;
}
