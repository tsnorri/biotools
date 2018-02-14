/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iostream>
#include <map>
#include <regex>
#include <string>
#include "cmdline.h"


#define always_assert(X) do { if (!(X)) { std::cerr << "Assertion failed." << std::endl; abort(); }} while (false)


// Create an std::string_view from a pair of string iterators.
inline std::string_view sv_from_range(
	std::string const &str,
	std::string::const_iterator const &begin,
	std::string::const_iterator const &end
)
{
	auto const start_pos(std::distance(str.begin(), begin));
	auto const length(std::distance(begin, end));
	return std::string_view(str.data() + start_pos, length);
}


// Create an std::string_view from a regex match range.
inline std::string_view sv_from_range(
	std::string const &str,
	std::smatch::const_reference range
)
{
	return sv_from_range(str, range.first, range.second);
}


// Transparent comparison for std::string_views intended to be used with std::map,
// works by converting the std::string keys to string views.
struct less
{
	using is_transparent = std::true_type;

	bool operator()(std::string_view const &lhs, std::string_view const &rhs) const
	{
		return lhs < rhs;
	}
};


struct seq_data
{
	std::string seq;
	std::size_t length{};
};


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	// Don't use C-style IO.
	std::ios_base::sync_with_stdio(false);
	
	bool seen_first_gdna(false);
	std::string co_ordinate;
	std::string line;
	std::size_t lineno{};
	std::string const *first_seq_id{};
	std::string const *first_seq{};
	std::map <std::string, seq_data, less> sequences;
	std::regex start_re("^\\s([^ ]+)\\s+(.+)\\s*$", std::regex_constants::optimize);
	std::smatch match;

	while (std::getline(std::cin, line))
	{
		++lineno;

		// Check if the line may be parsed.
		if (std::regex_search(line, match, start_re))
		{
			// All eligible lines start with a space and the HLA type.
			// Remove the ones that have “gDNA” at the same place.
			auto const seq_id(sv_from_range(line, match[1]));
			if ("gDNA" == seq_id)
			{
				if (!seen_first_gdna)
				{
					seen_first_gdna = true;
					co_ordinate = match[2];
				}
			}
			else
			{
				// Find the current sequence or create a new one.
				auto it(sequences.find(seq_id));
				if (sequences.cend() == it)
					it = sequences.emplace(seq_id, seq_data{}).first;
				auto &sd(it->second);
				auto const listed_seq(sv_from_range(line, match[2]));
				
				// Check the current sequence is the first one.
				if (!first_seq_id)
				{
					first_seq_id = &it->first;
					first_seq = &it->second.seq;
				}
				
				bool const is_first(it->first == *first_seq_id);
				
				std::size_t i(sd.seq.size());
				for (auto const c : listed_seq)
				{
					switch (c)
					{
						case ' ':
						case '|':
							goto loop_end;
							
						case 'A':
						case 'C':
						case 'G':
						case 'T':
						case '*':
							sd.seq.append(1, c);
							++sd.length;
							break;
							
						case '-':
							always_assert(!is_first);
							sd.seq.append(1, (*first_seq)[i]);
							++sd.length;
							break;
							
						case '.':
							sd.seq.append(1, '-');
							break;
							
						default:
							always_assert(false);
					}
					++i;
					
				loop_end:
					;
				}
			}
		}
		
		if (0 == lineno % 1000)
			std::cerr << "Handled " << lineno << " lines…" << std::endl;
	}

	auto const alignment_len(first_seq->size());
	
	if (output_format_arg_a2m == args_info.output_format_arg)
	{
		// Convert starts to N.
		for (auto &kv : sequences)
		{
			for (auto &c : kv.second.seq)
			{
				if ('*' == c)
					c = 'N';
			}
		}
		
		// Output.
		std::cout << "; Co-ordinate: " << co_ordinate << '\n';
		std::cout << "; Alignment length: " << alignment_len << '\n';
		for (auto const &kv : sequences)
		{
			always_assert(alignment_len == kv.second.seq.size());
			std::cout << '>' << kv.first << ' ' << "len:" << kv.second.length << '\n';
			std::cout << kv.second.seq << '\n';
		}
	}
	else
	{
		std::cout << "# Co-ordinate: " << co_ordinate << '\n';
		std::cout << "# Alignment length: " << alignment_len << '\n';
		for (auto const &kv : sequences)
		{
			always_assert(alignment_len == kv.second.seq.size());
			std::cout << kv.first << '\t' << kv.second.length << '\t' << kv.second.seq << '\n';
		}
		std::cout << std::flush;
	}

	return 0;
}
