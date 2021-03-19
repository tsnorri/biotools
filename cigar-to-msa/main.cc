/*
 * Copyright (c) 2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <fstream>
#include <iostream>
#include <vector>
#include "cmdline.h"


namespace {
	
	void open_file(char const *path, std::ifstream &stream)
	{
		stream.exceptions(std::ifstream::badbit);
		stream.open(path);
		if (!stream.is_open())
		{
			std::cerr << "Unable to open file at path " << path << '\n';
			std::exit(EXIT_FAILURE);
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	// Open the inputs.
	std::ifstream query_stream;
	std::ifstream reference_stream;
	std::ifstream cigar_stream;
	
	open_file(args_info.query_arg, query_stream);
	open_file(args_info.reference_arg, reference_stream);
	open_file(args_info.cigar_arg, cigar_stream);
	
	// Prepare the output buffers.
	std::vector <char> query_output;
	std::vector <char> reference_output;
	
	
	// Read CIGAR until EOF.
	// Currently we do not check that the other inputs are no longer than the CIGAR string.
	std::size_t output_position{};
	std::size_t count{};
	char op{};
	while (cigar_stream >> count >> op)
	{
		switch (op)
		{
			case 'M': // Match or mismatch, consumes both
			case '=': // Match, consumes both
			case 'X': // Mismatch, consumes both
				query_output.resize(output_position + count, 0);
				reference_output.resize(output_position + count, 0);
				
				query_stream.read(query_output.data() + output_position, count);
				reference_stream.read(reference_output.data() + output_position, count);
				break;
			
			case 'I': // Insertion to the reference, consumes query
			case 'S': // Soft clipping, consumes query
				query_output.resize(output_position + count, 0);
				reference_output.resize(output_position + count, '-');
				
				query_stream.read(query_output.data() + output_position, count);
				break;
			
			case 'D': // Deletion from the reference, consumes reference
			case 'N': // Skipped region from the reference (intron for mRNA, otherwise undefined), consumes reference
				query_output.resize(output_position + count, '-');
				reference_output.resize(output_position + count, 0);
				
				reference_stream.read(reference_output.data() + output_position, count);
				break;
			
			case 'H': // Hard clipping, consumes nothing
			case 'P': // Padding, consumes nothing
				query_output.resize(output_position + count, '-');
				reference_output.resize(output_position + count, '-');
				break;
			
			default:
				throw std::runtime_error("Unexpected CIGAR operation");
				break;
		}
		
		output_position += count;
	}
	
	// Write the results to stdout as A2M.
	std::cout << ">reference\n";
	std::copy(reference_output.begin(), reference_output.end(), std::ostream_iterator<char>(std::cout));
	std::cout << '\n';
	
	std::cout << ">query\n";
	std::copy(query_output.begin(), query_output.end(), std::ostream_iterator<char>(std::cout));
	std::cout << '\n';
	
	std::cout << std::flush;
	
	return EXIT_SUCCESS;
}
