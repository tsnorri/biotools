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
	
	
	// Read from the stream, replace the contents of dst.
	void read_string(std::istream &stream, std::size_t const count, std::string &dst)
	{
		dst.resize(count);
		stream.read(dst.data(), count);
	}
	
	
	// Read from the stream, append to dst.
	void append_to_string(std::istream &stream, std::size_t const count, std::string &dst)
	{
		auto const current_pos(dst.size());
		dst.resize(current_pos + count);
		stream.read(dst.data() + current_pos, count);
	}
	
	
	enum class segment_record_type : std::uint8_t
	{
		UNKNOWN = 0,
		MATCH,
		MISMATCH,
		INSERTION,
		DELETION 
	};
	
	
	struct segment_record
	{
		std::string chromosome_identifier;
		std::string ref;
		std::string alt;
		std::size_t position{};
		segment_record_type record_type{segment_record_type::UNKNOWN};
		
		void output_record() const
		{
			// Do not output matching segments.
			if (segment_record_type::MATCH == record_type)
				return;
			
			// The reference part may be empty but that implies that alt is also empty.
			if (ref.empty())
			{
				assert(alt.empty());
				return;
			}
			
			// Output the VCF record.
			// CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample
			std::cout << chromosome_identifier << '\t' << (1 + position) << "\t.\t" << ref << '\t';
			
			if (segment_record_type::DELETION == record_type)
				std::cout << "<DEL>";
			else
				std::cout << alt;
				
			std::cout << "\t.\tPASS\t.\tGT\t1\n";
		}
	};
	
	
	class input_handler
	{
	protected:
		segment_record	m_prev_record;
			
		std::ifstream	m_query_stream;
		std::ifstream	m_reference_stream;
		std::ifstream	m_cigar_stream;
		
		std::size_t		m_reference_position{};
		
	public:
		void open_inputs(char const *query_path, char const *reference_path, char const *cigar_path)
		{
			open_file(query_path, m_query_stream);
			open_file(reference_path, m_reference_stream);
			open_file(cigar_path, m_cigar_stream);
		}
		
		
		void output_header(char const *sample_name) const
		{
			std::cout << "##fileformat=VCFv4.3\n";
			std::cout << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
			std::cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
			std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << '\n';
		}
		
		
		void process(char const *chr_id, char const *sample_name)
		{
			m_prev_record.chromosome_identifier = chr_id;
			
			output_header(sample_name);
			
			// Read CIGAR until EOF.
			// Currently we do not check that the other inputs are no longer than the CIGAR string.
			process_first_insertions();
			process_remaining();
		}
		
	protected:
		void read_query(std::size_t const count)
		{
			read_string(m_query_stream, count, m_prev_record.alt);
		}
		
		
		void read_reference(std::size_t const count)
		{
			read_string(m_reference_stream, count, m_prev_record.ref);
		}
		
		
		void append_to_query(std::size_t const count)
		{
			append_to_string(m_query_stream, count, m_prev_record.alt);
		}
		
		
		void append_to_reference(std::size_t const count)
		{
			append_to_string(m_reference_stream, count, m_prev_record.ref);
		}
		
		
		void handle_cigar_operation(std::size_t const count, char const op)
		{
			assert(count);
			switch (op)
			{
				case '=': // Match, consumes both
					m_prev_record.output_record();
					read_query(count);
					read_reference(count);
					m_prev_record.position = m_reference_position;
					m_prev_record.record_type = segment_record_type::MATCH;
					m_reference_position += count;
					break;
					
				case 'X': // Mismatch, consumes both
					m_prev_record.output_record();
					read_query(count);
					read_reference(count);
					m_prev_record.position = m_reference_position;
					m_prev_record.record_type = segment_record_type::MISMATCH;
					m_reference_position += count;
					break;
					
				case 'D': // Deletion from the reference, consumes reference
					m_prev_record.output_record();
					read_reference(count);
					m_prev_record.alt.clear();
					m_prev_record.position = m_reference_position;
					m_prev_record.record_type = segment_record_type::DELETION;
					m_reference_position += count;
					break;
					
				case 'I': // Insertion to the reference, consumes query
				{
					if (segment_record_type::INSERTION != m_prev_record.record_type)
					{
						// Remove the last character from the previous record.
						assert(!m_prev_record.ref.empty());
						std::string ref_last(1, m_prev_record.ref.back());
						m_prev_record.ref.pop_back();
						
						std::string alt_last;
						if (!m_prev_record.alt.empty())
						{
							alt_last = m_prev_record.alt.back();
							m_prev_record.alt.pop_back();
						}
						
						m_prev_record.output_record();
						
						m_prev_record.ref = ref_last;
						m_prev_record.alt = alt_last;
						m_prev_record.position = m_reference_position - 1;
						m_prev_record.record_type = segment_record_type::INSERTION;
					}
					
					// Append the insertion.
					append_to_query(count);
					break;
				}
					
				case 'M': // Match or mismatch, consumes both
					throw std::runtime_error("Only extended CIGAR is handled");
					break;
					
				default:
					throw std::runtime_error("Unexpected CIGAR operation");
					break;
			}
		}
		
		
		void process_first_insertions()
		{
			std::size_t count{};
			char op{};
			
			if (m_cigar_stream >> count >> op)
			{
				// Check if the first record is an insertion. If so,
				// continue appending to it until we get another type of record.
				assert(count);
				if ('I' == op)
				{
					append_to_query(count);
					while (m_cigar_stream >> count >> op)
					{
						assert(count);
						switch (op)
						{
							case 'I':
								// Continue appending while there are insertions.
								append_to_query(count);
								break;
							
							case 'X': // Mismatch, consumes both
							case '=': // Match, consumes both
								// Read one character from both and output the record.
								// Then continue with the remaining part of the operation.
								append_to_query(1);
								append_to_reference(1);
								m_prev_record.record_type = segment_record_type::MISMATCH;
								
								++m_reference_position;
								--count;
								if (count)
									handle_cigar_operation(count, op);
								goto end_initial_loop;
								
							case 'D': // Deletion from the reference, consumes reference
								// Read one character from the reference and output the record.
								// Then continue with the remaining part of the operation.
								append_to_reference(1);
								m_prev_record.record_type = segment_record_type::MISMATCH;
								
								++m_reference_position;
								--count;
								if (count)
									handle_cigar_operation(count, op);
								goto end_initial_loop;
								
							case 'M': // Match or mismatch, consumes both
								throw std::runtime_error("Only extended CIGAR is handled");
								break;
					
							default:
								throw std::runtime_error("Unexpected CIGAR operation");
								break;
						}
					}
				}
				else
				{
					handle_cigar_operation(count, op);
				}
			}
			
		end_initial_loop:
			;
		}
		
		
		void process_remaining()
		{
			std::size_t count{};
			char op{};
			
			while (m_cigar_stream >> count >> op)
				handle_cigar_operation(count, op);
			
			m_prev_record.output_record();
		}
	};
}


// FIXME: Handle M operations. This could be done by splitting the segments to matching and non-matching parts, and ensuring that the last reference part has at least one character.
int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	input_handler handler;
	handler.open_inputs(
		args_info.query_arg,
		args_info.reference_arg,
		args_info.cigar_arg
	);
	
	handler.process(args_info.chr_id_arg, args_info.sample_name_arg);
	
	return EXIT_SUCCESS;
}
