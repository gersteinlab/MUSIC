#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ctype.h>
#include "ms_signal_track_tools.h"
#include "ms_annot_region_tools.h"
#include "ms_genomics_coords.h"
#include "ms_utils.h"
#include "ms_nomenclature.h"
#include <string.h>
#include <algorithm>
#include "ms_ansi_string.h"
#include "ms_annot_region_tools.h"
#include "ms_mapped_read_tools.h"

using namespace std; 

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

bool sort_read_lines(char* read1, char* read2)
{
	return(t_string::sort_strings(read1, read2));
}

#define MAX(x, y) ((x)>(y))?(x):(y)
#define MIN(x, y) ((x)<(y))?(x):(y)

#define L_CHROM (250*1000*1000)

vector<char*>* sort_bucket_read_lines(char* bucket_fp)
{
	// Load the reads.
	vector<char*>* bucket_read_lines = buffer_file(bucket_fp);	
	vector<int>* read_starts = new vector<int>();
	vector<t_read_line_sorting_info*>* sorting_info_list = new vector<t_read_line_sorting_info*>();
	for(int i_read = 0; i_read < (int)bucket_read_lines->size(); i_read++)
	{
		int cur_read_start = 0;
		sscanf(bucket_read_lines->at(i_read), "%*s %*s %d", &cur_read_start);

		t_read_line_sorting_info* cur_line_info = new t_read_line_sorting_info();
		cur_line_info->start = cur_read_start;
		cur_line_info->read_line = bucket_read_lines->at(i_read);

		sorting_info_list->push_back(cur_line_info);
	} // i_read loop.

	sort(sorting_info_list->begin(), sorting_info_list->end(), sort_read_line_info);
	vector<char*>* sorted_bucket_read_lines = new vector<char*>();

	for(int i_read = 0; i_read < (int)sorting_info_list->size(); i_read++)
	{
		sorted_bucket_read_lines->push_back(sorting_info_list->at(i_read)->read_line);

		delete sorting_info_list->at(i_read);
	} // i_read loop.

	delete(sorting_info_list);
	delete(read_starts);

	delete(bucket_read_lines);

	return(sorted_bucket_read_lines);
}

// Following is for sorting the mapped reads offline.
bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2)
{
	return(info1->start < info2->start);
}

void count_preprocessed_reads(char* mapped_reads_fp, int& n_F, int& n_R)
{
	n_F = 0;
	n_R = 0;

	FILE* f_mapped_reads = open_f(mapped_reads_fp, "r");
	char cur_line[1000];
	while(1)
	{
		if(fgets(cur_line, 1000, f_mapped_reads) == NULL)
		{
			break;
		}

		char cur_strand_char;
		if(sscanf(cur_line, "%*s %c", &cur_strand_char) != 1)
		{
			break;
		}

		if(cur_strand_char == 'F')
		{
			n_F++;
		}

		if(cur_strand_char == 'R')
		{
			n_R++;
		}
	} // file reading loop.
	fclose(f_mapped_reads);
}

// Generic preprocessing function for mapped read files.
void preprocess_mapped_reads_file(char* mrf_fp, char* parsed_reads_op_dir, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_id)
{
    // Divide SAM output with respect to chromosomes.
    FILE* f_mrf = NULL;
	t_file_buffer* mrf_file_buffer = NULL;
	if(strcmp(mrf_fp, "stdin") == 0)
	{
		f_mrf = stdin;
	}
	else
	{
		f_mrf = open_f(mrf_fp, "r");
	}

	if(f_mrf == NULL && mrf_file_buffer == NULL)
	{
		fprintf(stderr, "mapped read file pointer and file buffer are both NULL for %s\n", mrf_fp);
		return;
	}

    //char cur_line[100000];
    int n_frags = 0;
    //int n_total_frags = 0;

    vector<FILE*>* frag_f_ptrs = new vector<FILE*>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);

	vector<char*>* chr_ids = NULL;
	if(check_file(chr_ids_fp))
	{
		chr_ids = buffer_file(chr_ids_fp);

		fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

		// Open the files for appending.
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char new_fn[1000];
			sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
			fprintf(stderr, "Opening %s for pooling.\n", new_fn);
			if(!check_file(new_fn))
			{
				fprintf(stderr, "Could not open %s\n", new_fn);
				open_f(chr_ids_fp, "w");
				exit(0);
			}
			else
			{
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
			}
		} // i_chr loop.
	}
	else
	{
		// The chromosomes will be added now.
		chr_ids = new vector<char*>();
	}

	while(1)
	{
		//char* cur_line = getline(f_mrf);
		char* cur_line = NULL;
		if(f_mrf != NULL)
		{
			cur_line = getline(f_mrf);
		}

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char; 
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
									read_id,
									chrom, 
									chr_index, sequenced_length, 
									strand_char, 
									mapping_quality_str);

		// Make sure that the line is valid.
		if(chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the 
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if(i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));
				i_chr = t_string::get_i_str(chr_ids, chrom);

				char new_fn[10000];
				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chrom);

				// Does the file exist? If so, use the file, do not overwrite.
				frag_f_ptrs->push_back(open_f(new_fn, "w"));

				fprintf(stderr, "Added %s\n", chrom);
			}

			FILE* cur_frag_file = frag_f_ptrs->at(i_chr);

			if(cur_frag_file == NULL)
			{
					//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				if(dump_read_id)
				{
					fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
				}
				else
				{
					fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
				}
				n_frags++;
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for(int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_f_ptrs->size(); i_f++)
	{
		fclose(frag_f_ptrs->at(i_f));
	}

	// Unload/close the mapped read file.
	if(f_mrf != NULL)
	{
		fclose(f_mrf);
	}
	else if(mrf_file_buffer != NULL)
	{
		unload_file(mrf_file_buffer);
	}
}

void count_mapped_reads_per_file(char* mrf_fp, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	double n_total_reads)
{
    // Divide SAM output with respect to chromosomes.
    FILE* f_mrf = NULL;
	t_file_buffer* mrf_file_buffer = NULL;
	if(strcmp(mrf_fp, "stdin") == 0)
	{
		f_mrf = stdin;
	}
	else
	{
		f_mrf = open_f(mrf_fp, "r");
	}

	if(f_mrf == NULL && mrf_file_buffer == NULL)
	{
		fprintf(stderr, "mapped read file pointer and file buffer are both NULL for %s\n", mrf_fp);
		return;
	}

    //char cur_line[2000];
    //int n_frags = 0;
    //int n_total_frags = 0;
	n_total_reads = 0;

	fprintf(stderr, "Counting the mapped reads from %s\n", mrf_fp);
	while(1)
	{
		//char* cur_line = getline(f_mrf);
		char* cur_line = NULL;
		if(f_mrf != NULL)
		{
			cur_line = getline(f_mrf);
		}

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		//char chrom[1000];
		//int chr_index;
		//int sequenced_length;
		//char strand_char; 
		//char mapping_quality_str[1000];
		//preprocess_mapped_read_line(cur_line, 
		//							chrom, 
		//							chr_index, sequenced_length, 
		//							strand_char, 
		//							mapping_quality_str);
		char phred_quality_str[1000];
		int flag;
		if(cur_line[0] == '@')
		{
		}
		else
		{
			flag = 0;
			//if(sscanf(cur_line, "%*s %*d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 1)
			if(sscanf(cur_line, "%*s %d %*s %*d %*d %*s %*s %*d %*d %*s %s", &flag, phred_quality_str) == 2)
			{
				// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
				//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
				//chr_index = translate_coord(chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

				// Check the flag and determine the strand.
				/*strand_char = 'F';
				if(flag & 0x10)
				{
					strand_char = 'R';
				}*/

				// Sanity check. Is this fragment mapped?
				if(flag & 0x04)
				{
					// The read is not mapping.
					//chrom[0] = 0;
				}
				else
				{
					n_total_reads++;
				}
			}
			else
			{
				// Could not parse the line.
			}
		}

		if( (((int)n_total_reads) % 1000000) == 0)
		{
			fprintf(stderr, "At %lf. read.\n", n_total_reads);
		}

		delete [] cur_line;
	} // file reading loop.

	fprintf(stderr, "%lf reads.\n", n_total_reads);

	// Unload/close the mapped read file.
	if(f_mrf != NULL)
	{
		fclose(f_mrf);
	}
	else if(mrf_file_buffer != NULL)
	{
		unload_file(mrf_file_buffer);
	}
}

void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	int chr_start_index;
	int chr_end_index;
	char strand_sign;

	if(sscanf(cur_line, "%s %d %d %*s %*d %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	{
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		//chr_start_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		//chr_end_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		chr_start_index = translate_coord(chr_start_index, TAGALIGN_COORDS::start_base, CODEBASE_COORDS::start_base);
		chr_end_index = translate_coord(chr_end_index, TAGALIGN_COORDS::end_base, CODEBASE_COORDS::end_base);

		// Set quality to all matches.
		sprintf(cigar_str, "%dM", chr_end_index-chr_start_index+1);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = chr_end_index-chr_start_index+1;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int _chr_index;
	char strand;

	//if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 6)
	// X       14705460        35M     +
	if(sscanf(cur_line, "%s %d %s %c", chrom, &_chr_index, mapping_quality_str, &strand) == 4)
	{
		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		chr_index = _chr_index;

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand == '-')
		{
			strand_char = 'R';
		}

			chr_index = _chr_index;
			sequenced_length = 0;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char fragment[100000];
	char phred_quality_str[100000];

	//t_string_tokens* cur_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	//if(sscanf(cur_line, "%s %d %s %d %*s %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, cigar_str, fragment, phred_quality_str) == 7)
	//if(cur_tokens->size() >= 11)
	if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
	{
		//t_string::copy(read_id, cur_tokens->at(0)->str());
		//flag = atoi(cur_tokens->at(1)->str());
		//t_string::copy(chrom, cur_tokens->at(2)->str());
		//_chr_index = atoi(cur_tokens->at(3)->str());
		//t_string::copy(cigar_str, cur_tokens->at(5)->str());
		//t_string::copy(fragment, cur_tokens->at(9)->str());
		//t_string::copy(phred_quality_str, cur_tokens->at(10)->str());

		_chr_index = atoi(_chr_index_str);
		flag = atoi(flag_str);

		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if(flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			chr_index = _chr_index;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}

	//t_string::clean_tokens(cur_tokens);
}

void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char cur_fragment[100];
	char quality_str[100];
	int _chr_index;
	char _strand_char;

	if(sscanf(cur_line, "%s %s %s %*d %*d %*d %s %d %c", read_id, cur_fragment, quality_str, chrom, &_chr_index, &_strand_char) == 6)
	{
		chr_index = _chr_index;
		sequenced_length = strlen(cur_fragment);
		sprintf(mapping_quality_str, "%dM", sequenced_length);

		strand_char = 'F';
		if(_strand_char == '-')
		{
			strand_char = 'R';
		}
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	char nucs[1000];
	char strand_sign_str[100];
	char chr_start_index_str[100];
	//if(sscanf(cur_line, "%s %c %s %d %s", read_id, &strand_sign, chrom, &chr_start_index, nucs) == 4)
	if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t]", read_id, strand_sign_str, chrom, chr_start_index_str, nucs) == 5)
    {
		strand_sign = strand_sign_str[0];
		chr_start_index = atoi(chr_start_index_str);

		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		chr_start_index = translate_coord(chr_start_index, BOWTIE_COORDS::start_base, CODEBASE_COORDS::start_base);

		sprintf(mapping_quality_str, "%dM", (int)strlen(nucs));

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = strlen(nucs);
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;
	//if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	if(sscanf(cur_line, "%s %d %d %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		sprintf(mapping_quality_str, "%dM", chr_end_index-chr_start_index);
		sequenced_length = chr_end_index-chr_start_index;

		chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED5_read_line(char* cur_line, 
	char* read_id, 
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;
	if(sscanf(cur_line, "%s %d %d %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		sprintf(mapping_quality_str, "%dM", chr_end_index-chr_start_index);
		sequenced_length = chr_end_index-chr_start_index;

		chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED6_read_line(char* cur_line, 
	char* read_id, 
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;
	if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		sprintf(mapping_quality_str, "%dM", chr_end_index-chr_start_index);
		sequenced_length = chr_end_index-chr_start_index;

		chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
	}
	else
	{
		chrom[0] = 0;
	}
}

double get_n_mapped_nucs(vector<t_mapped_fragment*>* fragments)
{
	double n_mapped_nucs = 0.0;

	for(int i_frag = 0; i_frag < (int)fragments->size(); i_frag++)
	{
		n_mapped_nucs += fragments->at(i_frag)->sequenced_fragment_length;
	} // i_frag loop.

	return(n_mapped_nucs);
}

int get_l_signal_per_reads(char* reads_fp, int l_ext_tag)
{
	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	int l_profile = 0;
	FILE* f_sorted_reads = open_f(reads_fp, "r");
	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
		}	

		n_processed_reads++;

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int cur_l_profile = chr_index + 1000 + l_ext_tag;
		if(cur_l_profile > l_profile)
		{
			l_profile = cur_l_profile;
		}
	} // file loading loop.

	fclose(f_sorted_reads);

	return(l_profile);
} // get_l_signal_per_reads option.

// In order to ignore extended tag length, set it to non-positive value.
void buffer_per_nucleotide_profile_no_buffer(char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initilize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	// Non-positive tag extension lengths are ignored and falls back to using the length in the CIGAR string entry.
	if(l_extended_tag <= 0)
	{
		fprintf(stderr, "Ignoring tag extension.\n");
	}

	char strand_char = 0;
	//char cur_fragment[100000];
	char mapping_map_str[100000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = NULL;
	if(t_string::compare_strings(sorted_read_fp, "stdin"))
	{
		f_sorted_reads = stdin;
	}
	else
	{
		f_sorted_reads = open_f(sorted_read_fp, "r");
	}

	while(1)
	{
		char* cur_read = getline(f_sorted_reads);

		if(cur_read == NULL)
		{
			break;
		}
		else if((n_processed_reads % (1000*1000)) == 0)
		{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Processing %d. read.                        \r", n_processed_reads);
}
		}	

		n_processed_reads++;

		// We need a check on the current read line to make sure what we have is a valid: index must be strictly numbers; mapping map string is validated
		// below.
		char chr_index_str[1000];

		if(sscanf(cur_read, "%s %c %s", mapping_map_str, &strand_char, chr_index_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		// Validate the string read as chromosome index.
		int char_i = 0;
		while(chr_index_str[char_i] != 0)
		{
			bool char_is_a_num = (chr_index_str[char_i] >= '0' && chr_index_str[char_i] <= '9');
			if(!char_is_a_num)
			{
				fprintf(stderr, "Chromosome index must be a number: %s\n", cur_read);
				exit(0);
			}

			char_i++;
		} // char_i loop.

		// This may be pointing to an error in read preprocessing.
		chr_index = atoi(chr_index_str);
		if(chr_index > l_buffer)
		{
			fprintf(stderr, "%s: The read mapped out of buffer.\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int right_most_match_pos = 0;
		int left_most_match_pos = 1000*1000*1000;

		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			if(is_matching)
			{
				// Should we ignore the tag extension?
				if(l_extended_tag <= 0)
				{
					// Increase the height: The indexing is already 1 based, there is no conversion needed.
					if(signal_profile_buffer != NULL)
					{
						for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
						{
							signal_profile_buffer[i_cur]++;
						} // i_cur loop.
					}

					// Update the strand signals, if requested.
					if(forward_strand_signal != NULL)
					{
						if(strand_char == 'F')
						{
							// Update the forward strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								forward_strand_signal[i_cur]++;
							} // i_cur loop.
						}
						else
						{
							// Update the reverse strand signal.
							for(int i_cur = chr_index; i_cur <= chr_index+l_cur_entry-1; i_cur++)
							{
								reverse_strand_signal[i_cur]++;
							} // i_cur loop.
						}
					} // strand signal check.
				}
				else // tag extension validation check.
				{
					left_most_match_pos = MIN(left_most_match_pos, chr_index);
					right_most_match_pos = (chr_index + l_cur_entry - 1);
				} // extension length check.
			} // match block check.

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
			}
		} // mapping map string processing loop.

		// If tag extension is requested, utilize the leftmost and rightmost matching position for the mapped read.
		if(l_extended_tag > 0)
		{
			int ext_start = 0;
			int ext_end = 0;
			if(strand_char == 'F')
			{
				ext_start = left_most_match_pos;
			}
			else
			{
				//ext_start = (chr_index + l_cur_entry - 1) - (l_extended_tag - 1);
				ext_start = right_most_match_pos - (l_extended_tag - 1);
			}

			// Check for negative starts.
			if(ext_start < 0)
			{
				ext_start = 1;
			}

			ext_end = ext_start + l_extended_tag - 1;

			// Update profiles for the strands.
			if(forward_strand_signal != NULL)
			{
				if(strand_char == 'F')
				{
					// Update the forward strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						forward_strand_signal[i_cur]++;
					} // i_cur loop.
				}
				else
				{
					// Update the reverse strand signal.
					for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
					{
						reverse_strand_signal[i_cur]++;
					} // i_cur loop.
				}
			} // strand signal check.

			// Update the profile.
			if(signal_profile_buffer != NULL)
			{
				for(int i_cur = ext_start; i_cur <= ext_end; i_cur++)
				{
					signal_profile_buffer[i_cur]++;
				} // i_cur loop.
			} // profile existence check.
		} // signal update check.

		delete [] cur_read;
	} // file reading loop.

	if(t_string::compare_strings(sorted_read_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_sorted_reads);
	}

	// Get the length of the data.
	l_data = l_buffer;
	while(l_data > 0)
	{
		if(signal_profile_buffer != NULL && signal_profile_buffer[l_data] > 0.0)
		{
			break;
		}
		else if(forward_strand_signal != NULL && forward_strand_signal[l_data] > 0.0)
		{
			break;
		}
		else
		{
			l_data--;
		}
	} // l_data setting loop.
}

int fragment_5p_accessor(void* obj_ptr)
{
	t_mapped_fragment* frag_obj_ptr = (t_mapped_fragment*)obj_ptr;

	if(frag_obj_ptr->strand_char == 'F')
	{
		return(frag_obj_ptr->base_index);		
	}
	else if(frag_obj_ptr->strand_char == 'R')
	{
		return(frag_obj_ptr->base_index);
	}
	else
	{
		fprintf(stderr, "The strand char for fragment object is %c @ %s(%d)\n", frag_obj_ptr->strand_char, __FILE__, __LINE__);
		exit(0);
		return(-1);
	}
}

int fragment_3p_accessor(void* obj_ptr)
{
	t_mapped_fragment* frag_obj_ptr = (t_mapped_fragment*)obj_ptr;

	if(frag_obj_ptr->strand_char == 'F')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else if(frag_obj_ptr->strand_char == 'R')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else
	{
		fprintf(stderr, "The strand char for fragment object is %c @ %s(%d)\n", frag_obj_ptr->strand_char, __FILE__, __LINE__);
		exit(0);
		return(-1);
	}
}

int* get_n_reads_per_window(int n_wins, vector<t_mapped_fragment*>* frag_list)
{
	int* n_reads_per_window = new int[n_wins + 2];
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		n_reads_per_window[i_win] = 0;
	} // i_win loop

	for(int i_frag = 0; i_frag < (int)frag_list->size(); i_frag++)
	{
		int cur_i_win = frag_list->at(i_frag)->base_index / MEG_BASE;

		n_reads_per_window[cur_i_win]++;
	} // i_nuc loop.

	return(n_reads_per_window);
}

enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
	int i = 0;

	is_read_spliced = false;

	int state = VAL;
	while(mapping_map_str[i] != 0)
	{		
		if(state == VAL)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
			// MIDNSHPX=
			if(mapping_map_str[i] == 'M' ||
				mapping_map_str[i] == 'I' ||
				mapping_map_str[i] == 'D' ||
				mapping_map_str[i] == 'N' ||
				mapping_map_str[i] == 'S' ||		
				mapping_map_str[i] == 'H' ||
				mapping_map_str[i] == 'P' ||
				mapping_map_str[i] == 'X' ||
				mapping_map_str[i] == '=')
			{
				state = TYPE;

				if(mapping_map_str[i] != 'M')
				{
					is_read_spliced = true;
				}
			}
			else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				// State is still VAL.
			}
			else
			{
				return(false);
			}
		}
		else if(state == TYPE)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
			// A number is expected.
			if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				state = VAL;
			}
			else
			{
				return(false);
			}
		}

		// Move to next character.
		i++;
	}

	return(true);
}

void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										//t_string* cur_entry_length_str,
										int& l_cur_entry,
										char& entry_type_char)
{	
	// Clean the length string.
	//cur_entry_length_str->empty();
	l_cur_entry = 0;

	// Get the next entry in the cigar string.
	while(mapping_map_str[i_mapp_map] != 0)
	{
		if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
		{
			break;
		}
		//cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
		l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
		i_mapp_map++;
	}

	is_matching = false;
	if(mapping_map_str[i_mapp_map] == 'M')
	{
		//fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
		is_matching = true;
	}
	else
	{
		//fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
	}	

	entry_type_char = mapping_map_str[i_mapp_map];

	// Move over the current entry identifier.
	i_mapp_map++;
}

// Check if the 
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'H' ||
	//if(entry_char == 'P' ||
	//if(entry_char == 'X' ||
	if(entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}

bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'D' || 
	//	entry_char == 'M' ||
	//	entry_char == 'N' ||
	//	entry_char == 'H')
	if(entry_char == 'D' || 
		entry_char == 'M' ||
		entry_char == 'N')
	{
		return(true);
	}

	return(false);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}
