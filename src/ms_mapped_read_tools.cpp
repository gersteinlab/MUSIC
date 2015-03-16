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

void load_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = open_f(sam_fp, "r");
	}

	char cur_seq_id[1000];
	char cur_seq_nucs[1000];
	char phred_quality_str[1000];
	while(1)
	{
		char* cur_sam_line = getline(f_sam);
		if(cur_sam_line == NULL)
		{
			break;
		}

		// Parse the sam read line.
		if(sscanf(cur_sam_line, "%s %*d %*s %*d %*d %*s %*s %*d %*d %s %s", cur_seq_id, cur_seq_nucs, phred_quality_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_sam_line);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = t_string::copy_me_str(cur_seq_id);
		new_sequenced_read->nucs = t_string::copy_me_str(cur_seq_nucs);
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = t_string::copy_me_str(phred_quality_str);
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	// If the sam file is not stdin, close it.
	if(strcmp(sam_fp, "stdin") != 0)
	{
		fclose(f_sam);
	}
}

// FASTQ LOADING INTERFACE:
// Following are for fastq loading from mapped read files.
void load_sequenced_reads_per_fastq(char* fastq_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_fastq = NULL;
	if(strcmp(fastq_fp, "stdin") == 0)
	{
		f_fastq = stdin;
	}
	else
	{
		f_fastq = open_f(fastq_fp, "r");
	}

	while(1)
	{
		char* cur_seq_id = getline(f_fastq);
		if(cur_seq_id == NULL)
		{
			break;
		}

		char* cur_seq_nucs = getline(f_fastq);
		if(cur_seq_nucs == NULL)
		{
			fprintf(stderr, "Could not read the sequence for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_id = getline(f_fastq);
		if(cur_seq_qual_id == NULL)
		{
			fprintf(stderr, "Could not read the quality id for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_str = getline(f_fastq);
		if(cur_seq_qual_str == NULL)
		{
			fprintf(stderr, "Could not read the quality str for %s.\n", cur_seq_id);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = cur_seq_id;
		new_sequenced_read->nucs = cur_seq_nucs;
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = cur_seq_qual_str;
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	if(strcmp(fastq_fp, "stdin") != 0)
	{
		fclose(f_fastq);
	}
}

//void get_per_strand_read_stats(vector<t_annot_region*>* regions,
//	char* preprocessed_reads_dir)
//{
//	t_sorted_annot_region_lists* restructured_regions = restructure_annot_regions(regions);
//
//	char chr_ids_fp[1000];
//	sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
//	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
//
//	for(int i_reg_chr = 0; i_reg_chr < (int)restructured_regions->chr_ids->size(); i_reg_chr++)
//	{
//		// Load the fragments for the current chromosome.
//		int i_read_chr = t_string::get_i_str(chr_ids, restructured_regions->chr_ids->at(i_reg_chr));
//
//		if(i_read_chr < (int)chr_ids->size())
//		{
//			char mapped_reads_fp[1000];
//			sprintf(mapped_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_read_chr));
//
//			// Load the fragments for the current chromosome.
//			vector<t_mapped_fragment*>* fore_strand_fragments = new vector<t_mapped_fragment*>();
//			vector<t_mapped_fragment*>* rev_strand_fragments = new vector<t_mapped_fragment*>();
//			load_fragments(mapped_reads_fp, 
//				fore_strand_fragments, rev_strand_fragments, 
//				0);
//
//			// Sort the fragments before setting the 3p ends.
//			sort(fore_strand_fragments->begin(), fore_strand_fragments->end(), sort_mapped_fragments);
//			sort(rev_strand_fragments->begin(), rev_strand_fragments->end(), sort_mapped_fragments);
//
//            vector<int>* cur_chr_fore_fragment_3p_posns = new vector<int>();
//            for(int i = 0; i < (int)fore_strand_fragments->size(); i++)
//            {
//				cur_chr_fore_fragment_3p_posns->push_back(fore_strand_fragments->at(i)->base_index + fore_strand_fragments->at(i)->sequenced_fragment_length - 1);
//            } // i loop.
//
//            vector<int>* cur_chr_rev_fragment_3p_posns = new vector<int>();
//            for(int i = 0; i < (int)rev_strand_fragments->size(); i++)
//            {
//				//cur_chr_rev_fragment_3p_posns->push_back(rev_strand_fragments->at(i)->base_index);
//				cur_chr_rev_fragment_3p_posns->push_back(rev_strand_fragments->at(i)->base_index + rev_strand_fragments->at(i)->sequenced_fragment_length - 1);
//				//rev_strand_fragments->at(i)->base_index -= (rev_strand_fragments->at(i)->sequenced_fragment_length - 1);
//            } // i loop.
//
//			// Go over all the regions in this chromosome.
//			vector<t_annot_region*>* cur_chr_regions = restructured_regions->regions_per_chrom[i_reg_chr];
//			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
//			{
//				double cur_reg_n_mapped_fore_fragments = 0.0;
//				double cur_reg_n_mapped_fore_nucs = 0.0;
//				get_read_statistics_per_region(fore_strand_fragments, 
//													cur_chr_fore_fragment_3p_posns,
//													cur_chr_regions->at(i_reg), 
//													cur_reg_n_mapped_fore_fragments, 
//													cur_reg_n_mapped_fore_nucs);
//
//				double cur_reg_n_mapped_rev_fragments = 0.0;
//				double cur_reg_n_mapped_rev_nucs = 0.0;
//				get_read_statistics_per_region(rev_strand_fragments, 
//													cur_chr_rev_fragment_3p_posns, 
//													cur_chr_regions->at(i_reg), 
//													cur_reg_n_mapped_rev_fragments, 
//													cur_reg_n_mapped_rev_nucs);
//
//				// Set the read stats to the region.
//				double* cur_cnts = new double[2];
//				cur_cnts[0] = cur_reg_n_mapped_fore_nucs;
//				cur_cnts[1] = cur_reg_n_mapped_rev_nucs;
//				cur_chr_regions->at(i_reg)->data = cur_cnts;
//			} // i_reg loop.
//
//			// Free memory for the current chromosome.
//			delete_fragments(fore_strand_fragments);
//			delete_fragments(rev_strand_fragments);
//			delete(cur_chr_fore_fragment_3p_posns);
//			delete(cur_chr_rev_fragment_3p_posns);
//		} // read chr check.
//	} // i_reg_chr loop.
//}

void delete_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads)
{
	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		delete [] sequenced_reads->at(i_r)->id;
		delete [] sequenced_reads->at(i_r)->nucs;
		delete [] sequenced_reads->at(i_r)->quality_str;

		delete sequenced_reads->at(i_r);
	} // i_r loop.

	delete sequenced_reads;
}

void dump_fastq(vector<t_sequenced_read*>* sequenced_reads, char* op_fastq_fp)
{
	FILE* f_op = open_f(op_fastq_fp, "w");

	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		fprintf(f_op, "@%s\n%s\n+%s\n%s\n", 
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->nucs,
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->quality_str);
	} // i_r loop.

	fclose(f_op);
}

void dump_phred_quality_distribution(vector<t_sequenced_read*>* sequenced_reads, char* op_fp)
{
	vector<vector<char>*>* distributions_per_read_posn = new vector<vector<char>*>();
	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		if(i_r % 100000 == 0)
		{
			fprintf(stderr, "Processing %d. read.           \r", i_r);
		}

		// Process the reads.
		int l_cur_read = strlen(sequenced_reads->at(i_r)->nucs);

		// Add the entries if it is necessary to extend the posns.
		if(l_cur_read >= (int)distributions_per_read_posn->size())
		{
			for(int i = (int)distributions_per_read_posn->size(); i < l_cur_read; i++)
			{
				fprintf(stderr, "Adding %d. position to the distribution.\n", i);
				distributions_per_read_posn->push_back(new vector<char>());
			} // i loop.
		}

		// Process all the entries.
		for(int i = 0; i < l_cur_read; i++)
		{
			distributions_per_read_posn->at(i)->push_back(sequenced_reads->at(i_r)->quality_str[i]);
		} // i loop.
	} // i_r loop.

	// Dump the mean and variance for the qualities.
	FILE* f_quals = open_f(op_fp, "w");
	for(int i = 0; i < (int)distributions_per_read_posn->size(); i++)
	{
		vector<char>* cur_qual_sample = distributions_per_read_posn->at(i);

		double qual_sum = 0.0;
		for(int i_q = 0; i_q < (int)cur_qual_sample->size(); i_q++)
		{
			qual_sum += cur_qual_sample->at(i_q);
		} // i_q loop.

		double mean = qual_sum / (int)cur_qual_sample->size();
		double var = 0.0;
		for(int i_q = 0; i_q < (int)cur_qual_sample->size(); i_q++)
		{
			var += (cur_qual_sample->at(i_q) - mean) * (cur_qual_sample->at(i_q) - mean);
		} // i_q loop.

		var /= (cur_qual_sample->size()-1);

		fprintf(f_quals, "%d\t%lf\t%lf\n", i, mean, var);
	} // i loop.
	fclose(f_quals);
}

void load_mapped_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_ELAND(char* eland_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_tagAlign(char* tagalign_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_bowtie(char* bowtie_fp, vector<t_sequenced_read*>* sequenced_reads);
// ABOVE IS THE FASTQ LOADING INTERFACE:

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
	if(sscanf(cur_line, "%[^'\t'] %[^'\t'] %[^'\t'] %[^'\t'] %*[^'\t'] %[^'\t'] %*[^'\t'] %*[^'\t'] %*[^'\t'] %[^'\t'] %[^'\t']", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
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

//void preprocess_BED_read_line(char* cur_line, 
//	char* chrom, 
//	int& chr_index, int& sequenced_length, 
//	char& strand_char, 
//	char* mapping_quality_str)
//{
//	// Skip the comment and headers.
//	if(cur_line[0] == '@')
//	{
//		chrom[0] = 0;
//		chr_index = 0;
//		return;
//	}
//
//	int _chr_start;
//	int _chr_end;
//	char strand;
//
//	if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &_chr_start, &_chr_end, &strand) == 6)
//	{
//		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
//		int chr_start = translate_coord(_chr_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
//		int chr_end = translate_coord(_chr_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
//
//		// Check the flag and determine the strand.
//		strand_char = 'F';
//		if(strand_char == '-')
//		{
//			strand_char = 'R';
//		}
//
//		chr_index = chr_start;
//		sequenced_length = chr_end-chr_start+1;
//	}
//	else
//	{
//		chrom[0] = 0;
//	}
//}

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
	if(sscanf(cur_line, "%[^'\t'] %[^'\t'] %[^'\t'] %[^'\t'] %[^'\t']", read_id, strand_sign_str, chrom, chr_start_index_str, nucs) == 5)
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

void buffer_per_nucleotide_preprocessed_read_profile_no_buffer(char* sorted_read_fp,
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int max_l_read,
	int l_buffer, int& l_data)
{
	for(int i = 0; i <= l_buffer; i++)
	{
		if(signal_profile_buffer != NULL)
		{
			signal_profile_buffer[i] = 0.0;
		}

		// If the per strand signal generation is requested, initialize the per strand signal.
		if(forward_strand_signal != NULL)
		{
			forward_strand_signal[i] = 0.0;
			reverse_strand_signal[i] = 0.0;
		}
	} // i loop.

	char strand_char = 0;
	//char cur_fragment[10000];
	char mapping_map_str[10000];
	int chr_index;
	int n_processed_reads = 0;
	FILE* f_sorted_reads = open_f(sorted_read_fp, "r");
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

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		int read_start = 0;
		int read_end = 0;

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
				if(read_start == 0)
				{
					read_start = chr_index;
				}
			} // CIGAR string matching check

			read_end = chr_index + l_cur_entry - 1;

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

		// Increase the height: The indexing is already 1 based, there is no conversion needed.
		if(signal_profile_buffer != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			for(int i_cur = read_start; i_cur <= read_end; i_cur++)
			{
				signal_profile_buffer[i_cur]++;
			} // i_cur loop.
		}

		// Update the strand signals, if requested.
		if(forward_strand_signal != NULL &&
			max_l_read >= (read_end-read_start+1))
		{
			if(strand_char == 'F')
			{
				// Update the forward strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					forward_strand_signal[i_cur]++;
				} // i_cur loop.
			}
			else
			{
				// Update the reverse strand signal.
				for(int i_cur = read_start; i_cur <= read_end; i_cur++)
				{
					reverse_strand_signal[i_cur]++;
				} // i_cur loop.
			}
		} // strand signal check.

		delete(cur_entry_length_str);
		delete [] cur_read;
	} // file reading loop.
	fclose(f_sorted_reads);

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

		if(sscanf(cur_read, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_read);
			exit(0);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		// Parse the cigar string to get the fragments.
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		// Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
		bool splicing_tag_extension_check_pass = true;
		if(is_read_spliced && l_extended_tag > 0)
		{
			fprintf(stderr, "Spliced read with tag extension: %s\r", mapping_map_str);
			splicing_tag_extension_check_pass = false;
		}

		while(splicing_tag_extension_check_pass &&
			mapping_map_str_valid && 
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
					// Extend the read, update the profile.
					int ext_start = 0;
					int ext_end = 0;
					if(strand_char == 'F')
					{
						ext_start = chr_index;
					}
					else
					{
						ext_start = (chr_index + l_cur_entry - 1) - (l_extended_tag - 1);
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
					}
				} // extension length check.
			}

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

		delete(cur_entry_length_str);
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

void delete_fragments(vector<t_mapped_fragment*>* fragment_list)
{
	for(int i_frag = 0; i_frag < (int)fragment_list->size(); i_frag++)
	{
		delete(fragment_list->at(i_frag));
	}

	fragment_list->clear();
	delete(fragment_list);
}

void delete_fragments(t_mapped_fragment** fragment_list)
{
	int i_frag = 0;
        while(fragment_list[i_frag] != NULL)
        {
                free(fragment_list[i_frag]);
		i_frag++;
        }

        delete[] fragment_list;
}

//void load_fragments(char* mapped_reads_fp, 
//	vector<t_mapped_fragment*>* fore_strand_frags, vector<t_mapped_fragment*>* rev_strand_frags, 
//	int max_n_pcr_amplified_reads)
//{
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//	fprintf(stderr, "Loading fragments from %s.\n", mapped_reads_fp);
//}
//
//	vector<t_mapped_read*>* fore_reads = new vector<t_mapped_read*>();
//	vector<t_mapped_read*>* rev_reads = new vector<t_mapped_read*>();
//
//	// Add the fragments: If the fragmens do not exist, the algorithm returns.
//	//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
//	load_reads(mapped_reads_fp, fore_reads, rev_reads, max_n_pcr_amplified_reads);
//	get_mapped_fragments_per_mapped_reads(fore_reads, fore_strand_frags);
//	get_mapped_fragments_per_mapped_reads(rev_reads, rev_strand_frags);
//
//	// Delete the memory for the mapped reads, that is not needed any more.
//	for(int i_f_r = 0; i_f_r < (int)fore_reads->size(); i_f_r++)
//	{
//		delete_mapped_read(fore_reads->at(i_f_r));
//	} // i_f_r loop.
//	delete fore_reads;
//
//	for(int i_r_r = 0; i_r_r < (int)rev_reads->size(); i_r_r++)
//	{
//		delete_mapped_read(rev_reads->at(i_r_r));
//	} // i_f_r loop.
//	delete rev_reads;
//
//if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
//{
//	fprintf(stderr, "Loaded %ld forward and %ld reverse fragments.\n", fore_strand_frags->size(), rev_strand_frags->size());
//}
//}
//
//bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str)
//{
//        int mismatches = 0;
//        for(int i = 0; i < (int)strlen(fragment); i++)
//        {
//                char upper_subseq = toupper(chr_subseq[i]);
//                char upper_frag = toupper(fragment[i]);
//
//                if((upper_subseq != upper_frag) &&
//                        (upper_subseq != 'N') &&
//                        (upper_frag != 'N'))
//                {
//                        mismatches++;
//                        //printf("Fragments are not matching for %s at %d (%s). Enriched fragment is:\n%s\n%s\n", chr_fps->at(i_chr), chr_index, quality_str, cur_seq, chr_subseq);
//                        //getc(stdin);
//                        //exit(0);
//                }
//        } // i loop.
//
//        if(mismatches == 0)
//        {
//		return true;
///*
//                if(strcmp(quality_str, "U0") == 0)
//                {
//                        return(true);
//                }
//                else
//                {
//                        return(false);
//                }
//*/
//        }
//        else if(mismatches == 1)
//        {
//                if(strcmp(quality_str, "U1") == 0)
//                {
//                        return(true);
//                }
//                else
//                {
//                        return(false);
//                }
//        }
//        else if(mismatches == 2)
//        {
//                if(strcmp(quality_str, "U2") == 0)
//                {
//                        return(true);
//                }
//                else
//                {
//                        return(false);
//                }
//        }
//        else
//        {
//		return(false);
//        }
//}

bool sort_mapped_fragments_per_3p(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	int frag1_3p = fragment_3p_accessor(frag1);
	int frag2_3p = fragment_3p_accessor(frag2);

	return(frag1_3p < frag2_3p);
}

bool sort_mapped_fragments(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	if(frag1->base_index < frag2->base_index)
	{
		return(frag1->base_index < frag2->base_index);
	}
	else if(frag1->base_index == frag2->base_index)
	{
		if((frag1->base_index + frag1->sequenced_fragment_length) < (frag2->base_index + frag2->sequenced_fragment_length))
		{
			return(true);
		}
		else
		{
			return(false);
		}
	}
	else
	{
		return(false);
	}
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

// Count the statistics of mapped fragments for a given regions using binary search.
void get_read_statistics_per_region(vector<t_mapped_fragment*>* fragments, 
									vector<int>* all_fragment_3p_posns, 
									t_annot_region* region, 
									double& n_mapped_reads, 
									double& n_mapped_nucs)
{
	n_mapped_nucs = 0;
	n_mapped_reads = 0;

	// For the given region, count the reads with a binary search. Make sure that the fragments are sorted.
	int i_leftmost_read = locate_posn_per_sorted_posn_list(region->start, all_fragment_3p_posns, 0, all_fragment_3p_posns->size()-1);
	//int i_leftmost_read = locate_posn_per_sorted_obj_list(region->start, (vector<void*>*)all_fragment_3p_posns, 0, all_fragment_3p_posns->size()-1, fragment_3p_accessor);

	//fprintf(stderr, "Found leftmost read (%d-%d) @ %d. read for region %s:%d-%d\n", fragments->at(i_leftmost_read)->base_index, fragments->at(i_leftmost_read)->base_index + fragments->at(i_leftmost_read)->sequenced_fragment_length, i_leftmost_read, region->chrom, region->start, region->end);

	// Starting from the leftmost read, count the overlapping reads.
	while(i_leftmost_read < (int)fragments->size() && 
		fragments->at(i_leftmost_read)->base_index <= region->end)
	{
		// Get the overlap of the current read with the current region.
		int ol_start = MAX(fragments->at(i_leftmost_read)->base_index, region->start);
		int ol_end = MIN(fragments->at(i_leftmost_read)->base_index + fragments->at(i_leftmost_read)->sequenced_fragment_length-1, region->end);

		if(ol_end >= ol_start)
		{
			n_mapped_reads++;
			n_mapped_nucs += (ol_end - ol_start + 1);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
			fprintf(stderr, "%s:%d-%d: Read @ %d, Overlap: %d-%d Totals: %d, %d\n", region->chrom, region->start, region->end, fragments->at(i_leftmost_read)->base_index, ol_start, ol_end, (int)n_mapped_nucs, (int)n_mapped_reads);
		}

		i_leftmost_read++;
	} // i_leftmost_read loop.
}

//void load_reads_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
//	vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
//	int max_n_pcr_amplified_reads)
//{
//	char cur_dir_chr_ids_fp[1000];
//	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
//	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);
//
//	if(cur_dir_chr_ids == NULL)
//	{
//		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
//		getc(stdin);
//		return;
//	}
//
//	// Load the fragment in this chromosome from all the directories.
//	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
//	{
//		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//
//		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));
//
//		char cur_dir_cur_chr_reads_fp[1000];
//		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));
//
//		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
//		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();
//
//		// Add the fragments: If the fragmens do not exist, the algorithm returns.
//		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
//		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);
//
//		// Add the loaded fragments.
//		fore_strand_reads_per_chr->push_back(cur_chr_fore_reads);
//		rev_strand_reads_per_chr->push_back(cur_chr_rev_reads);
//	} // i_chr loop.
//}
//
//void load_fragments_per_dir(char* mapped_reads_dir, 
//	vector<char*>* chr_ids, 
//	vector<vector<t_mapped_fragment*>*>* fore_strand_fragments_per_chr, 
//	vector<vector<t_mapped_fragment*>*>* rev_strand_fragments_per_chr, 
//	int max_n_pcr_amplified_reads)
//{
//	char cur_dir_chr_ids_fp[1000];
//	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
//	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);
//
//	if(cur_dir_chr_ids == NULL)
//	{
//		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
//		getc(stdin);
//		return;
//	}
//
//	// Load the fragment in this chromosome from all the directories.
//	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
//	{
//		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//
//		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));
//
//		char cur_dir_cur_chr_reads_fp[1000];
//		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));
//
//		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
//		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();
//
//		// Add the fragments: If the fragmens do not exist, the algorithm returns.
//		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
//		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);
//
//		vector<t_mapped_fragment*>* cur_chr_fore_frags = new vector<t_mapped_fragment*>(); 
//		get_mapped_fragments_per_mapped_reads(cur_chr_fore_reads, cur_chr_fore_frags);
//
//		vector<t_mapped_fragment*>* cur_chr_rev_frags = new vector<t_mapped_fragment*>(); 
//		get_mapped_fragments_per_mapped_reads(cur_chr_rev_reads, cur_chr_rev_frags);
//
//		// Add the loaded fragments.
//		fore_strand_fragments_per_chr->push_back(cur_chr_fore_frags);
//		rev_strand_fragments_per_chr->push_back(cur_chr_rev_frags);
//
//		// Delete the memory for the mapped reads, that is not needed any more.
//		delete_mapped_reads(cur_chr_fore_reads);
//		delete_mapped_reads(cur_chr_rev_reads);
//	} // i_chr loop.
//}

vector<t_mapped_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_mapped_fragment*>* fore_frag_list, vector<t_mapped_fragment*>* rev_frag_list, int enrichment_mapped_fragment_length)
{
	vector<t_mapped_fragment*>* combined_frags = new vector<t_mapped_fragment*>();

	if(fore_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)fore_frag_list->size(); i_frag++)
		{
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));
			new_fragment_node->base_index = fore_frag_list->at(i_frag)->base_index;

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > fore_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				new_fragment_node->sequenced_fragment_length = fore_frag_list->at(i_frag)->sequenced_fragment_length;
			}


			new_fragment_node->strand_char = fore_frag_list->at(i_frag)->strand_char;

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// For the reverse fragment list, we need to also set the base index since it may be move further to left.
	if(rev_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)rev_frag_list->size(); i_frag++)
		{
			// If this enrichment fragment length is set to 0, this means that there is no fragment extension.
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > rev_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				// The enrichment fragment length is larger than the sequenced length, translate the fragment base.
				new_fragment_node->base_index = (rev_frag_list->at(i_frag)->base_index + rev_frag_list->at(i_frag)->sequenced_fragment_length - enrichment_mapped_fragment_length);
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				// The enrichment fragment length is smaller than the sequenced fragment length, do not do any translation of the fragment base.
				new_fragment_node->base_index = rev_frag_list->at(i_frag)->base_index;
				new_fragment_node->sequenced_fragment_length = rev_frag_list->at(i_frag)->sequenced_fragment_length;
			}

			new_fragment_node->strand_char = 'F'; // Change the place of these to forward strand.

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// Sort one last time.
	sort(combined_frags->begin(), combined_frags->end(), sort_mapped_fragments);

	return(combined_frags);
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

void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments)
{
	int i_mapp_map = 0;
	t_string* cur_entry_length_str = new t_string();
	bool is_matching = false;
	char entry_type_char;
	char strand_char = mapped_read->strand;
	//int chr_index = (mapped_read->strand=='F')?(mapped_read->base_index):(mapped_read->base_index-mapped_read->span+1);
	int chr_index = (mapped_read->base_index);
	char* mapping_map_str = mapped_read->mapping_str;

	//fprintf(stderr, "Processing cigar string: %s (%d, %c)\n", mapping_map_str, chr_index, strand_char);
	bool is_read_spliced = false;
	bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

	// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
	while(mapping_map_str_valid && 
		mapping_map_str[i_mapp_map] != 0)
	{
		int l_cur_entry = 0;
		get_next_entry_per_mapp_map_string(mapping_map_str,
											i_mapp_map, 
											is_matching,
											l_cur_entry,
											entry_type_char);

		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.		

		//int l_fragment = strlen(cur_fragment);
		if(is_matching)
		{
			//int l_fragment = get_l_fragment_per_cigar(quality_str);
			if(strand_char == 'F')
			{
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;
		
				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
			}
			else if(strand_char == 'R')
			{
				// Allocate and initialize a fragment and add it to the reverse strand fragment list.			
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				//new_fragment->base_index = chr_index + l_cur_entry - 1;
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;

				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
				//rev_strand_frags->push_back(new_fragment);
			} // reverse strand check.
		} // maching check.

		// Update the base for the current entry.
		// Must check whether to update the mapping posn: Update only for D and M entries.
		//if(entry_type_char == 'D' || 
		//	entry_type_char == 'M' ||
		//	entry_type_char == 'N' ||
		//	entry_type_char == 'H')
		if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
		{
			chr_index += l_cur_entry;
		}
	} // mapping map string processing loop.

	delete cur_entry_length_str;
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

void exclude_reads_per_regions_per_chr(char* read_chr_id, 
	vector<t_mapped_read*>* cur_chr_reads,
	vector<t_mapped_read*>* no_overlap_reads,
	vector<t_annot_region*>* regions_2_exclude)
{
	if(regions_2_exclude == NULL)
	{
		for(int i_r = 0; i_r < (int)cur_chr_reads->size(); i_r++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_r));
		} // i_r loop.

		return;
	}

	// Restructure the regions.
	t_sorted_annot_region_lists* restructured_regions = restructure_annot_regions(regions_2_exclude);
	
	int n_total_excluded_reads = 0;
	
	int n_excluded_reads = 0;

	int i_reg_chr = t_string::get_i_str(restructured_regions->chr_ids, read_chr_id);

	if(i_reg_chr < (int)restructured_regions->chr_ids->size())
	{
		fprintf(stderr, "Excluding the reads in %s\n", restructured_regions->chr_ids->at(i_reg_chr));			

		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			int cur_left_base_posn = (cur_chr_reads->at(i_read)->base_index);
			int cur_right_base_posn = (cur_chr_reads->at(i_read)->base_index + cur_chr_reads->at(i_read)->span - 1);

			int left_reg_i = locate_posn_per_sorted_obj_list(cur_left_base_posn, (vector<void*>*)(restructured_regions->regions_per_chrom[i_reg_chr]), 0, restructured_regions->regions_per_chrom[i_reg_chr]->size()-1, region_5p_accessor);

			bool cur_read_overlaps = false;
			while(left_reg_i > 0 &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end > cur_left_base_posn)
			{
				left_reg_i--;
			}

			while(left_reg_i < (int)restructured_regions->regions_per_chrom[i_reg_chr]->size() &&
				restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start < cur_right_base_posn)
			{
				int ol_start = MAX(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, cur_left_base_posn);
				int ol_end = MIN(restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end, cur_right_base_posn);

				if(ol_end >= ol_start)
				{
					cur_read_overlaps = true;
					break;
				}

				left_reg_i++;
			} // go over the region for the read.

			if(cur_read_overlaps)
			{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
				if(cur_chr_reads->at(i_read)->strand == 'R')
				{
					fprintf(stderr, "Read: %s:%d(%d)(%c) overlaps with %s:%d-%d\n", read_chr_id, cur_chr_reads->at(i_read)->base_index, cur_chr_reads->at(i_read)->span, cur_chr_reads->at(i_read)->strand, 
						restructured_regions->chr_ids->at(i_reg_chr), restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->start, restructured_regions->regions_per_chrom[i_reg_chr]->at(left_reg_i)->end);
					getc(stdin);
				}
}

				// Update # of excluded reads.
				n_excluded_reads++;
			}
			else
			{
				no_overlap_reads->push_back(cur_chr_reads->at(i_read));
			}
		} // i_read loop.
	} // i_reg_chr search check.
	else
	{
		fprintf(stderr, "Skipping the regions in %s\n", read_chr_id);
		for(int i_read = 0; i_read < (int)cur_chr_reads->size(); i_read++)
		{
			no_overlap_reads->push_back(cur_chr_reads->at(i_read));
		} // i_read loop.
	} // read pile chromosome check.

	fprintf(stderr, "Excluded %d reads.\n", n_excluded_reads);

	n_total_excluded_reads += n_excluded_reads;

	delete_restructured_annot_regions(restructured_regions);

	fprintf(stderr, "Excluded %d reads in total.\n", n_total_excluded_reads);
}

void exclude_reads_per_regions(vector<char*>* read_chr_ids, 
	vector<vector<t_mapped_read*>*>* reads_per_chrs, 
	vector<vector<t_mapped_read*>*>* no_overlap_reads_per_chrs,
	vector<t_annot_region*>* regions_2_exclude)
{
	for(int i_read_chr = 0; i_read_chr < (int)read_chr_ids->size(); i_read_chr++)
	{
		vector<t_mapped_read*>* cur_chr_no_overlap_reads = new vector<t_mapped_read*>();
		no_overlap_reads_per_chrs->push_back(cur_chr_no_overlap_reads);

		exclude_reads_per_regions_per_chr(read_chr_ids->at(i_read_chr), 
			reads_per_chrs->at(i_read_chr),
			cur_chr_no_overlap_reads,
			regions_2_exclude);
	} // i_read_chr loop.
}

// Merge all the fragments of all the reads.
void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		add_mapped_fragments_per_mapped_read(mapped_reads->at(i_r), mapped_fragments);
	} // i_r loop.
}

void delete_mapped_reads(vector<t_mapped_read*>* mapped_reads)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		delete_mapped_read(mapped_reads->at(i_r));
	} // i_r loop.

	delete(mapped_reads);
}

void delete_mapped_read(t_mapped_read* mapped_read)
{
	delete [] mapped_read->mapping_str;
	delete(mapped_read);
}

/*
Prune reads:
Note that the trick here is to deal with the reads that have exact same pattern of mapping. We do not care about the
strand, this should be taken care of before the function is called.

Note that the pruning must be done at the read level, not at the fragment level.
*/
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
	vector<t_mapped_read*>* pruned_forward_reads, 
	vector<t_mapped_read*>* pruned_reverse_reads)
{
	// If the pruning is not requested, return all the reads.
	if(n_max_reps_per_posn == 0)
	{
		fprintf(stderr, "Skipping pruning.\n");
		for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		} // i_r loop.

		return;
	}

	// Sort the mapped reads with respect to their 5' posn.
	sort(mapped_reads->begin(), mapped_reads->end(), sort_mapped_reads_per_5p);

    // First get rid of the extra fragments on forward strand.
	int* rep_cnts = new int[mapped_reads->size() + 2];
	memset(rep_cnts, 0, sizeof(int) * (mapped_reads->size() + 1));

	int prev_read_left_posn = 0;
    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		int cur_read_left_posn = (mapped_reads->at(i_r)->base_index);
		
		if(i_r > 0 &&
			prev_read_left_posn == cur_read_left_posn)
		{
			rep_cnts[i_r] = rep_cnts[i_r-1] + 1;
		}
		else // This is a new fragment set its copy number to 1.
		{
			rep_cnts[i_r] = 1;
		}

		// Update prev. left posn.
		prev_read_left_posn = cur_read_left_posn;
    } // i_frag loop.

    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		if(rep_cnts[i_r] <= n_max_reps_per_posn)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		}
		else // This is a new fragment set its copy number to 1.
		{
			// Delete the pruned reads, otherwise they will be lost.
			delete_mapped_read(mapped_reads->at(i_r));

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Pruning %d repetition read @ %d, %c.\n", rep_cnts[i_r], mapped_reads->at(i_r)->base_index, mapped_reads->at(i_r)->strand);
			getc(stdin);
}
		}
    } // i_r loop.

	delete [] rep_cnts;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	fprintf(stderr, "Pruned to %ld forward, %ld reverse strand fragments.\n", pruned_forward_reads->size(), pruned_reverse_reads->size());
}

bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2)
{
	int frag1_5p = read_5p_accessor(read1);
	int frag2_5p = read_5p_accessor(read2);

	return(frag1_5p < frag2_5p);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}

//int read_3p_accessor(void* obj_ptr)
//{
//	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;
//
//	if(frag_obj_ptr->strand == 'F')
//	{
//		return(frag_obj_ptr->base_index+frag_obj_ptr->span-1);
//	}
//	else if(frag_obj_ptr->strand == 'R')
//	{
//		return(frag_obj_ptr->base_index);
//	}
//	else
//	{
//		fprintf(stderr, "The strand char for fragment object is %s @ %s(%d)\n", frag_obj_ptr->strand, __FILE__, __LINE__);
//		exit(0);
//		return(-1);
//	}
//}

