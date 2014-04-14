#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>

#include "ms_ansi_cli.h"
#include "ms_gsl_fft_filter_utils.h"
#include "ms_mapped_read_tools.h"
#include "ms_profile_normalization.h"
#include "ms_utils.h"
#include "ms_signal_enrichment_utils.h"
#include "ms_peak_calling_utils.h"
#include "ms_signal_track_tools.h"
#include "ms_ansi_string.h"
#include "ms_utils.h"
#include "ms_annot_region_tools.h"
#include "ms_genomics_coords.h"
#include "ms_nomenclature.h"
#include "ms_xlog_math.h"
#include "ms_genome_sequence_tools.h"

bool __DUMP_PEAK_MESSAGES__ = false;

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
Read Preprocessing:\n\
	-preprocess [File format (\"SAM\"/\"eland\"/\"bowtie\"/\"BED\")] [Mapped reads file path (\"stdin\" for piped input)] [Output directory]\n\
	-sort_reads [Reads directory] [Output directory]\n\
	-remove_duplicates [Sorted reads directory] [Max # of duplicates per position] [Output directory]\n\
Peak Selection:\n\
	-get_multiscale_ERs [Options/Values]\n\
P-val window length processing:\n\
	-get_per_window_p_vals_vs_FC [Reads directory]\n", argv[0]);
		exit(0);
	}
	
	/************************************************************************
	TODO: Combine all the preprocessing (preprocess, sort, dedup) under 
	-preprocess option.
	***********************************************************************/
	if(strcmp(argv[1], "-preprocess") == 0)
	{
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -preprocess [File format (\"SAM\"/\"ELAND\"/\"bowtie\"/\"tagAlign\"/\"BED\")] [Mapped reads file path (\"stdin\" for piped input)] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* format_str = argv[2];
		char* mapped_reads_fp = argv[3];
		char* op_dir = argv[4];

		printf("Preprocessing:\ninput format: %s\nchip_seq_eland_op_fp: %s\nparsed_reads_op_dir: %s\n", 
			format_str, mapped_reads_fp, op_dir);

		// Do not do validation while dumping the binary read files.
		if(strcmp(format_str, "ELAND") == 0)
		{
			// Read the ELAND file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_ELAND_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_ELAND_read_line, false);
		}
		else if(strcmp(format_str, "SAM") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_SAM_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_SAM_read_line, false);
		}
		else if(strcmp(format_str, "tagAlign") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_tagAlign_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_tagAlign_read_line, false);
		}
		else if(strcmp(format_str, "bowtie") == 0)
		{
			// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_bowtie_read_line, false);
		}
		else if(strcmp(format_str, "BED") == 0)
		{
			// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(mapped_reads_fp, op_dir, preprocess_BED5_read_line, false);
		}
		else
		{
			printf("Unknown format for the mapped read file name, use SAM/ELAND/bowtie/tagAlign.\n");
			exit(0);
		}
	} // -preprocess option.
	else if(strcmp(argv[1], "-sort_reads") == 0)
	{
		if(argc != 4)
		{
			fprintf(stderr, "%s -sort_preprocessed_reads [Preprocessed read directory] [Sorted reads output directory]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];
		char* sorted_reads_op_dir = argv[3];

		int bucket_size = 50*1000*1000;

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);

		// Process per chromosome.
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);

		char sorted_chr_ids_fp[1000];
		sprintf(sorted_chr_ids_fp, "%s/chr_ids.txt", sorted_reads_op_dir);
		FILE* f_sorted_chr_ids = open_f(sorted_chr_ids_fp, "w");
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(f_sorted_chr_ids, "%s\n", chr_ids->at(i_chr));
		} // i_chr loop.
		fclose(f_sorted_chr_ids);

		// After sorting, dump the reads into a new directory named xxxx_sorted/.
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));
			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			// Load the entries, store them in a buffer.
			//FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");
			vector<FILE*>* bucket_f_list = new vector<FILE*>();
			vector<int>* bucket_starts = new vector<int>();
			vector<int>* bucket_sizes = new vector<int>();
			while(1)
			{
				// Get the read and determine which bucket it belongs to.
				char* cur_line = getline(f_cur_chr_reads);
				if(cur_line == NULL)
				{
					break;
				}

				char cur_cigar_str[1000];
				char cur_strand;
				int cur_start;
				sscanf(cur_line, "%s %c %d", cur_cigar_str, &cur_strand, &cur_start);
				int bucket_start = (int)(floor((double)cur_start / bucket_size) * bucket_size);

				bool bucket_found = false;
				for(int i_st = 0; i_st < (int)bucket_starts->size(); i_st++)
				{
					if(bucket_start == bucket_starts->at(i_st))
					{
						bucket_found = true;
						fprintf(bucket_f_list->at(i_st), "%s\n", cur_line);
						bucket_sizes->at(i_st)++;
					}
				} // i_st loop.

				// Need a new bucket?
				if(bucket_found == false)
				{
					bucket_starts->push_back(bucket_start);
					char new_bucket_fp[1000];
					//sprintf(new_bucket_fp, "%s/bucket_%d.txt", preprocessed_reads_dir, bucket_start);
					sprintf(new_bucket_fp, "%s/bucket_%d.txt", sorted_reads_op_dir, bucket_start);

if(__DUMP_PEAK_MESSAGES__)
					fprintf(stderr, "Adding bucket %s\n", new_bucket_fp);

					FILE* f_bucket = open_f(new_bucket_fp, "w");
					bucket_f_list->push_back(f_bucket);
					bucket_sizes->push_back(0);

					bucket_found = false;
					for(int i_st = 0; i_st < (int)bucket_starts->size(); i_st++)
					{
						if(bucket_start == bucket_starts->at(i_st))
						{
							bucket_found = true;	
							fprintf(bucket_f_list->at(i_st), "%s\n", cur_line);
							bucket_sizes->at(i_st)++;
						}
					} // i_st loop.				
				}

				delete [] cur_line;
			} // file reading loop.

			fclose(f_cur_chr_reads);

			for(int i_st = 0; i_st < (int)bucket_f_list->size(); i_st++)
			{
				fclose(bucket_f_list->at(i_st));
			} // i_st loop.
			delete(bucket_f_list);

			char cur_chr_sorted_reads_fp[1000];
			sprintf(cur_chr_sorted_reads_fp, "%s/%s_mapped_reads.txt", sorted_reads_op_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_sorted_reads = open_f(cur_chr_sorted_reads_fp, "w");
			sort(bucket_starts->begin(), bucket_starts->end());
			for(int i_buck = 0; i_buck < (int)bucket_starts->size(); i_buck++)
			{
				// Load the reads in the current bucket, sort them, dump them.
				char cur_bucket_fp[1000];
				//sprintf(cur_bucket_fp, "%s/bucket_%d.txt", preprocessed_reads_dir, bucket_starts->at(i_buck));
				sprintf(cur_bucket_fp, "%s/bucket_%d.txt", sorted_reads_op_dir, bucket_starts->at(i_buck));
				fprintf(stderr, "Sorting reads in bucket %s.\n", cur_bucket_fp);
				vector<char*>* cur_sorted_bucket_read_lines = sort_bucket_read_lines(cur_bucket_fp);

				for(int i_l = 0; i_l < (int)cur_sorted_bucket_read_lines->size(); i_l++)
				{
					fprintf(f_cur_chr_sorted_reads, "%s\n", cur_sorted_bucket_read_lines->at(i_l));

					// Free memory.	
					delete [] cur_sorted_bucket_read_lines->at(i_l);
				} // i_l loop.
				
				delete cur_sorted_bucket_read_lines;
			} // i_buck loop.
			fclose(f_cur_chr_sorted_reads);

			delete bucket_starts;
			delete bucket_sizes;
		} // i_chr loop.
	} // -sort_reads
	else if(strcmp(argv[1], "-remove_duplicates") == 0)
	{
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -remove_duplicates [Sorted preprocessed read directory] [Max # of duplicates per position] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* sorted_preprocessed_reads_dir = argv[2];
		int max_n_amp_reads = atoi(argv[3]);
		char* pruned_reads_dir = argv[4];

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", sorted_preprocessed_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);

		char pruned_chr_ids_fp[1000];
		sprintf(pruned_chr_ids_fp, "%s/chr_ids.txt", pruned_reads_dir);
		FILE* f_pruned_chr_ids = open_f(pruned_chr_ids_fp, "w");
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(f_pruned_chr_ids, "%s\n", chr_ids->at(i_chr));
		} // i_chr loop.
		fclose(f_pruned_chr_ids);

		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Pruning %s\n", chr_ids->at(i_chr));
			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", sorted_preprocessed_reads_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			char cur_chr_pruned_reads_fp[1000];
			sprintf(cur_chr_pruned_reads_fp, "%s/%s_mapped_reads.txt", pruned_reads_dir, chr_ids->at(i_chr));
			FILE* f_cur_chr_pruned_reads = open_f(cur_chr_pruned_reads_fp, "w");

			int n_processed_reads = 0;
			int n_pruned_reads = 0;
			int prev_read_start = 0;
			int n_amp_reads = 0;
			while(1)
			{
				char* cur_line = getline(f_cur_chr_reads);
				if(cur_line == NULL)
				{
					break;
				}

if(__DUMP_PEAK_MESSAGES__)
{
				if(n_processed_reads % 1000000 == 0)
				{
					fprintf(stderr, "Processing %d. read.              \r", n_processed_reads);
				}
}

				int cur_read_start = 0;
				if(sscanf(cur_line, "%*s %*s %d", &cur_read_start) != 1)
				{
					fprintf(stderr, "Could not parse: %s\n", cur_line);
					exit(0);
				}

				// Check if the read start updated.
				if(cur_read_start == prev_read_start)
				{
					n_amp_reads++;
				}
				else
				{
					prev_read_start = cur_read_start;
					n_amp_reads = 1;
				}

				if(n_amp_reads <= max_n_amp_reads)
				{
					n_pruned_reads++;
					fprintf(f_cur_chr_pruned_reads,  "%s\n", cur_line);
				}

				n_processed_reads++;
			} // file reading loop.
			fclose(f_cur_chr_reads);
			fclose(f_cur_chr_pruned_reads);

			fprintf(stderr, "Processed %d reads, pruned to %d reads.\n", n_processed_reads, n_pruned_reads);
		} // i_chr loop.
	} // -remove_duplicates
	else if(strcmp(argv[1], "-count_preprocessed_reads") == 0)
	{
		if(argc != 3)
		{
			fprintf(stderr, "USAGE: %s -count_preprocessed_reads [Preprocessed reads directory]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];

		// Load all the reads per chromosome: Load all the chromosomes.
		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.txt", preprocessed_reads_dir);

		// Load the chromosome id's.
		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<int>* n_reads_per_chromosome = new vector<int>();

		int n_total_reads = 0;
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char cur_chr_reads_fp[1000];
			sprintf(cur_chr_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			//fprintf(stderr, "Counting the reads in %s (%s)\n", chr_ids->at(i_chr), cur_chr_reads_fp);

			//int n_cur_chr_reads = 0;
			FILE* f_cur_chr_reads = open_f(cur_chr_reads_fp, "r");

			int n_cur_chr_reads = get_n_lines(f_cur_chr_reads);

			fclose(f_cur_chr_reads);

			n_reads_per_chromosome->push_back(n_cur_chr_reads);

			// Update total # of reads.
			n_total_reads += n_cur_chr_reads;

			//fprintf(stderr, "%s: %d reads. (%d reads in total)\n", chr_ids->at(i_chr), n_cur_chr_reads, n_total_reads);
			fprintf(stdout, "%s\t%d\n", chr_ids->at(i_chr), n_cur_chr_reads);
		} // i_chr loop.

		fprintf(stdout, "Total\t%d\n", n_total_reads);
	} 
	else if(strcmp(argv[1], "-estimate_fragment_length_per_mapability") == 0)
	{
		if(argc != 7)
		{
			fprintf(stderr, "USAGE: %s -estimate_fragment_length_per_mapability [Preprocessed reads directory] [Mapability signal directory] [Min separation] [Max separation] [Mapability read length]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];
		char* mapability_signal_dir = argv[3];
		int min_separation = atoi(argv[4]);
		int max_separation = atoi(argv[5]);
		int l_mapability_read = atoi(argv[6]);

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
			exit(0);
		}

		double* per_separ_merged_cnts = new double[max_separation - min_separation + 2];
		memset(per_separ_merged_cnts, 0, sizeof(double) * (max_separation - min_separation + 2));
		per_separ_merged_cnts -= min_separation;
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Estimating fragment length for %s\n", chr_ids->at(i_chr));

			char cur_chrom_map_signal_fp[1000];
			sprintf(cur_chrom_map_signal_fp, "%s/%s.bin", mapability_signal_dir, chr_ids->at(i_chr));
			int l_map_profile = 0;
			double* cur_map_signal = load_per_nucleotide_binary_profile(cur_chrom_map_signal_fp, l_map_profile);
			for(int i_map = 1; i_map < l_map_profile; i_map++)
			{
				cur_map_signal[i_map] /= 2*l_mapability_read;
			} // i_map loop.

			char cur_chrom_reads_fp[1000];
			sprintf(cur_chrom_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			// Get the forward and reverse profiles.
			int l_buffer = get_l_signal_per_reads(cur_chrom_reads_fp, 200);
			fprintf(stderr, "Buffer length is %d\n", l_buffer);
			int* fore_profile = new int[l_buffer+10];
			memset(fore_profile, 0, sizeof(int) * (l_buffer+2));
			int* rev_profile = new int[l_buffer+10];
			memset(rev_profile, 0, sizeof(int) * (l_buffer+2));

			// Open and load the mapped reads file, update the profiles.
			FILE* f_mapped_reads = open_f(cur_chrom_reads_fp, "r");
			while(1)
			{
				int chr_index;
				char cur_cigar[1000];
				char cur_strand;
				if(fscanf(f_mapped_reads, "%s %c %d", cur_cigar, &cur_strand, &chr_index) == 3)
				{
					int l_cur_entry;
					int i_mapp_map = 0;
					bool is_matching = 0;
					char entry_type_char;

					bool is_read_spliced = false;
					//bool mapping_map_str_valid = validate_mapping_map_str(cur_cigar, is_read_spliced);

					if(is_read_spliced)
					{
						// Skip the spliced reads for the analysis.
					}
					else
					{
						get_next_entry_per_mapp_map_string(cur_cigar,
															i_mapp_map, 
															is_matching,
															l_cur_entry,
															entry_type_char);


						// Make sure that the mapability around this read is useable.
						double total_mapability_sig = 0.0;
						for(int i = chr_index; i <= chr_index + l_cur_entry; i++)
						{
							total_mapability_sig += cur_map_signal[i];
						} // i loop.

						double avg_mapability_sig = total_mapability_sig / (l_cur_entry + 1);

						if(avg_mapability_sig < 1.05 &&
							is_matching == true)
						{
							if(cur_strand == 'R')
							{
								rev_profile[chr_index + l_cur_entry] = 1;
							}
							else if(cur_strand == 'F')
							{
								fore_profile[chr_index] = 1;
							}
							else
							{
								fprintf(stderr, "Invalid strand character: %c\n", cur_strand);
								exit(0);
							}
						}
					} // Spliced read check.
				}
				else
				{
					break;
				}
			} // file reading loop.
			fclose(f_mapped_reads);

			// Start pushing the forward strand on backward strand signal.
			//int* overlaps_per_separation = new int[max_separation - min_separation + 1];
			char per_separ_overlap_cnts_fp[1000];
			sprintf(per_separ_overlap_cnts_fp, "per_separation_counts_%s.txt", chr_ids->at(i_chr));
			FILE* f_per_separ_overlap = open_f(per_separ_overlap_cnts_fp, "w");
			for(int cur_sep = min_separation; cur_sep <= max_separation; cur_sep++)
			{
				fprintf(stderr, "Processing separation of %d nucleotides.               \r", cur_sep);
				double cur_separ_merged_cnt = 0;
				for(int i = 1; i <= l_buffer; i++)
				{
					if(i+cur_sep < l_buffer)
					{
						//overlaps_per_separation[cur_sep-min_separation] += (fore_profile[i+cur_sep] * rev_profile[i]);
						if(fore_profile[i] > 0 || rev_profile[i+cur_sep] > 0)
						{
							cur_separ_merged_cnt++;
						}
					}
				} // i loop.

				fprintf(f_per_separ_overlap, "%d\t%lf\n", cur_sep, cur_separ_merged_cnt);
				per_separ_merged_cnts[cur_sep] += cur_separ_merged_cnt;
			} // cur_sep loop.

			fclose(f_per_separ_overlap);

			delete [] fore_profile;
			delete [] rev_profile;
			delete [] cur_map_signal;
		} // i_chr loop.

		// Dump the aggregate merged counts.
		int estimated_l_frag = 0;
		double min_separation_merged_cnts = 1000*1000*1000;
		FILE* f_aggregate_merged_cnts = open_f("aggregate_per_separation_counts.txt", "w");
		for(int cur_sep = min_separation; cur_sep <= max_separation; cur_sep++)
		{
			fprintf(f_aggregate_merged_cnts, "%d\t%lf\n", cur_sep, per_separ_merged_cnts[cur_sep]);

			if(min_separation_merged_cnts > per_separ_merged_cnts[cur_sep])
			{
				min_separation_merged_cnts = per_separ_merged_cnts[cur_sep];
				estimated_l_frag = cur_sep;
			}
		} // cur_sep loop.
		fclose(f_aggregate_merged_cnts);

		fprintf(stderr, "Estimated fragment length is %d\n", estimated_l_frag);
	} // -estimate_fragment_length_per_mapability option
	else if(strcmp(argv[1], "-estimate_fragment_length") == 0)
	{
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -estimate_fragment_length [Preprocessed reads directory] [Min separation] [Max separation]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];
		int min_separation = atoi(argv[3]);
		int max_separation = atoi(argv[4]);

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
			exit(0);
		}

		double* per_separ_merged_cnts = new double[max_separation - min_separation + 2];
		memset(per_separ_merged_cnts, 0, sizeof(double) * (max_separation - min_separation + 2));
		per_separ_merged_cnts -= min_separation;
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Estimating fragment length for %s\n", chr_ids->at(i_chr));

			char cur_chrom_reads_fp[1000];
			sprintf(cur_chrom_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			// Get the forward and reverse profiles.
			int l_buffer = get_l_signal_per_reads(cur_chrom_reads_fp, 200);
			fprintf(stderr, "Buffer length is %d\n", l_buffer);
			int* fore_profile = new int[l_buffer+10];
			memset(fore_profile, 0, sizeof(int) * (l_buffer+2));
			int* rev_profile = new int[l_buffer+10];
			memset(rev_profile, 0, sizeof(int) * (l_buffer+2));

			// Open and load the mapped reads file, update the profiles.
			FILE* f_mapped_reads = open_f(cur_chrom_reads_fp, "r");
			while(1)
			{
				int chr_index;
				char cur_cigar[1000];
				char cur_strand;
				if(fscanf(f_mapped_reads, "%s %c %d", cur_cigar, &cur_strand, &chr_index) == 3)
				{
					int l_cur_entry;
					int i_mapp_map = 0;
					bool is_matching = 0;
					char entry_type_char;

					bool is_read_spliced = false;
					//bool mapping_map_str_valid = validate_mapping_map_str(cur_cigar, is_read_spliced);

					if(is_read_spliced)
					{
						// Skip the spliced reads for the analysis.
					}
					else
					{
						get_next_entry_per_mapp_map_string(cur_cigar,
															i_mapp_map, 
															is_matching,
															l_cur_entry,
															entry_type_char);


						if(is_matching == true)
						{
							if(cur_strand == 'R')
							{
								rev_profile[chr_index + l_cur_entry] = 1;
							}
							else if(cur_strand == 'F')
							{
								fore_profile[chr_index] = 1;
							}
							else
							{
								fprintf(stderr, "Invalid strand character: %c\n", cur_strand);
								exit(0);
							}
						}
					} // Spliced read check.
				}
				else
				{
					break;
				}
			} // file reading loop.
			fclose(f_mapped_reads);

			// Start pushing the forward strand on backward strand signal.
			//int* overlaps_per_separation = new int[max_separation - min_separation + 1];
			char per_separ_overlap_cnts_fp[1000];
			sprintf(per_separ_overlap_cnts_fp, "per_separation_counts_%s.txt", chr_ids->at(i_chr));
			FILE* f_per_separ_overlap = open_f(per_separ_overlap_cnts_fp, "w");
			for(int cur_sep = min_separation; cur_sep <= max_separation; cur_sep++)
			{
				fprintf(stderr, "Processing separation of %d nucleotides.               \r", cur_sep);
				double cur_separ_merged_cnt = 0;
				for(int i = 1; i <= l_buffer; i++)
				{
					if(i+cur_sep < l_buffer)
					{
						//overlaps_per_separation[cur_sep-min_separation] += (fore_profile[i+cur_sep] * rev_profile[i]);
						if(fore_profile[i] > 0 && rev_profile[i+cur_sep] > 0)
						{
							cur_separ_merged_cnt++;
						}
					}
				} // i loop.

				fprintf(f_per_separ_overlap, "%d\t%lf\n", cur_sep, cur_separ_merged_cnt);
				per_separ_merged_cnts[cur_sep] += cur_separ_merged_cnt;
			} // cur_sep loop.

			fclose(f_per_separ_overlap);

			delete [] fore_profile;
			delete [] rev_profile;
		} // i_chr loop.

		// Dump the aggregate merged counts.
		int estimated_l_frag = 0;
		double min_separation_merged_cnts = 1000*1000*1000;
		FILE* f_aggregate_merged_cnts = open_f("aggregate_per_separation_counts.txt", "w");
		for(int cur_sep = min_separation; cur_sep <= max_separation; cur_sep++)
		{
			fprintf(f_aggregate_merged_cnts, "%d\t%lf\n", cur_sep, per_separ_merged_cnts[cur_sep]);

			if(min_separation_merged_cnts > per_separ_merged_cnts[cur_sep])
			{
				min_separation_merged_cnts = per_separ_merged_cnts[cur_sep];
				estimated_l_frag = cur_sep;
			}
		} // cur_sep loop.
		fclose(f_aggregate_merged_cnts);

		fprintf(stderr, "Estimated fragment length is %d\n", estimated_l_frag);
	} // -estimate_fragment_length
	else if(strcmp(argv[1], "-cross_correlation") == 0)
	{
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -cross_correlation [Preprocessed reads directory] [Min separation] [Max separation]\n", argv[0]);
			exit(0);
		}

		char* preprocessed_reads_dir = argv[2];
		int min_separation = atoi(argv[3]);
		int max_separation = atoi(argv[4]);

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", preprocessed_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load chromosome id's from %s\n", chr_ids_fp);
			exit(0);
		}

		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Estimating fragment length for %s\n", chr_ids->at(i_chr));

			char cur_chrom_reads_fp[1000];
			sprintf(cur_chrom_reads_fp, "%s/%s_mapped_reads.txt", preprocessed_reads_dir, chr_ids->at(i_chr));

			// Get the forward and reverse profiles.
			int l_buffer = get_l_signal_per_reads(cur_chrom_reads_fp, 200);
			//fprintf(stderr, "Buffer length is %d\n", l_buffer);
			double* fore_profile = new double[l_buffer+10];
			memset(fore_profile, 0, sizeof(double) * (l_buffer+2));
			double* rev_profile = new double[l_buffer+10];
			memset(rev_profile, 0, sizeof(double) * (l_buffer+2));

			int l_data;
			buffer_per_nucleotide_profile_no_buffer(cur_chrom_reads_fp, 1, NULL, fore_profile, rev_profile, l_buffer, l_data);

			char per_separ_corrs_fp[1000];
			sprintf(per_separ_corrs_fp, "per_separation_correlations_%s.txt", chr_ids->at(i_chr));
			FILE* f_per_separ_corrs = open_f(per_separ_corrs_fp, "w");
			for(int cur_sep = min_separation; cur_sep <= max_separation; cur_sep++)
			{
				fprintf(stderr, "Processing separation %d nucleotides.            \r", cur_sep);

				double cur_win_corr = 0.0;
				//get_correlation(&(rev_profile[cur_sep]), fore_profile, l_data - cur_sep, cur_win_corr);
				for(int i = 1; i <= l_buffer; i++)
				{
					if(i+cur_sep < l_buffer)
					{
						if(fore_profile[i] > 0 || rev_profile[i+cur_sep] > 0)
						{
							cur_win_corr += fore_profile[i] * rev_profile[i+cur_sep];
						}
					}
				} // i loop.

				fprintf(f_per_separ_corrs, "%d\t%lf\n", cur_sep, cur_win_corr);
			} // cur_sep loop.

			fclose(f_per_separ_corrs);

			delete [] fore_profile;
			delete [] rev_profile;
		} // i_chr loop.
	} // -cross_correlation.
	else if(strcmp(argv[1], "-write_BedGraph_per_preprocessed_reads") == 0)
	{
	} // -write_BedGraph_per_preprocessed_reads
	else if(strcmp(argv[1], "-get_multiscale_ERs") == 0)
	{
		// We need at least 5 parameters.
		if(argc < 5)
		{
			fprintf(stderr, "USAGE: %s -get_multiscale_ERs [Option/Values]\n\
							-chip [ChIP reads directory]\n\
							-control [control reads directory]\n\
							-mapp [multi-mapability profiles directory]\n\
							-begin_l [First scale smoothing window length (1000)]\n\
							-end_l [Last scale smoothing window length (16000)]\n\
							-step [Multiplicative window length step (1.5)]\n\
							-l_mapp [Read length of multi-mapability profiles]\n\
							-mapp_thr [Multi-mapability signal threshold used in correction (1.2)]\n\
							-l_frag [Fragment length (200)]\n\
							-l_c [Mapability correction window length (2000)]\n\
							-gamma [Min threshold for unsmoothed/smoothed (4)]\n", argv[0]);
			exit(0);
		}

		double base_scale_l_win = 1000;
		double end_scale_l_win = 16000; 
		double log_step = 1.5;
		int l_p_val_norm_win = 1750;
		int l_mapability_filtering_win = 2000;
		double max_normalized_mapability_signal = 1.2;
		double filtered_vs_non_filtered_max_scaler = 4;
		double signal_profile_scaler = 1;
		int l_fragment = 200;

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "-");
		bool ret = false;
		char* chip_reads_dir = cli->get_value_by_option("-chip", ret);
		if(!ret)
		{
			fprintf(stderr, "Need chip.\n");
			exit(0);
		}

		ret = false;
		char* control_reads_dir = cli->get_value_by_option("-control", ret);
		if(!ret)
		{
			fprintf(stderr, "Need control.\n");
			exit(0);
		}

		ret = false;
		char* mapability_signal_dir = cli->get_value_by_option("-mapp", ret);
		if(!ret)
		{
			mapability_signal_dir = NULL;
		}

		ret = false;
		char* l_read_mapability_signal_str = cli->get_value_by_option("-l_mapp", ret);
		int l_read_mapability_signal = 0;
		if(ret)
		{
			l_read_mapability_signal = atoi(l_read_mapability_signal_str);
		}
		else if(mapability_signal_dir != NULL)
		{
			fprintf(stderr, "Could not read the mapability read length.\n");
			exit(0);
		}

		ret = false;
		char* l_frag_str = cli->get_value_by_option("-l_frag", ret);
		if(ret)
		{
			l_fragment = atoi(l_frag_str);
		}

		ret = false;
		char* base_scale_l_win_str = cli->get_value_by_option("-begin_l", ret);
		if(ret)
		{
			base_scale_l_win = atof(base_scale_l_win_str);
		}

		ret = false;
		char* end_scale_l_win_str = cli->get_value_by_option("-end_l", ret);
		if(ret)
		{
			end_scale_l_win = atof(end_scale_l_win_str);
		}

		ret = false;
		char* log_step_str = cli->get_value_by_option("-step", ret);
		if(ret)
		{
			log_step = atof(log_step_str);
		}

		ret = false;
		char* l_p_val_norm_win_str = cli->get_value_by_option("-l_p", ret);
		if(ret)
		{
			l_p_val_norm_win = atoi(l_p_val_norm_win_str);
		}

		// Mapability correction window length.
		ret = false;
		char* l_mapability_filtering_win_str = cli->get_value_by_option("-l_c", ret);
		if(ret)
		{
			l_mapability_filtering_win = atoi(l_mapability_filtering_win_str);
		}

		ret = false;
		char* max_normalized_mapability_signal_str = cli->get_value_by_option("-mapp_thr", ret);
		if(ret)
		{
			max_normalized_mapability_signal = atof(max_normalized_mapability_signal_str);
		}

		ret = false;
		char* filtered_vs_non_filtered_max_scaler_str = cli->get_value_by_option("-gamma", ret);
		if(ret)
		{
			filtered_vs_non_filtered_max_scaler = atof(filtered_vs_non_filtered_max_scaler_str);
		}

		get_peaks(chip_reads_dir,
			control_reads_dir,
			l_fragment,
			mapability_signal_dir,
			l_read_mapability_signal,
			base_scale_l_win,
			end_scale_l_win,
			log_step,
			l_p_val_norm_win,
			l_mapability_filtering_win,
			max_normalized_mapability_signal,
			filtered_vs_non_filtered_max_scaler,
			signal_profile_scaler,    // Signal profile scaler.
			0.05,   // RD based end pruning p-value.
			.5,	// Forward/Reverse signal evenness.
			P_VAL_PER_WIN_MEAN_SIGNAL, // p-value selection.
			false,  // do_replace_profiles_w_smoothed_profiles
			false,  // do_BJ_correction_on_minima
			true,   // do_filter_minima_per_length_min_base_scale_l_win: For broad peaks, this is set to true to enforce that there are no small weird peaks.
			false,	// do_post_peak_p_value_minimization (p-val min)
			.5); // The flip probability.
	} // -get_multi_scale_ERs
	else if(strcmp(argv[1], "-get_per_window_p_vals_vs_FC") == 0)
	{
		// This option tries to find the p-value normalization window length that has the smallest FDR in comparison with the poisson based thresholding.
		// -analyze_p_value_normalization_window_lengths\n
		if(argc != 5)
		{
			fprintf(stderr, "%s -get_per_window_p_vals_vs_FC [ChIP reads directory] [Control reads directory] [Scaling factor (fragment length)]\n", argv[0]);
			exit(0);
		}

		char* chip_reads_dir = argv[2];
		char* control_reads_dir = argv[3];
		int l_fragment = atoi(argv[4]);

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", chip_reads_dir);
		vector<char*>* chr_ids = buffer_file(chr_ids_fp);

		// Go over all the chromosomes.
		//vector<t_per_win_stat_info*>* all_window_stat_info = new vector<t_per_win_stat_info*>();

		double* log_factorials = buffer_log_factorials(100*1000);

		FILE* f_win_stats = open_f("window_stats.txt", "w");
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Processing chromosome %s\n", chr_ids->at(i_chr));

			// Generate the current signal profile.
			int l_buffer = 300*1000*1000;
			int l_profile = 0;
			double* buffered_signal_profile = new double[l_buffer + 2];	
			char cur_chr_chip_reads_fp[1000];
			sprintf(cur_chr_chip_reads_fp, "%s/%s_mapped_reads.txt", chip_reads_dir, chr_ids->at(i_chr));
			buffer_per_nucleotide_profile_no_buffer(cur_chr_chip_reads_fp, l_fragment, 
				buffered_signal_profile, NULL, NULL,
				l_buffer, l_profile);

			double* signal_profile = new double[l_profile + 2];
			for(int i = 1; i <= l_profile; i++)
			{
				signal_profile[i] = buffered_signal_profile[i];
			} // i loop.
			delete [] buffered_signal_profile;

			// Generate the current control profile.
			l_buffer = 300*1000*1000;
			int l_control = 0;
			double* buffered_control_profile = new double[l_buffer + 2];
			char cur_chr_control_reads_fp[1000];
			sprintf(cur_chr_control_reads_fp, "%s/%s_mapped_reads.txt", control_reads_dir, chr_ids->at(i_chr));
			buffer_per_nucleotide_profile_no_buffer(cur_chr_control_reads_fp, l_fragment, 
				buffered_control_profile, NULL, NULL, 
				l_buffer, l_control);

			double* control_profile = new double[l_control + 2];
			for(int i = 1; i <= l_control; i++)
			{
				control_profile[i] = buffered_control_profile[i];
			} // i loop.
			delete [] buffered_control_profile;

			// Scale the control with respect to signal.
			double per_win_2DOF_lls_scaling_factor = 0;
			double per_win_1DOF_lls_scaling_factor = 0;
			double total_sig_scaling_factor = 0;
			get_per_window_scaling_factors_for_profile1_per_profile2(control_profile, l_control, signal_profile, l_profile, 10000,
																		per_win_2DOF_lls_scaling_factor,
																		per_win_1DOF_lls_scaling_factor,
																		total_sig_scaling_factor);

			double scaling_factor = per_win_1DOF_lls_scaling_factor;
			fprintf(stderr, "Scaling control with factor of %lf\n", scaling_factor);

			// Go over the whole control signal.
			for(int i = 0; i < l_control; i++)
			{
				control_profile[i] *= scaling_factor;
			} // i loop.

			for(int l_win = 500; l_win <= 10000; l_win += 250)
			{
				int cur_win_start = 1;
				while(cur_win_start < l_profile)
				{
					fprintf(stderr, "Processing l_win = %d (%d, %d)             \r", l_win, cur_win_start, l_profile);
					int cur_l_win = l_win;
					if(cur_win_start+cur_l_win > l_profile)
					{
						cur_l_win = l_profile - cur_win_start - 1;
					}

					// Over the current thresholding window, process all the p-value windows and count the false positives.
					// Using the current window size, compute the estimate of false positives.
					//t_per_win_stat_info* new_win = new t_per_win_stat_info();
					//all_window_stat_info->push_back(new_win);
					//new_win->chr_id = i_chr;
					int win_start = cur_win_start;
					int win_end = win_start + cur_l_win;

					// Check if there is good coverage of high signal in the current window.
					double total_signal = 0.0;
					double total_control = 0.0;
					for(int i = win_start; i < win_end; i++)
					{
						if(i < l_profile)
						{
							total_signal += signal_profile[i];
						}

						if(i < l_control)
						{
							total_control += control_profile[i];
						}
					} // i loop.
						
					// Make sure that the region is in the profile length limits.
					double log_p_val = 0;
					if(win_end < l_profile && 
						win_end < l_control)
					{
						log_p_val = get_binomial_pvalue_per_region(signal_profile, control_profile, win_start, win_end, log_factorials, l_fragment, true, cur_l_win);
					}
					else
					{
						log_p_val = 0.0;
					}

					if(log_p_val != 0 &&
						total_signal / (win_end-win_start) > 1)
					{
						fprintf(f_win_stats, "%d\t%d\t%d\t%lf\t%lf\t%lf\n", 
							i_chr, 
							win_start, 
							win_end, 
							log_p_val,
							total_signal,
							total_control);
					}

					// Update the thresholding window start.
					cur_win_start += l_win;
				} // cur_win_start loop.
			} // cur_l_win loop.	

			delete [] signal_profile;
			delete [] control_profile;
		} // i_chr loop.
		fclose(f_win_stats);
	} // -get_per_win_p_vals_FC option.
	else if(strcmp(argv[1], "-preprocess_FASTA") == 0)
	{
		if(argc != 4)
		{
			fprintf(stderr, "USAGE: %s -preprocess_FASTA [FASTA file path] [Output firectory]\n", argv[0]);
			exit(0);
		}

		char* fasta_fp = argv[2];
		char* bin_dir = argv[3];

		printf("Binarizing %s.\n", fasta_fp);
		FILE* f_fasta = open_f(fasta_fp, "r");

		char* cur_line = NULL;
		char* cur_entry_buffer = new char[250 * 1000 * 1000];
		int cur_entry_i = 0;
		char cur_entry_id[1000];

		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", bin_dir);
		FILE* f_chr_ids = open_f(chr_ids_fp, "w");

		while(1)
		{
			// Process the current buffer.
			cur_line = getline(f_fasta);

			if(cur_line == NULL)
			{
				// File ended, dump the last entry if there are values in it.
				if(cur_entry_i > 1)
				{
					char cur_entry_bin_fp[1000];
					normalize_chr_id(cur_entry_id);
					sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);

					fprintf(f_chr_ids, "%s\n", cur_entry_id);

					FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

					// Open the new binary file.
					fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

					// Dump the sequence buffer.
					fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
					fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

					// Close current file.
					fclose(f_bin);
				}
				break;
			}
			else if(cur_line[0] == '>')
			{
				if(cur_entry_i > 1)
				{
					char cur_entry_bin_fp[1000];
					normalize_chr_id(cur_entry_id);
					sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);

					fprintf(f_chr_ids, "%s\n", cur_entry_id);

					FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

					// Open the new binary file.
					fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

					// Dump the sequence buffer.
					fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
					fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

					// Close current file.
					fclose(f_bin);
				}

				// Update the id, reset the counter.
				strcpy(cur_entry_id, &(cur_line[1]));
				cur_entry_i = 0;
			}
			else
			{
				// Concatenate the current sequence line.
				int l_cur_line = t_string::string_length(cur_line);
				for(int i = 0; i < l_cur_line; i++)
				{
					cur_entry_buffer[cur_entry_i] = cur_line[i];
					cur_entry_i++;
				} // i loop.
			}

			delete [] cur_line;
		} // file reading loop.

		fclose(f_fasta);
		fclose(f_chr_ids);
		delete [] cur_entry_buffer;
	} // -preprocess_FASTA option.	
	else if(strcmp(argv[1], "-fragment_sequence_2_stdout") == 0)
	{
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -fragment_sequence_2_stdout [Binarized sequence file path] [Fragment length] [Output prefix]\n", argv[0]);
			exit(0);
		}

		char* binarized_seq_file_path = argv[2];
		int l_frag = atoi(argv[3]);
		char* op_prefix = argv[4];

		fprintf(stderr, "Fragmenting %s, using fragment length of %d to stdout.\n", 
			binarized_seq_file_path, l_frag);

		int l_seq = 0;
		char* bin_seq_data = load_binary_sequence_file(binarized_seq_file_path, l_seq);
		fprintf(stderr, "Loaded %d nucleotides.\n", l_seq);

		//int cur_i_bucket = 0;
		//int n_reads_per_cur_bucket = 0;
		//sprintf(cur_op_fp, "%s_%d.fasta", op_prefix, cur_i_bucket);
		//FILE* f_op = open_f(cur_op_fp, "w");
		char* cur_frag = new char[l_frag+2];
		for(int frag_start_i = 0; frag_start_i < l_seq-l_frag; frag_start_i++)
		{
			if(frag_start_i % (1000*1000) == 0)
			{
				fprintf(stderr, "Start: %d                  \r", frag_start_i);
			}

			// Get the current fragment.
			for(int i_f = frag_start_i; i_f < frag_start_i+l_frag; i_f++)
			{
				cur_frag[i_f-frag_start_i] = (char)(bin_seq_data[i_f]);
			} // i_f loop.

			// Dump the fragment in fastq format.
			fprintf(stdout, ">%s_F_%d\n%s\n", op_prefix, frag_start_i, cur_frag);
		} // frag_start_i loop.

		fprintf(stderr, "Dumping the reverse strand fragments.\n");

		// Dump the reverse complement reads.
		char* cur_reverse_comp = get_reverse_complement(bin_seq_data, 0, l_seq-1);
		delete [] bin_seq_data;
		bin_seq_data = cur_reverse_comp;
		for(int frag_start_i = 0; frag_start_i < l_seq-l_frag; frag_start_i++)
		{
			if(frag_start_i % (1000*1000) == 0)
			{
				fprintf(stderr, "Start: %d                  \r", frag_start_i);
			}

			// Get the current fragment.
			for(int i_f = frag_start_i; i_f < frag_start_i+l_frag; i_f++)
			{
				cur_frag[i_f-frag_start_i] = (char)(bin_seq_data[i_f]);
			} // i_f loop.

			// Dump the fragment in fastq format.
			fprintf(stdout, ">%s_R_%d\n%s\n", op_prefix, frag_start_i, cur_frag);
		} // frag_start_i loop.
	} // -fragment_sequence_2_stdout option.
	else if(strcmp(argv[1], "-get_multimapability_signal_per_mapped_reads") == 0)
	{
		// At this point, 
		if(argc != 5)
		{
			fprintf(stderr, "USAGE: %s -get_multimapability_signal_per_mapped_reads [Mapped reads file path (stdin)] [Output file path] [Read length]\n", argv[0]);
			exit(0);
		}

		char* mapped_reads_fp = argv[2];
		char* op_fp = argv[3];
		int l_read = atoi(argv[4]);

		int init_l = 300*1000*1000;
		double* signal_profile_buffer = new double[init_l+1];
		for(int i = 0; i <= init_l; i++)
		{
			signal_profile_buffer[i] = 0.0;
		} // i loop.

		// Buffer the current profile: Do not impose the extension in this step.
		int l_data = 0;
		buffer_per_nucleotide_profile_no_buffer(mapped_reads_fp, 0, signal_profile_buffer, NULL, NULL, init_l, l_data);

		// Normalize  the profile.
		for(int i = 1; i <= l_data; i++)
		{
			signal_profile_buffer[i] /= 2*l_read;
		} // i loop.

		// Convert the profile to char, multiply by 100.
		unsigned char* char_signal_buffer = new unsigned char[l_data + 2];
		for(int i = 0; i <= l_data; i++)
		{
			// Make sure that 
			if(ceil(signal_profile_buffer[i] * 100) > 255)
			{
				char_signal_buffer[i] = 255;
			}
			else
			{
				char_signal_buffer[i] = (unsigned char)(ceil(signal_profile_buffer[i] * 100));
			}
		} // i loop.

		// Dump the profile.
		dump_per_nucleotide_uchar_binary_profile(char_signal_buffer, l_data, op_fp);

		// Free memory.
		delete [] signal_profile_buffer;
		delete [] char_signal_buffer;
	} // -get_multimapability_signal_per_mapped_reads

	/*
Cross-correlation/Fragment length estimation:\n\
	-estimate_fragment_length [Reads directory] [Minimum separation] [Maximum separation]\n\
	-estimate_fragment_length_per_mapability [Reads directory] [Multi-mapability signal directory] [Minimum separation] [Maximum separation] [Mapability read length]\n\
	-cross_correlation [Reads directory] [Minimum separation] [Maximum separation]\n\
Signal generation:\n\
	-write_BedGraph_per_preprocessed_reads [Reads directory] [Fragment length]\n\
Multi-Mapability Signal Generation:\n\
	-preprocess_FASTA [FASTA file path] [Output firectory]\n\
	-fragment_sequence_2_stdout [Binarized sequence file path] [Fragment length]\n\
	-get_multimapability_signal_per_mapped_reads [Mapped reads file path (stdin)] [Output file path]\n
	*/
}

