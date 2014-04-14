#include "ms_signal_track_tools.h"
#include "ms_annot_region_tools.h"
#include "ms_genomics_coords.h"
#include "ms_xlog_math.h"
#include "ms_utils.h"
#include "ms_ansi_string.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

bool __DUMP_SIGNAL_TRACK_MSGS__ = false;

void get_profile_extrema(double* profile, int l_profile, double& prof_min, double& prof_max)
{
	prof_min = 1000*1000;
	prof_max = -1000*1000;
	for(int i = 1; i <= l_profile; i++)
	{
		if(profile[i] > prof_max)
		{
			prof_max = profile[i];
		}

		if(profile[i] < prof_min)
		{
			prof_min = profile[i];
		}
	} // i loop.
}

double* get_zero_indexed_per_one_indexed_data(double* one_indexed_data, int l_profile)
{
	double* zero_indexed_data = new double[l_profile];
	for(int i = 1; i <= l_profile; i++)
	{
		zero_indexed_data[i-1] = one_indexed_data[i];
	} // i loop.

	return(zero_indexed_data);
}

double* get_one_indexed_per_zero_indexed_data(double* zero_indexed_data, int l_profile)
{
	double* one_indexed_data = new double[l_profile+2];
	for(int i = 0; i < l_profile; i++)
	{
		one_indexed_data[i+1] = zero_indexed_data[i];
	} // i loop.

	return(one_indexed_data);
}

double* copy_profile(double* signal_profile, int l_profile)
{
	double* copy_prof = new double[l_profile+2];

	for(int i_sig = 1; i_sig <= l_profile; i_sig++)
	{
		copy_prof[i_sig] = signal_profile[i_sig];
	} // i_sig loop.

	return(copy_prof);
}

double* extract_one_indexed_profile_per_profile(double* signal_profile_buffer, int l_profile, int start, int end, int& l_extracted_profile)
{
	double* extracted_prof = new double[end - start + 3];	
	for(int i_sig = start; i_sig <= end; i_sig++)
	{
		extracted_prof[i_sig - start] = 0;
	} // i_sig loop.
	
	// Count the extracted profile.
	l_extracted_profile = 1;
	for(int i_sig = start; i_sig <= MIN(l_profile, end); i_sig++)
	{
		extracted_prof[l_extracted_profile] = signal_profile_buffer[i_sig];
		l_extracted_profile++;
	} // i_sig loop.

	// Profiles are one based.
	l_extracted_profile--;

	return(extracted_prof);
}

// The signal profile is 1 based, consistent with the codebase indexing.
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, char* chrom, char* op_fp)
{
	FILE* f_op = NULL;
	if(check_file(op_fp))
	{
		fprintf(stderr, "%s exists, concatting.\n", op_fp);
		f_op = open_f(op_fp, "a");
	}
	else
	{
		f_op = open_f(op_fp, "w");
	}

	// Get the bedgraph for the current profile.
	int i_nuc = 1; 
	double cur_height = signal_profile_buffer[i_nuc];
	int cur_block_start_i = i_nuc;
	i_nuc++;
	while(1)
	{
		// Find the point where the height changes: The end goes till it's equal to the profile length since the profile is 1-based.
		while(i_nuc <= l_profile)
		{
			// Wait till there is a change in the height, which marks the start of a new block.
			if(cur_height != signal_profile_buffer[i_nuc])
			{
				break;
			}

			i_nuc++;
		} // i_nuc loop.

		// At this point, either this is the end of the profile, or there was a change in the height, either way, this was the end of the current block. Definitely dump it.
		if(cur_height != signal_profile_buffer[i_nuc])
		{
			// Dump the current block.
			fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
				translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
				cur_height);

			// Update the new height and new start.
			cur_height = signal_profile_buffer[i_nuc];
			
			// Current position starts the next block.
			cur_block_start_i = i_nuc; 
		}

		// If the above block end was the end of the whole profile, we are done, otherwise continue to the next block.
		if(i_nuc > l_profile)
		{
			break;
		}

		//i_nuc++;
	} // i_nuc loop.

	//fprintf(f_op, "%s\t%d\t%d\t%lf\n", chrom, 
	//			translate_coord(cur_block_start_i, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//			translate_coord(i_nuc-1, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
	//			cur_height);

	fclose(f_op);
}

void dump_per_nucleotide_binary_profile_per_bedgraph(char* bgr_fp, bool dump_binary, char* op_fp)
{
	FILE* f_bgr = NULL;
	if(t_string::compare_strings(bgr_fp, "stdin"))
	{
		f_bgr = stdin;
	}
	else
	{
		f_bgr = open_f(bgr_fp, "r");
	}

	// Initialize the signal profile.
	int l_max = 300*1000*1000;
	double* signal_profile = new double[l_max+1];
	for(int i_nuc = 0; i_nuc <= l_max; i_nuc++)
	{
		signal_profile[i_nuc] = 0.0;
	} // i_nuc loop.

	// Go over all the input file and process all the lines.
	fprintf(stderr, "Dumping the profile to %s.\n", op_fp);
	while(1)
	{
		char* cur_bgr_line = getline(f_bgr);
		if(cur_bgr_line == NULL)
		{
			break;
		}

		char cur_chrom[1000];
		int start;
		int end;
		double cur_sig;
		if(sscanf(cur_bgr_line, "%s %d %d %lf", cur_chrom, &start, &end, &cur_sig) != 4)
		{
			fprintf(stderr, "Could not parse bgr file line: %s\n", cur_bgr_line);
			exit(0);
		}

		int trans_start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		int trans_end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		for(int i_nuc = trans_start; i_nuc <= trans_end; i_nuc++)
		{
			if(i_nuc > l_max)
			{
				fprintf(stderr, "Cannot set the position, greater than l_max: %d (%s)\n", i_nuc, cur_bgr_line);
				exit(0);
			}
			signal_profile[i_nuc] = cur_sig;
		} // i_nuc loop.

		delete [] cur_bgr_line;
	} // file reading loop.	

	// Close the bedgraph file.
	if(!t_string::compare_strings(bgr_fp, "stdin"))
	{
		fclose(f_bgr);
	}

	// Get the end of the signal.
	int l_data = l_max;
	while(signal_profile[l_data] == 0.0)
	{
		l_data--;
	} // i_nuc loop.

	fprintf(stderr, "Signal length is %d, dumping the per nucleotide profile.\n", l_data);

	if(dump_binary)
	{
		dump_per_nucleotide_binary_profile(signal_profile, l_data, op_fp);
		return;
	}

	// Dump the per nucleotide signal profile.
	FILE* f_op = NULL;
	
	fprintf(stderr, "Dumping ASCII.\n");
	f_op = open_f(op_fp, "w");
	
	fprintf(f_op, "%d\n",  l_data);
	
	for(int i_nuc = 1; i_nuc <= l_data; i_nuc++)
	{
		if(dump_binary)
		{
			fwrite(&(signal_profile[i_nuc]), sizeof(double), 1, f_op);
		}
		else
		{
			fprintf(f_op, "%lf ",  signal_profile[i_nuc]);
		}
	} // i_nuc loop.

	fclose(f_op);
}

void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(unsigned char), l_profile+1, f_op);

	fclose(f_op);
}

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	unsigned char* signal_profile_buffer = new unsigned char[l_profile+2];

if(__DUMP_SIGNAL_TRACK_MSGS__)
	fprintf(stderr, "Loading %d data values.\n", l_profile);

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(char), l_profile+1, f_prof);
	
	fclose(f_prof);

	return(signal_profile_buffer);
}

void dump_per_nucleotide_binary_profile(double* signal_profile, int l_profile, char* op_fp)
{
	// Dump the per nucleotide signal profile.
	FILE* f_op = open_f(op_fp, "wb");

	// Write data length to first couple bytes.
	fwrite(&(l_profile), sizeof(int), 1, f_op);

	// Dump the data: Dump 0 based data.
	fwrite(&(signal_profile[1]), sizeof(double), l_profile+1, f_op);

	fclose(f_op);
}

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile)
{
	FILE* f_prof = open_f(binary_per_nucleotide_profile_fp, "rb");

	// Read the profile length.
	int l_data = 0;
	fread(&l_data, sizeof(int), 1, f_prof);
	l_profile = l_data;

	// Read the data.
	double* signal_profile_buffer = new double[l_profile+2];

if(__DUMP_SIGNAL_TRACK_MSGS__)
	fprintf(stderr, "Loading %d data values.\n", l_profile);

	// Following is to use the codebase indexing: 1 based indexing.
	fread(&(signal_profile_buffer[1]), sizeof(double), l_profile+1, f_prof);
	
	fclose(f_prof);

	return(signal_profile_buffer);
}

void exclude_regions_from_signal_profiles(double* signal_profile, int l_profile, vector<t_annot_region*>* regions_2_exclude, double* pruned_signal_profile, int& l_pruned_profile)
{
	sort(regions_2_exclude->begin(), regions_2_exclude->end(), sort_regions);

	bool* exclusion_profile = new bool[l_profile+2];
	for(int i = 1; i <= l_profile; i++)
	{
		exclusion_profile[i] = false;
	} // i loop.

	for(int i_reg = 0; i_reg < (int)regions_2_exclude->size(); i_reg++)
	{
		for(int i = regions_2_exclude->at(i_reg)->start; i <= regions_2_exclude->at(i_reg)->end; i++)
		{
			if(i <= l_profile)
			{
				exclusion_profile[i] = true;
			}
		} // i loop.
	} // i_reg loop.

	l_pruned_profile = 1;

	for(int i = 1; i <= l_profile; i++)
	{
		if(exclusion_profile[i])
		{
		}
		else
		{
			pruned_signal_profile[l_pruned_profile] = signal_profile[i];
			l_pruned_profile++;
		}
	} // i loop.

	// Must move the pruned profile index back by one.
	l_pruned_profile--;

	fprintf(stderr, "Pruned the signal of %d values to %d values.\n", l_profile, l_pruned_profile);
	delete [] exclusion_profile;
}

void get_log_plus_one_profile(double* signal_profile, double base, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_log_value = xlog(signal_profile[i]+1.0)/xlog(base);
		signal_profile[i] = cur_log_value;
	} // i loop.
}

void floorize_profile(double* signal_profile, int l_profile)
{
	for(int i = 1; i <= l_profile; i++)
	{
		double cur_floor_val = floor(signal_profile[i]);
		signal_profile[i] = cur_floor_val;
	} // i loop.
}

double* get_zero_profile(int l_profile)
{
	double* cur_profile = new double[l_profile+2];
	for(int i = 0; i <= l_profile; i++)
	{
		cur_profile[i] = 0.0;
	} // i loop.

	return(cur_profile);
}

