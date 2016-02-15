#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "ms_gsl_fft_filter_utils.h"
#include "ms_utils.h"
#include "ms_ansi_string.h"
#include "ms_genomics_coords.h"
#include "ms_annot_region_tools.h"
#include "ms_signal_track_tools.h"
#include "ms_min_max_utils.h"
#include <algorithm>
#include "ms_xlog_math.h"
#include "string.h"

using namespace std;

bool __DUMP_FILTER_MSGS__ = false;

struct t_ext_double_array
{
	double* buffer;	
	int n_vals;
	int buffer_length;
};

t_ext_double_array* alloc_ext_array()
{
	t_ext_double_array* ext_array = new t_ext_double_array();
	ext_array->buffer = new double[1000];
	ext_array->buffer_length = 1000;
	ext_array->n_vals = 0;

	return(ext_array);
}

void add_val(t_ext_double_array* ext_array, double val)
{
	// Check if the buffer memory will be overwhelmed by adding the new value.
	if(ext_array->n_vals+1 >= ext_array->buffer_length)
	{
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Extending the array: %d -> %d\n", ext_array->buffer_length, ext_array->buffer_length * 2);
		int l_ext_ext_buffer = ext_array->buffer_length * 2;
		double* new_ext_buffer = new double [l_ext_ext_buffer];

		// Copy the values.
		for(int i = 0; i < ext_array->n_vals; i++)
		{
			new_ext_buffer[i] = ext_array->buffer[i];
		} // i loop.

		// Copy the new buffer the extended array.		
		delete [] ext_array->buffer;
		ext_array->buffer_length = l_ext_ext_buffer;
		ext_array->buffer = new_ext_buffer;
	}

	// Add the value and return.
	ext_array->buffer[ext_array->n_vals] = val;
	ext_array->n_vals++;
}

void get_next_2_exp(int val, int& larger_exp_val, int& expon)
{
	int cur_val = 1;
	int cur_exp = 0;
	while(cur_val < val)
	{
		cur_exp++;
		cur_val *= 2;
	}

	expon = cur_exp;
	larger_exp_val = cur_val;
}

int sort_doubles_descending(const void* p1, const void* p2)
{
	double val1 = *((double*)p1);
	double val2 = *((double*)p2);

	if(val1 < val2)
	{
		return(-1);
	}
	else if(val1 == val2)
	{
		return(0);
	}
	else
	{
		return(1);
	}
}

double* mean_filter_data(double* track_data,
						int l_track_data, 
						int l_averaging_win,
						int skip_value) // This is the value to be skipped from the window values.
{
	// Do copying from start to end.
    double* cur_filtered_track = new double[l_track_data + 2];

	int half_l_averaging = l_averaging_win / 2;
		
	// Initialize the current pdf.
	double cur_total_win_val = 0;
	bool val_inited = false;

	int prev_avg_start = -1;
	int prev_avg_end = -1;

	// Go over all the positions as the middle of the filtering window.
	for(int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
	{
		cur_filtered_track[cur_win_mid] = 0;

		int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
		int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
		if(!val_inited)
		{
			for(int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				cur_total_win_val += track_data[i];
			} // i loop.

			val_inited = true;
		}
		else
		{
			// Exclude the signal values that got out of the current window.
			for(int i = prev_avg_start; i < cur_avg_start; i++)
			{
				cur_total_win_val -= track_data[i];
			} // i loop.

			// Add the signal values that entered the current window.
			for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
			{
				cur_total_win_val += track_data[i];
			} // i loop.
		}

		// Following does a sanity check for the computed values.
if(__DUMP_FILTER_MSGS__)
{
		double check_val = 0;
		for(int i = cur_avg_start; i <= cur_avg_end; i++)
		{
			check_val += track_data[i];
		} // i loop.

		if(fabs(check_val - cur_total_win_val) > 2)
		{
			fprintf(stderr, "Check failed: %d: %lf, %lf\n", cur_win_mid, check_val, cur_total_win_val);
			exit(0);
		}
}

		cur_filtered_track[cur_win_mid] = cur_total_win_val / (cur_avg_end - cur_avg_start + 1);

		prev_avg_start = cur_avg_start;
		prev_avg_end = cur_avg_end;
	} // main signal filtering loop.

	return(cur_filtered_track);
}

// Median filter: Generates 1 based filtered tracks.
double* median_filter_data(double* track_data,
	int l_track_data, 
	int l_averaging_win,
	int skip_value) // This is the value to be skipped from the window values.
{
	// Do copying from start to end.
    double* cur_filtered_track = new double[l_track_data + 2];

	int n_signal_wins = 0;
	int half_l_averaging = l_averaging_win / 2;

	// Set the maximum value for the histogram.
	int MAX_VAL = -1000*1000;
	int MIN_VAL = 1000*1000;
	for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
	{
		if(MAX_VAL < (int)(track_data[i_sig]))
		{
			MAX_VAL = (int)(track_data[i_sig]);
		}

		if(MIN_VAL > (int)(track_data[i_sig]))
		{
			MIN_VAL = (int)(track_data[i_sig]);
		}
	} // i_sig loop.

	// Make sure that this is the maximum.
	MAX_VAL += 1000;
	MIN_VAL -= 1000;
		
	// Initialize the current pdf.
	int* cur_win_pdf = NULL;
	int cur_win_max = 0;
	int cur_win_min = 1000*1000;
	double* cur_win_vals = new double[l_averaging_win + 2];
	int prev_avg_start = 0;
	int prev_avg_end = 0;

	// Go over all the positions as the middle of the filtering window.
	for(int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
	{
		int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
		int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
		if(cur_win_pdf == NULL)
		{
			cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
			memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
			//t_string::set_byte_buffer(cur_win_pdf, sizeof(int) * (MAX_VAL - MIN_VAL + 1), 0);
			cur_win_pdf -= MIN_VAL;

			// Generate the pdf, get the minimum and maximum in the current window.
			for(int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				if(track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if(cur_win_max < track_data[i])
					{
						cur_win_max = (int)track_data[i];
					}
					
					if(cur_win_min > track_data[i])
					{
						cur_win_min = (int)track_data[i];
					}
				} // skip_Value check.
			} // i loop.
		} // cur_win_pdf is NULL check.
		else
		{
			// Remove the old values from the pdf, add the new values.
			for(int i = prev_avg_start; i < cur_avg_start; i++)
			{
				if(track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]--;
				}
			} // i loop.

			// Update the min and max only for the values that are new in this window.
			for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
			{
				if(track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if(cur_win_max < track_data[i])
					{
						cur_win_max = (int)track_data[i];
					}
					
					if(cur_win_min > track_data[i])
					{
						cur_win_min = (int)track_data[i];
					}
				} // skip_value check.
			}

			// Sanity Check: The total # of points must be equal to the window length.
			int n_total_pts = 0;
			for(int i = cur_win_min; i <= cur_win_max; i++)
			{
				if(i != skip_value)
				{
					n_total_pts += cur_win_pdf[i];
				}
			} // i loop.
		} // cur_win_pdf is NULL check.

		// Count the total number of points without the skip value.
		int n_total = 0;
		for(int i = cur_win_min; i <= cur_win_max; i++)
		{
			if(i != skip_value)
			{
				n_total += cur_win_pdf[i];
			}
		} // i loop.
		
		// Generate the window median.
		int cur_win_median = 0;
		int n_cur_total = 0;
		for(int i = cur_win_min; i <= cur_win_max; i++)
		{
			if(i != skip_value)
			{
				n_cur_total += cur_win_pdf[i];
				if(n_cur_total > n_total/2)
				{
					// We found the median, can break out of the loop.
					cur_win_median = i;
					break;
				}
			}
		} // i loop.

		// Track the minimum and maximum from the histogram to update them for the current window.
		int updated_win_min = 0;			
		for(int i = cur_win_min; i <= cur_win_max; i++)
		{
			if(i != skip_value)
			{
				if(cur_win_pdf[i] > 0)
				{
					updated_win_min = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		int updated_win_max = 0;
		for(int i = cur_win_max; i >= cur_win_min; i--)
		{
			if(i != skip_value)
			{
				if(cur_win_pdf[i] > 0)
				{
					updated_win_max = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		// Set the median.
		//int median_per_pdf = cur_win_median;
		cur_filtered_track[cur_win_mid] = cur_win_median;

		// Update the previous averaging window limits.
		prev_avg_start = cur_avg_start;
		prev_avg_end = cur_avg_end;
		cur_win_min = updated_win_min;
		cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
		// Get the median via qsort and compare as a sanity check.
		int n_valid_vals = 0;
		for(int i = cur_avg_start; i <= cur_avg_end; i++)
		{
			if(track_data[i] != skip_value)
			{
				cur_win_vals[n_valid_vals] = track_data[i];
				n_valid_vals++;
			}
		} // i loop.
		qsort(cur_win_vals, n_valid_vals, sizeof(double), sort_doubles_descending);
		//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

		if(cur_win_vals[n_valid_vals/2] != median_per_pdf)
		{
			fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[n_valid_vals/2], median_per_pdf);
			for(int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				fprintf(stderr, "%lf ", track_data[i]);
			} // i loop.
			fprintf(stderr, "\n");
			getc(stdin);
		}
#endif // _QSORT_CHECK_

		n_signal_wins++;
	} // main signal filtering loop.

	// Free pdf memory.
	delete [] cur_win_vals;
	delete [] (cur_win_pdf+MIN_VAL);

	return(cur_filtered_track);
}

vector<double*>* multiscale_avg_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
{
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);
}

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Initialize the current pdf.
			//int* cur_win_pdf = NULL;
			//int cur_win_max = 0;
			//int cur_win_min = 1000*1000;
			//double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			double cur_win_total = -123123;
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_total == -123123)
				{
					// Get the total by a loop for the current window.
					cur_win_total = 0.0;
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_total += track_data[i];
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_total -= track_data[i];
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_total += track_data[i];
					}
				} // cur_win_pdf is NULL check.

				// Set the median.
				cur_filtered_track[cur_win_mid] = cur_win_total / l_averaging_win;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;

#undef _MANUAL_TOTAL_CHECK_
#ifdef _MANUAL_TOTAL_CHECK_
				// Get the median via qsort and compare as a sanity check.
				double manual_total = 0.0;
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					manual_total += track_data[i];
				} // i loop.

				if(manual_total != cur_win_total)
				{
					fprintf(stderr, "Sanity check failed at manual total check, %s(%d)\n", __FILE__, __LINE__);
					exit(0);
				}
#endif // _MANUAL_TOTAL_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				fprintf(stderr, "Closing file\n");

				// Clean track memory.
				delete [] _cur_filtered_track;
				delete [] cur_filtered_track;
			}
			else if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size(); i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size(); i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
				delete [] cur_filtered_track;
			}
			else
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
} // multiscale_avg_filter_data

vector<double*>* multiscale_median_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	bool return_filtered_tracks,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	const char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %.2f (%.2f)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			//for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			for(int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];
				} // i loop.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int cur_n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					cur_n_total += cur_win_pdf[i];
					if(cur_n_total > n_total/2)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);

if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);

				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				delete [] _cur_filtered_track;

if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Closing file\n");
			}
			
			if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
			}
			
			// Do we need the filtered track back?
			if(return_filtered_tracks)
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}
			else
			{
				// Clean track memory.
				delete [] cur_filtered_track;
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
}

double* mapability_aware_median_filter(double* signal_profile, int l_profile,
	double* scaled_mapability_profile, int l_mapability_profile,
	double max_mapable_signal_2_use_in_filter,
	int l_mapability_filtering_win)
{
	// Generate the signal profile to input to the filtering.
	double* input_signal_profile = new double[l_profile + 2];
	for(int i = 1; i <= l_profile; i++)
	{
		if(l_mapability_profile > i)
		{
			if(scaled_mapability_profile[i] > max_mapable_signal_2_use_in_filter)
			{
				input_signal_profile[i] = -1;
			}
			else
			{
				input_signal_profile[i] = signal_profile[i];
			}
		}
		else
		{
			input_signal_profile[i] = signal_profile[i];
		}
	} // i loop.

if(__DUMP_FILTER_MSGS__)
{
	fprintf(stderr, "Mapability aware filtering.\n");
}

	double* mapable_median_profile = median_filter_data(input_signal_profile,
													l_profile, 
													l_mapability_filtering_win,
													-1);

	double* filtered_signal = new double[l_profile + 5];

	// Compare the actual signal profile with the 
	for(int i = 1; i <= l_profile; i++)
	{
		filtered_signal[i] = MAX(mapable_median_profile[i], signal_profile[i]);
	} // i loop.

	delete [] mapable_median_profile;
	delete [] input_signal_profile;

	return(filtered_signal);
}
	
void get_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
														int l_track_data, 
														double scale_start, double scale_end, double log_scale_step, 
														vector<double>* scales_per_i_scale,
														vector<vector<t_annot_region*>*>* per_scale_minima_regions,
														double min_allowed_filtered_value_scale)
{
	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;

	// Define where which quantile we set as the signal.
	double quantile = 0.5;

	FILE* f_feature_signal_vals = NULL;

if(__DUMP_FILTER_MSGS__)
	f_feature_signal_vals = open_f("features_signal_values.bed", "w");

	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			i_scale++;
			per_scale_minima_regions->push_back(new vector<t_annot_region*>());
			continue;
		}
		else
		{
//if(__DUMP_FILTER_MSGS__)
//{
			//fprintf(stderr, "Processing scale %.2f (%.2f)\n", scale, scale_end);
			fprintf(stderr, "Generating ERs @ scale with %d/%d bins.                      \r", (int)scale, (int)scale_end);
//}

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];
			memset(cur_filtered_track, 0, (l_track_data+1)*sizeof(double));

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 100;
			MIN_VAL -= 100;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window: This must be 1 based.
			for(int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];

					if(n_total > (cur_avg_end - cur_avg_start + 1)*quantile)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

if(__DUMP_FILTER_MSGS__)
{
			fprintf(stderr, "Getting the extrema regions.\n");
}

			// Get the extrema regions for the current filtered regions.
			vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
			vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
			int* derivative_map = new int[l_track_data+2];
			memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
			get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

			// Sort the extremas.
			sort(minima->begin(), minima->end(), sort_extremas_per_posn);
			sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

			// Get the minima regions.
			vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();

if(__DUMP_FILTER_MSGS__)
{
			fprintf(stderr, "Generating filtered minima regions from %d minima.\n", (int)minima->size());
}

			for(int i_m = 0; 
				minima->size() > 0 &&
				i_m < (int)minima->size()-1; 
				i_m++)
			{
				t_annot_region* new_minima = get_empty_region();
				new_minima->chrom = t_string::copy_me_str("XX");
				//new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
				//new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
				new_minima->start = minima->at(i_m)->extrema_posn;
				new_minima->end = minima->at(i_m+1)->extrema_posn;
				new_minima->strand = '+';

				// The maximum in the original signal is found within the effective data that contributes to the filtering process.
				int i_eff_data_start = (1 > new_minima->start-(int)scale/2)?(1):(new_minima->start-(int)scale/2);
				int i_eff_data_end = (l_track_data < new_minima->end+scale/2)?(l_track_data):(new_minima->end+(int)scale/2);
				double original_max = 0.0;
				for(int i = i_eff_data_start; i <= i_eff_data_end; i++)
				{
					if(original_max < track_data[i])
					{
						original_max = track_data[i];
					}
				}

				double filtered_max = 0.0;
				for(int i = new_minima->start; i <= new_minima->end; i++)
				{
					if(filtered_max < cur_filtered_track[i])
					{
						filtered_max = cur_filtered_track[i];
					}
				} // i loop.

if(__DUMP_FILTER_MSGS__)
{
				fprintf(f_feature_signal_vals, "XX\t%d\t%d\t%lf\t%lf\t%lf\n", new_minima->start, new_minima->end, filtered_max, original_max, scale);
}

				// This is where the filtering is done.
				if(filtered_max >= original_max / min_allowed_filtered_value_scale)
				{
					cur_scale_minima_regions->push_back(new_minima);
				}
				else
				{					
					delete new_minima;
				}
			} // i_min loop.

if(__DUMP_FILTER_MSGS__)
{
			fprintf(stderr, "%ld minima regions are generated.\n", cur_scale_minima_regions->size());
}

			per_scale_minima_regions->push_back(cur_scale_minima_regions);

			delete_extrema_nodes(minima);
			delete_extrema_nodes(maxima);

if(__DUMP_FILTER_MSGS__)
{
			// Dump the current filtered track.
			char cur_filtered_data_op_fp[1000];
			sprintf(cur_filtered_data_op_fp, "filtered_%lf.bin", scale);
			dump_per_nucleotide_binary_profile(cur_filtered_track, l_track_data, cur_filtered_data_op_fp);
}

			// Free the current decomposition.
			delete [] derivative_map;
			delete [] cur_filtered_track;

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.
	fprintf(stderr, "Generated SSERs for %d scales.                                \n", (int)scales_per_i_scale->size());

if(__DUMP_FILTER_MSGS__)
		fclose(f_feature_signal_vals);
} // -get_filtered_maxima_regions_multiscale_filtered_data option.


bool get_origin_of_symmetry(double* filtered_data, int l_signal, int search_start_i, int& origin, int l_mirror_check)
{
	if(search_start_i < l_mirror_check)
	{
		search_start_i = l_mirror_check;
	}

	for(int i_c = l_mirror_check; i_c < l_signal-l_mirror_check; i_c++)
	{
		bool current_origin_check = true;
		for(int i = 0; 
			current_origin_check && i < l_mirror_check; 
			i++)
		{
			if(fabs(filtered_data[i_c+i] - filtered_data[i_c-i]) > 0.00001)
			{
				current_origin_check = false;
			}
		} // i loop.

		if(current_origin_check)
		{
			fprintf(stderr, "Found origin of symmetry @ %d\n", i_c);
			origin = i_c;
			return(true);
		}
	} // i_c loop.

	return(false);
}

double* sum_tracks(vector<double*>* multitrack_data, int l_signal)
{
	int n_tracks = (int)multitrack_data->size();
	fprintf(stderr, "Summing %d tracks of %d data points.\n", n_tracks, l_signal);
	double* aggregate_track = new double[l_signal];

	for(int i = 0; i < l_signal; i++)
	{
		double cur_aggregated_val = 0.0;

		// Go over all the tracks and sum their values at this position.
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			cur_aggregated_val += (multitrack_data->at(i_t))[i];
		} // i_t loop.

		aggregate_track[i] = cur_aggregated_val;
	} // i loop.

	return(aggregate_track);
}
