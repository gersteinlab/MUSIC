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
#include <algorithm>
#include "ms_xlog_math.h"
#include <string.h>

using namespace std;

bool __DUMP_FILTER_MSGS__ = false;


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
	fprintf(stderr, "Mapability aware filtering.\n");

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
	for(int i_sig = 0; i_sig < l_track_data; i_sig++)
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
	for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
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
				cur_win_vals[i-cur_avg_start] = track_data[i];
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

bool sort_extremas_per_posn(t_extrema_node* node1, t_extrema_node* node2)
{
	return(node1->extrema_posn < node2->extrema_posn);
}

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale,
	double zero_deriv)
{
	// Go over all the data, identify the plateaus.
	//

	////double* abs_deriv = new double[l_signal];
	//vector<double>* abs_deriv = new vector<double>(l_signal, 0.0);

	//// Following can be implemented with 
	////abs_deriv[0] = 0;
	//abs_deriv->push_back(0);
	//for(int i = 1; i < l_signal; i++)
	//{
	//	abs_deriv->push_back(fabs(data[i] - data[i-1]));
	//} // i loop.

	//sort(abs_deriv->begin(), abs_deriv->end());

	//zero_deriv = abs_deriv->at((int)(l_signal * 0.8));

	//delete abs_deriv;

	//t_1D_hist* deriv_hist = get_1D_dist(abs_deriv, l_signal, 50);

	//int cumul_sum = 0;
	//int i_satisfying_bin = 0;
	//for(int i_bin = 0; i_bin < l_signal; i_bin++)
	//{
	//	cumul_sum += deriv_hist->counts[i_bin];

	//	// Put 20% of the positions into fast changing positions.
	//	if(cumul_sum > (0.8 * l_signal))
	//	{
	//		i_satisfying_bin = i_bin;
	//	}
	//} // i loop.

	//delete [] abs_deriv;
	//zero_deriv = deriv_hist->data_starts[i_satisfying_bin];
	//fprintf(stderr, "zero_deriv = %.5f\n", zero_deriv);	
	//getc(stdin);

	vector<t_extrema_node*>* extrema_nodes = new vector<t_extrema_node*>();

	double* deriv = new double[l_signal];
	vector<int>* plateau_locs = new vector<int>();
	// Following can be implemented with 
	deriv[0] = 0;
	for(int i = 1; i < l_signal; i++)
	{
		deriv[i] = data[i] - data[i-1];

		if(fabs(deriv[i]) < zero_deriv)
		{
			//fprintf(stderr, "Pleateau loc: %d\n", i);
			plateau_locs->push_back(i);
		}
		else
		{
		}
	} // i loop.

	// Form the plateau signal.
	vector<t_plateau*>* plateaus = new vector<t_plateau*>();

	int* plateau_signal = new int[l_signal];
	memset(plateau_signal, 0, sizeof(int) * l_signal);
	int n_pl_pts = (int)plateau_locs->size(); 
	int i = 0; 
	while(i < n_pl_pts)
	{
		int cur_start = plateau_locs->at(i);
		int cur_end = plateau_locs->at(i);

		while(i+1 < n_pl_pts &&
			plateau_locs->at(i+1) == plateau_locs->at(i)+1)
		{			
			cur_end = plateau_locs->at(i+1);

			i++;
		} // plateau finder loop.

		i++;

		t_plateau* new_plateau = new t_plateau();
		new_plateau->start = cur_start;
		new_plateau->end = cur_end;
		//fprintf(stderr, "Pleateau: %d-%d\n", cur_start, cur_end);

		plateaus->push_back(new_plateau);
	} // i loop.	

	// Go over all the plateaus and set the max/min ver plateau.
	int l_win = 1;
	for(int i_pl = 0; i_pl < (int)plateaus->size(); i_pl++)
	{
		if(plateaus->at(i_pl)->start > l_win &&
			plateaus->at(i_pl)->end+l_win < l_signal &&
			deriv[plateaus->at(i_pl)->start-l_win] > 0 &&
			deriv[plateaus->at(i_pl)->end+l_win] < 0)
		{
			t_extrema_node* new_ext_node = new t_extrema_node();
			new_ext_node->extrema_posn = (plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2;
			new_ext_node->extrema_type = EXTREMA_MAX;
			new_ext_node->higher_node = NULL;
			new_ext_node->lower_node = NULL;
			new_ext_node->scale = i_scale;
			new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
			extrema_nodes->push_back(new_ext_node);
			//fprintf(stderr, "Adding max\n");
			//p2n->push_back((plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2);
		}

		if(plateaus->at(i_pl)->start > l_win &&
			plateaus->at(i_pl)->end+l_win < l_signal &&
			deriv[plateaus->at(i_pl)->start-l_win] < 0 &&
			deriv[plateaus->at(i_pl)->end+l_win] > 0)
		{
			t_extrema_node* new_ext_node = new t_extrema_node();
			new_ext_node->extrema_posn = (plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2;
			new_ext_node->extrema_type = EXTREMA_MIN;
			new_ext_node->higher_node = NULL;
			new_ext_node->lower_node = NULL;
			new_ext_node->scale = i_scale;
			new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
			extrema_nodes->push_back(new_ext_node);
			//fprintf(stderr, "Adding max\n");
			//n2p->push_back((plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2);
		}
	} // i_pl loop.

	// Add the sudden changes.
	for(int i = 1; i < l_signal; i++)
	{
		if(i+1 < l_signal)
		{
			if(deriv[i] > zero_deriv &&
				deriv[i+1] < -1 * zero_deriv)
			{
				t_extrema_node* new_ext_node = new t_extrema_node();
				new_ext_node->extrema_posn = i;
				new_ext_node->extrema_type = EXTREMA_MAX;
				new_ext_node->higher_node = NULL;
				new_ext_node->lower_node = NULL;
				new_ext_node->scale = i_scale;
				new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
				extrema_nodes->push_back(new_ext_node);
				//p2n->push_back(i);
			}

			if(deriv[i] < -1 * zero_deriv &&
				deriv[i+1] > zero_deriv)
			{
				t_extrema_node* new_ext_node = new t_extrema_node();
				new_ext_node->extrema_posn = i;
				new_ext_node->extrema_type = EXTREMA_MIN;
				new_ext_node->higher_node = NULL;
				new_ext_node->lower_node = NULL;
				new_ext_node->scale = i_scale;
				new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
				extrema_nodes->push_back(new_ext_node);
				//n2p->push_back(i);
			}
		}
	} // i loop.

	if((int)extrema_nodes->size() == 0 ||
		(int)extrema_nodes->size() == 1)
	{
if(__DUMP_FILTER_MSGS__)
{
		fprintf(stderr, "Found %ld extrema.\n", extrema_nodes->size());
}
		return;
	}

	sort(extrema_nodes->begin(), extrema_nodes->end(), sort_extremas_per_posn);

	t_extrema_node* begin_ext_node = new t_extrema_node();
	begin_ext_node->extrema_posn = 0;
	begin_ext_node->extrema_type = EXTREMA_MIN;
	begin_ext_node->higher_node = NULL;
	begin_ext_node->lower_node = NULL;
	begin_ext_node->scale = i_scale;
	begin_ext_node->height_at_extrema = data[begin_ext_node->extrema_posn];

	t_extrema_node* end_ext_node = new t_extrema_node();
	end_ext_node->extrema_posn = l_signal-1;
	end_ext_node->extrema_type = EXTREMA_MIN;
	end_ext_node->higher_node = NULL;
	end_ext_node->lower_node = NULL;
	end_ext_node->scale = i_scale;
	end_ext_node->height_at_extrema = data[end_ext_node->extrema_posn];

	if(extrema_nodes->at(0)->extrema_type == EXTREMA_MAX)
	{
		begin_ext_node->extrema_type = EXTREMA_MIN;
	}
	else
	{
		begin_ext_node->extrema_type = EXTREMA_MAX;
	}

	if(extrema_nodes->back()->extrema_type == EXTREMA_MAX)
	{
		end_ext_node->extrema_type = EXTREMA_MIN;
	}
	else
	{
		end_ext_node->extrema_type = EXTREMA_MAX;
	}

	extrema_nodes->push_back(begin_ext_node);
	extrema_nodes->push_back(end_ext_node);

	sort(extrema_nodes->begin(), extrema_nodes->end(), sort_extremas_per_posn);

	// Set the derivative map.
	int i_cur_nuc = 0;
	for(int i_ext = 0; i_ext < (int)extrema_nodes->size(); i_ext++)
	{
		if(i_ext > 0 &&
			extrema_nodes->at(i_ext-1)->extrema_type == extrema_nodes->at(i_ext)->extrema_type)
		{
			fprintf(stderr, "Consecutive same extrema: %d, %d: %d\n", 
				extrema_nodes->at(i_ext-1)->extrema_posn, 
				extrema_nodes->at(i_ext)->extrema_posn, 
				extrema_nodes->at(i_ext)->extrema_type);
			exit(0);
		}

		int der_val = 0;
		if(extrema_nodes->at(i_ext)->extrema_type == EXTREMA_MAX)
		{
			der_val = 1;
		}
		else
		{
			der_val = -1;
		}

		for(int i_nuc = i_cur_nuc; i_nuc <= extrema_nodes->at(i_ext)->extrema_posn-1; i_nuc++)
		{
			derivative_sign_map[i_nuc] = der_val;
		} // i_nuc loop.

		derivative_sign_map[extrema_nodes->at(i_ext)->extrema_posn] = 0;
		i_cur_nuc = extrema_nodes->at(i_ext)->extrema_posn+1;
		//cur_type = extrema_nodes->at(i_ext)->extrema_type;
	} // i_ext loop.

if(__DUMP_FILTER_MSGS__)
{
	fprintf(stderr, "Found %ld extrema.\n", extrema_nodes->size());
}
	//for(int i = 0; i < 20; i++)
	//{
	//	fprintf(stderr, "%d: %d\n", extrema_nodes->at(i)->extrema_posn, extrema_nodes->at(i)->extrema_type);
	//}

	for(int i = 0; i < (int)extrema_nodes->size(); i++)
	{
		if(extrema_nodes->at(i)->extrema_type == EXTREMA_MAX)
		{
			maxes->push_back(extrema_nodes->at(i));
		}
		else
		{
			mins->push_back(extrema_nodes->at(i));
		}
	} // i loop.

	// Free memory. These are pretty large chunks of memory.
	delete [] deriv;
	delete [] plateau_signal;
	for(int i_p = 0; i_p < (int)plateaus->size(); i_p++)
	{
		delete(plateaus->at(i_p));
	} // i_p loop.
	delete plateau_locs;
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

	FILE* f_feature_signal_vals = NULL;
if(__DUMP_FILTER_MSGS__)
	f_feature_signal_vals = open_f("features_signal_values.bed", "w");

	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);


		//if(scale < scale_start)
		//{
		//	decomps->push_back(NULL);

		//	if(compute_extrema_regions)
		//	{
		//		if(per_scale_minima_regions != NULL)
		//		{
		//			
		//		}

		//		if(per_scale_maxima_regions != NULL)
		//		{
		//			per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
		//		}
		//	}

		//	i_scale++;
		//	continue;
		//}


		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			i_scale++;
			per_scale_minima_regions->push_back(new vector<t_annot_region*>());
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

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
			for(int i_sig = 0; i_sig < l_track_data; i_sig++)
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
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
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

					//if(n_total > l_averaging_win/2)
					if(n_total > ((cur_avg_end - cur_avg_start + 1))/2)
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
			fprintf(stderr, "Getting the extrema regions.\n");

			// Get the extrema regions for the current filtered regions.
			vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
			vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
			int* derivative_map = new int[l_track_data+2];
			memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
			double zero_deriv = pow(10.0, -10.0);
			get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0, zero_deriv);

			// Sort the extremas.
			sort(minima->begin(), minima->end(), sort_extremas_per_posn);
			sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

			// Get the minima regions.
			vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();

if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Generating filtered minima regions from %ld minima.\n", minima->size());

			for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
			{
				t_annot_region* new_minima = get_empty_region();
				new_minima->chrom = t_string::copy_me_str("XX");
				new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
				new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
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
				fprintf(f_feature_signal_vals, "XX\t%d\t%d\t%lf\t%lf\t%lf\n", new_minima->start, new_minima->end, filtered_max, original_max, scale);

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
			fprintf(stderr, "%ld minima regions are generated.\n", cur_scale_minima_regions->size());

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

if(__DUMP_FILTER_MSGS__)
		fclose(f_feature_signal_vals);
} // -get_filtered_maxima_regions_multiscale_filtered_data option.

void delete_extrema_nodes(vector<t_extrema_node*>* extrema_nodes)
{
	for(int i_n = 0; i_n < (int)extrema_nodes->size(); i_n++)
	{
		delete extrema_nodes->at(i_n);
	}

	delete extrema_nodes;
}


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

