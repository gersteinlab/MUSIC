#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "ms_signal_track_tools.h"
#include "ms_ansi_string.h"
#include "ms_utils.h"
#include "ms_annot_region_tools.h"
#include "ms_genomics_coords.h"
#include "ms_nomenclature.h"
#include "ms_xlog_math.h"
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include "ms_gsl_fft_filter_utils.h"
#include "ms_mapped_read_tools.h"
#include "ms_profile_normalization.h"
#include "ms_utils.h"
#include "ms_signal_enrichment_utils.h"
#include "ms_signal_enrichment_utils.h"

bool __DUMP_SIGNAL_ENRICHMENT_MSGS__ = false;

double* buffer_log_factorials(int n)
{
	double* factorials = new double[n+2];
	factorials[0] = xlog(1.0);
	for(int i = 1; i <= n; i++)
	{
		factorials[i] = xlog_mul(factorials[i-1], xlog(i));
	} // i loop.

	return(factorials);
}

/*
Given the counts, compute the p-value for the enrichment of the enriched value with respect to background.
*/
double get_binomial_p_val_per_enriched_background_values(int enriched_value, int background_value, double* _log_factorials)
{
	int grand_total_int = (int)(enriched_value + background_value);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if(_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int+3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(.5);
	double cur_region_log_p_val = xlog(0.0);
	for(int i = 0; i <= background_value; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Free memory if it is allocated in the function.
	if(_log_factorials == NULL)
	{
		delete [] log_factorials;
	}
	else
	{
	}

	return(cur_region_log_p_val);
}

double get_poisson_thresholds_per_win(double* signal_profile,
								int l_profile,
								int l_win,
								int win_start,
								double target_fdr,
								int min_thresh,
								int max_thresh)
{
	// Count the number of reads per window.
	double total_sig_in_cur_window = 0.0;
	for(int i = win_start; i <= win_start + l_win; i++)
	{
		total_sig_in_cur_window += signal_profile[i];
	} // i loop.

	// Get the total amount of signal in the current window.
	if(total_sig_in_cur_window == 0.0)
	{
		return(0.0);
	}

	if(l_win == 0)
	{
		return(0.0);
	}
	else
	{
		// This is the average of poisson distribution for the heights in this window.
		double avg_read_depth = total_sig_in_cur_window / l_win;
		double lambda = avg_read_depth;
		double log_lambda = xlog(lambda);

		double log_target_cdf = xlog(1.0 - target_fdr);

		// Update the log_cdf value for thresholds smaller than min_threshold.
		double current_log_cdf = -1.0 * lambda;
		double current_factor = -1.0 * lambda; // This is the probability value differentially updated in CDF computation.
		for(int thresh = 1; 
			thresh < min_thresh;
			thresh++)
		{
			double new_multiplier = xlog_div(log_lambda, xlog((double)thresh));
			current_factor = xlog_mul(current_factor, new_multiplier);
			current_log_cdf = xlog_sum(current_log_cdf, current_factor);
		}

		// For thresholds between min and max thresholds, compare the value of (log)CDF with the (log)target FDR.
		bool found_thresh = false;
		double cur_threshold = 0;
		for(int thresh = min_thresh; 
			thresh <= max_thresh && !found_thresh; 
			thresh++)
		{
			double new_multiplier = xlog_div(log_lambda, xlog((double)thresh));
			current_factor = xlog_mul(current_factor, new_multiplier);
			current_log_cdf = xlog_sum(current_log_cdf, current_factor);

			if(current_log_cdf > log_target_cdf)
			{
				found_thresh = true;
				cur_threshold = thresh; // Set the threshold in this window.
			}
		} // thresh loop.

		if(cur_threshold > 0)
		{
			return(cur_threshold);
		}
		else
		{
			return(0.0);
		}
	} // Check for positive number of reads in the window.
} // get_poisson_thresholds_per_win

double get_geometric_thresholds_per_win(double* signal_profile,
								int l_profile,
								int l_win,
								int win_start,
								double target_fdr,
								int min_thresh,
								int max_thresh)
{
	// Count the number of reads per window.
	double total_sig_in_cur_window = 0.0;
	for(int i = win_start; i <= win_start + l_win; i++)
	{
		total_sig_in_cur_window += signal_profile[i];
	} // i loop.

	// Get the total amount of signal in the current window.
	if(total_sig_in_cur_window == 0.0)
	{
		return(0.0);
	}

	if(l_win == 0)
	{
		return(0.0);
	}
	else
	{
		// This is the average of poisson distribution for the heights in this window.
		double avg_read_depth = total_sig_in_cur_window / l_win;

		if(floor(avg_read_depth) == 0)
		{
			return(1);
		}

		double p = 1/avg_read_depth;
		double log_one_minus_p = xlog(1-p);

		double log_target_cdf = xlog(1.0 - target_fdr);

		// For thresholds between min and max thresholds, compare the value of (log)CDF with the (log)target FDR.
		bool found_thresh = false;
		double cur_threshold = 0;
		for(int thresh = min_thresh; 
			thresh <= max_thresh && !found_thresh; 
			thresh++)
		{
			double current_log_cdf = xlog_sub(0.0, thresh * log_one_minus_p);

			if(current_log_cdf > log_target_cdf)
			{
				found_thresh = true;
				cur_threshold = thresh; // Set the threshold in this window.
			}
		} // thresh loop.

		if(cur_threshold > 0)
		{
			return(cur_threshold);
		}
		else
		{
			return(0.0);
		}
	} // Check for positive number of reads in the window.
} // get_geometric_thresholds_per_win

// Following is used to make the no-control based p-value estimates stable.
double* get_per_win_profile_per_nuc_profile(double* profile, int l_profile, int l_win)
{
	fprintf(stderr, "Generating per window profile for %d long windows.\n", l_win);

	double* per_win_profile = new double[l_profile/l_win + 50];
	int cur_win_start = 1;
	//int thresh_win_length = 1000*1000;
	//int cur_thresh_win_start = 1 - thresh_win_length;
	//double cur_thresh_win_thresh = 0;
	int i_win = 1;
	int n_windows_reset = 0;
	while(1)
	{
		int l_cur_win = l_win;
		if(cur_win_start+l_cur_win > l_profile)
		{
			l_cur_win = l_profile - cur_win_start - 1;
		}

		if(cur_win_start > l_profile)
		{
			break;
		}

		// Set the current value.
		per_win_profile[i_win] = 0;
		//int n_vals_above_thresh = 0;
		for(int i = cur_win_start; i <= cur_win_start+l_cur_win; i++)
		{
			per_win_profile[i_win] += profile[i];
			//if(profile[i_win] > cur_thresh_win_thresh)
			//{
			//	n_vals_above_thresh++;
			//}
		}
		
		cur_win_start += l_win;
		i_win++;
	} // cur_win_start loop.

	fprintf(stderr, "Resetted %d window's signal values.\n", n_windows_reset);
	return(per_win_profile);
}

double get_fdr_based_bckgrnd_window_value(double* signal_profile, int l_profile,
												double target_fdr,
												int thresh_win_length,
												int l_background_est_win,
												int l_pval_norm_win,
												double significant_log_p_val_cutoff,
												double scaling_factor)
{
	fprintf(stderr, "Generating per window profile for %d long windows.\n", l_background_est_win);

	// This is all the window values: We choose the 
	vector<double>* all_win_vals = new vector<double>();

	// This is the window values for which we are use as true negatives 
	vector<double>* true_negative_win_vals = new vector<double>();

	double* per_win_profile = new double[l_profile/l_background_est_win + 50];
	int cur_win_start = 1;
	int cur_thresh_win_start = 1 - thresh_win_length;
	double cur_thresh_win_thresh = 0;
	//int n_windows_reset = 0;
	while(1)
	{
		int l_cur_win = l_background_est_win;
		if(cur_win_start+l_cur_win > l_profile)
		{
			l_cur_win = l_profile - cur_win_start - 1;
		}

		// Get the current threshold for the window: Did the window start get out of the threshold window range?
		if((cur_thresh_win_start+thresh_win_length) < cur_win_start)
		{
			// Push the threshold forward.
			cur_thresh_win_start += thresh_win_length;

			if(cur_thresh_win_start + thresh_win_length > l_profile)
			{
				thresh_win_length = l_profile - cur_thresh_win_start;
			}

#ifdef __POISSON_BCKGRND__
			cur_thresh_win_thresh = get_poisson_thresholds_per_win(signal_profile,
																	l_profile,
																	thresh_win_length,
																	cur_thresh_win_start,
																	.05,
																	1,
																	5000);
#elif defined __GEOMETRIC_BCKGRND__
			cur_thresh_win_thresh = get_geometric_thresholds_per_win(signal_profile,
																	l_profile,
																	thresh_win_length,
																	cur_thresh_win_start,
																	.05,
																	1,
																	5000);

				//double cur_win_thresh = get_geometric_thresholds_per_win(signal_profile,
				//								l_profile,
				//								cur_l_win,
				//								cur_win_start,
				//								fdr,
				//								1,
				//								max_thresh);
#else
			#error "Did not select a background model."
#endif // Bckgrnd selection
		}

		if(cur_win_start > l_profile)
		{
			break;
		}

		// Set the current value.
		double cur_win_val = 0;
		int n_vals_above_thresh = 0;
		for(int i = cur_win_start; i <= cur_win_start+l_cur_win; i++)
		{
			cur_win_val += signal_profile[i];
			if(signal_profile[i] > cur_thresh_win_thresh)
			{
				n_vals_above_thresh++;
			}
		}

		// If therei is signal in this window, add this to the all window values.
		if(cur_win_val > 0)
		{
			all_win_vals->push_back(cur_win_val);

			// If there are less than 5% of the values in the current window above threshold, then this is assumed an unenriched window.
			if(n_vals_above_thresh <= (.50 * l_cur_win))
			{
				// Scale the current window value wrt the p-value normalization window length.
				cur_win_val /= ((double)l_background_est_win / (double)l_pval_norm_win);
				true_negative_win_vals->push_back(cur_win_val);
			}
		}

		cur_win_start += l_background_est_win;
	} // cur_win_start loop.

	fprintf(stderr, "Generated %ld windows and %ld true negative windows.\n", all_win_vals->size(), true_negative_win_vals->size());

	// Sort the window values.
	sort(all_win_vals->begin(), all_win_vals->end());

	fprintf(stderr, "Median value: %lf\nQuantile value: %lf\n", all_win_vals->at((int)all_win_vals->size()/2), all_win_vals->at((int)all_win_vals->size()/10));

	// Start from the smallest window value.
	double target_fdr_bckgrnd_win_val = 0.0;
	double prev_bckgrnd_win_val = 0.0;

	double cur_bckgrnd_win_val = all_win_vals->at(0);
	while(1)
	{
		int cur_bckgrnd_win_val_int = (int)(floor(cur_bckgrnd_win_val / scaling_factor));

		int n_significant_tn_wins = 0;
		
		// go over all the true negative windows and count the number of significant window.
		for(int i_tn_win = 0; i_tn_win < (int)true_negative_win_vals->size(); i_tn_win++)
		{
			int cur_tn_val_int = (int)(floor(true_negative_win_vals->at(i_tn_win) / scaling_factor));
			double cur_pval = get_binomial_p_val_per_enriched_background_values(cur_tn_val_int, cur_bckgrnd_win_val_int, NULL);

			// Is this window significant?
			if(cur_pval < significant_log_p_val_cutoff)
			{
				n_significant_tn_wins++;
			}
		} // i_tn_win loop.

		prev_bckgrnd_win_val = cur_bckgrnd_win_val;

		fprintf(stderr, "%lf -> %d TN wins significant (%lf)\n", cur_bckgrnd_win_val, n_significant_tn_wins, (double)n_significant_tn_wins / (int)true_negative_win_vals->size());

		// Are we below requested target fdr?
		if(n_significant_tn_wins <= target_fdr * (int)true_negative_win_vals->size())
		{
			target_fdr_bckgrnd_win_val = cur_bckgrnd_win_val;
			fprintf(stderr, "Reached target FDR @ window value %lf\n", target_fdr_bckgrnd_win_val);
			break;
		}

		cur_bckgrnd_win_val += (scaling_factor);
	} // i_w loop.

	delete all_win_vals;
	delete true_negative_win_vals;
	delete [] per_win_profile;

	return(target_fdr_bckgrnd_win_val);
}

double* get_per_win_median_profile_per_background_window_profile(double* per_bckgrnd_win_signal,
	int n_bckgrnd_wins,
	int l_background_est_win,
	int l_background_norm_win)
{
	fprintf(stderr, "Generating per window median profile for %d long windows.\n", l_background_norm_win);

	double* median_profile = new double[n_bckgrnd_wins+10];

	//int cur_win_i = 0;
	vector<double>* cur_sample_vals = new vector<double>();
	for(int i_win = 0; i_win < n_bckgrnd_wins; i_win++)
	{
		if(per_bckgrnd_win_signal[i_win] > 0)
		{
			cur_sample_vals->push_back(per_bckgrnd_win_signal[i_win]);
		}
	} // cur_win_i loop.

	if(cur_sample_vals->size() == 0)
	{
		fprintf(stderr, "Seems like there is no signal, exiting.\n");
		exit(0);
	}

	sort(cur_sample_vals->begin(), cur_sample_vals->end());

#undef __DUMP_PER_WIN_VALS__
#ifdef __DUMP_PER_WIN_VALS__
	FILE* f_per_win_vals = open_f("per_win_vals.txt", "w");
	for(int i_win = 0; i_win < (int)cur_sample_vals->size(); i_win++)
	{
		fprintf(f_per_win_vals, "%lf\n", cur_sample_vals->at(i_win));
	} // i_win loop.
	fclose(f_per_win_vals);
#endif // __DUMP_PER_WIN_VALS__
	// Use global value with the local median estimates to generate the background per window total signal.
	double global_val = cur_sample_vals->at((int)cur_sample_vals->size()/10);
	cur_sample_vals->clear();
	fprintf(stderr, "Global median signal per %d long window is %lf\n", l_background_norm_win, global_val);

	// Set all the values to global value.
	for(int i_win = 0; i_win <= n_bckgrnd_wins+2; i_win++)
	{
		median_profile[i_win] = global_val;
	} // cur_win_i loop.

#define __GLOBAL_SAMPLE_VALS__
#undef __LOCAL_WIN_SAMPLE_VALS__
#ifdef __GLOBAL_SAMPLE_VALS__
	delete cur_sample_vals;
	return(median_profile);
#elif defined __LOCAL_WIN_SAMPLE_VALS__
	int n_norm_wins_per_est_wins = l_background_est_win / l_background_norm_win;

	// Set the median value.
	for(int cur_win_i = 0; cur_win_i <= n_bckgrnd_wins; cur_win_i++)
	{
		cur_sample_vals->clear();
		int cur_reg_est_win_start_i = (cur_win_i > n_norm_wins_per_est_wins/2)?(cur_win_i - n_norm_wins_per_est_wins/2):(1);
		int cur_reg_est_win_end_i = ((cur_reg_est_win_start_i + n_norm_wins_per_est_wins) < n_bckgrnd_wins)?(cur_reg_est_win_start_i + n_norm_wins_per_est_wins):(n_bckgrnd_wins);

		for(int i_win = cur_reg_est_win_start_i; i_win < cur_reg_est_win_end_i; i_win++)
		{
			// Do not use the values that are exactly one.
			if(per_bckgrnd_win_signal[i_win] > 0)
			{
				cur_sample_vals->push_back(per_bckgrnd_win_signal[i_win]);
			}
		} // i_win loop.

		sort(cur_sample_vals->begin(), cur_sample_vals->end());

		// Set the median value; if it is smaller than global value, set it to global value.
		if(cur_sample_vals->size() == 0)
		{
			median_profile[cur_win_i] = global_val;
		}
		else
		{
			median_profile[cur_win_i] = cur_sample_vals->at((int)cur_sample_vals->size()/2);
		}
		if(global_val > median_profile[cur_win_i])
		{
			median_profile[cur_win_i] = global_val;
		}
	} // cur_win_i loop.
	delete(cur_sample_vals);

	return(median_profile);
#endif // __LOCAL_WIN_SAMPLE_VALS__
}

double get_binomial_pvalue_per_region_no_control(double* signal_profile, int l_profile, 
	int region_start, int region_end, 
	double per_win_bckgrnd_win_val,
	double* _log_factorials, 
	double scaling_factor,
	int l_norm_win,
	bool normalize) 
{
	// Count the signal values.
	double cur_reg_profile_sig = 0.0;	

	// Get the total profile signal in the current region.
	for(int i = region_start; i <= region_end; i++)
	{
		if(i < l_profile)
		{
			cur_reg_profile_sig += signal_profile[i];
		}
	} // i loop.

	// Normalize with the scaling factor.
	//int cur_reg_profile_sig_int = (int)(floor(cur_reg_profile_sig / scaling_factor));
	cur_reg_profile_sig = floor(cur_reg_profile_sig / scaling_factor);

	// Normalize the profile signal: If normalization is requested, use the background normalization window to do normalization.
	int cur_reg_profile_sig_int = (int)(cur_reg_profile_sig);
	if(normalize)
	{
		if((region_end-region_start+1) > l_norm_win)
		{
			double cur_profile_sig_val = floor(cur_reg_profile_sig / ((double)(region_end-region_start+1)/l_norm_win));
			cur_reg_profile_sig_int = (int)(cur_profile_sig_val);

			if(cur_reg_profile_sig_int < 0)
			{
				fprintf(stderr, "Overflow error in p-value computation, %lf, %d\n", cur_profile_sig_val, cur_reg_profile_sig_int);
				exit(0);
			}
		}
	}
		
	// Choose the minimum of the median estimates.
	double cur_reg_control_sig = per_win_bckgrnd_win_val;
	if(cur_reg_control_sig == 0)
	{
		fprintf(stderr, "The control signal is set to 0, which should not have happened.\n");
		exit(0);
	}

	// Scale with the scaling factor, e.g., read length.
	int cur_reg_control_sig_int = (int)(floor(cur_reg_control_sig / scaling_factor));

	// If the normalization is not requested, scale the control signal to the length of the region.
	if( !normalize || (normalize && ((region_end - region_start + 1) < l_norm_win)) )
	{
		cur_reg_control_sig_int *= (region_end - region_start + 1) / l_norm_win;
	}

	double cur_p_val = get_binomial_p_val_per_enriched_background_values(cur_reg_profile_sig_int, cur_reg_control_sig_int, NULL);

	return(cur_p_val);
}

double get_binomial_pvalue_per_region_no_control2(double* signal_profile, int l_profile, 
	int region_start, int region_end, 
	double* per_bckgrnd_per_win_median_signal,
	double* _log_factorials, 
	double scaling_factor,
	int l_background_est_win, 
	int l_background_norm_win,
	bool normalize) 
{
	// Count the signal values.
	double cur_reg_profile_sig = 0.0;	

	// Get the total profile signal in the current region.
	for(int i = region_start; i <= region_end; i++)
	{
		if(i < l_profile)
		{
			cur_reg_profile_sig += signal_profile[i];
		}
	} // i loop.

	// Normalize with the scaling factor.
	int cur_reg_profile_sig_int = (int)(floor(cur_reg_profile_sig / scaling_factor));

	// Normalize the profile signal: If normalization is requested, use the background normalization window to do normalization.
	if(normalize)
	{
		if((region_end-region_start+1) > l_background_norm_win)
		{
			cur_reg_profile_sig_int = (int)(floor(cur_reg_profile_sig_int / ((double)(region_end-region_start+1)/l_background_norm_win)));
		}
	}
		
	// Following chooses the background median estimate from the windows right next to the region of interest, this is necessary for long regions.
	//int n_bckgrnd_norm_wins = l_profile / l_background_norm_win;
	//int left_median_est_win_mid_i = (region_start - l_background_est_win/2) / l_background_norm_win;
	//int right_median_est_win_mid_i = (region_end + l_background_est_win/2) / l_background_norm_win;
	int median_profile_win_i = ((region_end + region_start)/2) / l_background_norm_win;

	// Get 3 values, self estimate, left window estimate, right window estimate.
	double self_win_est = per_bckgrnd_per_win_median_signal[median_profile_win_i];
	//double left_win_est = (left_median_est_win_mid_i>0)?(per_bckgrnd_per_win_median_signal[left_median_est_win_mid_i]):(1000*1000);
	//double right_win_est = (right_median_est_win_mid_i<n_bckgrnd_norm_wins)?(per_bckgrnd_per_win_median_signal[right_median_est_win_mid_i]):(1000*1000);

	// Choose the minimum of the median estimates.
	//double cur_reg_control_sig = MIN(MIN(self_win_est, left_win_est), right_win_est);
	double cur_reg_control_sig = self_win_est;
	if(cur_reg_control_sig == 0)
	{
		fprintf(stderr, "The control signal is set to 0, which should not have happened.\n");
		exit(0);
	}

	// Scale with the scaling factor, e.g., read length.
	int cur_reg_control_sig_int = (int)(floor(cur_reg_control_sig / scaling_factor));

	// If the normalization is not requested, scale the control signal to the length of the region.
	if( !normalize || (normalize && ((region_end - region_start + 1) < l_background_norm_win)) )
	{
		cur_reg_control_sig_int *= (region_end - region_start + 1) / l_background_norm_win;
	}

	double cur_p_val = get_binomial_p_val_per_enriched_background_values(cur_reg_profile_sig_int, cur_reg_control_sig_int, NULL);

	return(cur_p_val);
}

double get_per_win_percentile_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, 
	int region_start, int region_end, 
	double* _log_factorials, double scaling_factor,
	int l_win,
	double win_rank_per_region)
{
	// Divide the region into windows.
	//int n_wins = (region_end - region_start + 1) / l_norm_win;
	//double* per_win_p_val = new double[n_wins + 2];
	vector<double>* per_win_p_vals = new vector<double>();

	//int i_win = 0;
	int cur_start = region_start;
	int cur_end = 0;
	while(cur_start < region_end)
	{
		if(cur_start+l_win < region_end)
		{
			cur_end = cur_start+l_win;
		}
		else
		{
			// Following makes sure that the regions with size smaller than window size is processed as they are.
			cur_end = region_end;
		}

		double cur_win_p_val = get_binomial_pvalue_per_region(signal_track_data, control_track_data, cur_start, cur_end, NULL, 200,
			false, 0);

		per_win_p_vals->push_back(cur_win_p_val);

		// Move to the next window.
		cur_start += l_win;
	} // cur_start loop.

	// Sort the p-values.
	sort(per_win_p_vals->begin(), per_win_p_vals->end());

	int i_perc = (int)((double)((int)per_win_p_vals->size()) * win_rank_per_region);

	double p_val = per_win_p_vals->at(i_perc);
	delete per_win_p_vals;
	return(p_val);
}

// Following is to used for computing the p-values in punctate/TF peak calls where the p-value is computed from the largest window signal the region, which is not suitable for broad marks.
double get_binomial_pvalue_per_region_neighborhood_window_control_per_max_profile_signal(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win, int l_signal, int l_control, int n_vicinity_wins)
{
	// Get the total profile signal: 
	double total_profile_signal = 0;
	if(region_end - region_start + 1 > l_norm_win)
	{
		double cur_win_total = 0.0;
		for(int i = region_start; i < region_start + l_norm_win; i++)
		{
			cur_win_total += floor(signal_track_data[i]);
		} // i loop.

		total_profile_signal = cur_win_total;

		// Following loop computes the per window total signals.
		int cur_win_start = region_start;
		while(1)
		{
			// Remove the value that just left the window.
			cur_win_total -= floor(signal_track_data[cur_win_start]);

			// Add the newly entering value to the window.
			cur_win_total += floor(signal_track_data[cur_win_start+l_norm_win]);

			if(cur_win_total > total_profile_signal)
			{
				total_profile_signal = cur_win_total;
			}
			cur_win_start++;

			// Are we at the end of the region?
			if(cur_win_start+l_norm_win >= region_end)
			{
				break;
			}
		} // window loop.
	}
	else
	{
		for(int i = region_start; i <= region_end; i++)
		{
			total_profile_signal += floor(signal_track_data[i]);
		} // i loop.

		// Normalize the profile signal.
		//total_profile_signal = total_profile_signal / ((double)(region_end - region_start+1) / l_norm_win);
	}

	// Get the control signal.
	int l_reg = (region_end - region_start + 1);
	int ext_region_start = region_start - n_vicinity_wins * l_reg;
	int ext_region_end = region_end + n_vicinity_wins * l_reg;
	if(ext_region_start < 1)
	{
		ext_region_start = 1;
	}

	if(ext_region_end > l_control)
	{
		ext_region_end = l_control;
	}

	double total_control_signal = 0;

	double cur_win_total = 0.0;
	for(int i = ext_region_start; i < ext_region_start + l_norm_win; i++)
	{
		cur_win_total += floor(control_track_data[i]);
	} // i loop.

	total_control_signal = cur_win_total;

	// Following loop computes the per window total control signals.
	int cur_win_start = ext_region_start;
	while(1)
	{
		// Remove the value that just left the window.
		cur_win_total -= floor(control_track_data[cur_win_start]);

		// Add the newly entering value to the window.
		cur_win_total += floor(control_track_data[cur_win_start+l_norm_win]);

		if(cur_win_total > total_profile_signal)
		{
			total_profile_signal = cur_win_total;
		}
		cur_win_start++;

		// Are we at the end of the region?
		if(cur_win_start+l_norm_win >= ext_region_end)
		{
			break;
		}
	} // window loop.

	// Scale the control and profile signals.
	total_control_signal /= scaling_factor;
	total_profile_signal /= scaling_factor;

	if(floor(total_control_signal) == 0)
	{
		total_control_signal = 1.0;
	}

	// The total signal hit 0, return 1.
	if(floor(total_profile_signal) == 0)
	{
		//fprintf(stderr, "Hit 0.\n");
		// Return the most insignificant p-value.
		return(xlog(1.0));
	}

	int grand_total_int = (int)(total_profile_signal + total_control_signal);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if(_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int+3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(.5);
	double cur_region_log_p_val = xlog(0.0);
	for(int i = 0; i <= total_control_signal; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Free memory if it is allocated in the function.
	if(_log_factorials == NULL)
	{
		delete [] log_factorials;
	}

	return(cur_region_log_p_val);
} // get_binomial_pvalue_per_region_neighborhood_window_control_per_max_profile_signal

double get_binomial_pvalue_per_region_neighborhood_window_control(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win, int l_signal, int l_control, int n_vicinity_wins)
{
	// Get the total profile signal.
	double total_profile_signal = 0;
	for(int i = region_start; i <= region_end; i++)
	{
		total_profile_signal += floor(signal_track_data[i]);
	} // i loop.

	// If the region size is smaller than 2000 base pairs, get the control signal from a larger region.
	//double total_control_signal = 0;
	int eff_control_region_start = region_start;
	int eff_control_region_end = region_end;
	if((eff_control_region_end - eff_control_region_start+1) <= l_norm_win)
	{
		eff_control_region_start = (region_start+region_end)/2 - l_norm_win;
		eff_control_region_end = (region_start+region_end)/2 + l_norm_win;
	}

	if(eff_control_region_start < 0)
	{
		eff_control_region_start = 1;
	}

	if(eff_control_region_end > l_control)
	{
		eff_control_region_end = l_control;
	}

	if(eff_control_region_end ==eff_control_region_start)
	{
		return(0.0);
		//eff_control_region_end++;
	}

	// Look at the vicinity windows.
	//int n_vicinity_wins = 1;

	// IF the length of the region is smaller than normalization window, 
	//if((region_end-region_start) < l_norm_win)
	//{
	//	n_vicinity_wins = 1;
	//}

	double total_control_signal = 0;
	for(int cur_start = MAX(1, eff_control_region_start-n_vicinity_wins*(eff_control_region_end-eff_control_region_start)); 
		cur_start < MIN(l_control, eff_control_region_end + n_vicinity_wins*(eff_control_region_end-eff_control_region_start));
		cur_start += (eff_control_region_end-eff_control_region_start))
	{
		double cur_total_control_signal = 0;
		for(int i = cur_start; i <= MIN(l_control, cur_start+(eff_control_region_end-eff_control_region_start)); i++)
		{
			cur_total_control_signal += floor(control_track_data[i]);
		} // i loop.

		total_control_signal = MAX(total_control_signal, cur_total_control_signal);
	} // i loop.

	//// Get the control signal: Choose the maximum of the 3 windows.
	//double total_control_signal1 = 0;
	//for(int i = eff_control_region_start; i <= eff_control_region_end; i++)
	//{
	//	total_control_signal1 += floor(control_track_data[i]);
	//} // i loop.

	//// 2nd control value is the window to the right.
	//double total_control_signal2 = 0.0;
	//for(int i = eff_control_region_end; i <= MIN(l_control, 2*eff_control_region_end-eff_control_region_start); i++)
	//{
	//	total_control_signal2 += floor(control_track_data[i]);
	//} // i loop.

	//// 3rd control value is the window to the left.
	//double total_control_signal3 = 0.0;
	//for(int i = MAX(1, 2*eff_control_region_start-eff_control_region_end); i <= eff_control_region_start; i++)
	//{
	//	total_control_signal3 += floor(control_track_data[i]);
	//} // i loop.
	//total_control_signal = MAX(total_control_signal1, MAX(total_control_signal2, total_control_signal3));

	// Finally, normalize the control signal from effective window size to the original window size.
	total_control_signal = total_control_signal * (double)(region_end - region_start+1) / (double)(eff_control_region_end - eff_control_region_start+1);

	// Scale the control and profile signals.
	total_control_signal /= scaling_factor;
	total_profile_signal /= scaling_factor;

	// Normalize both signals with respect to normalization region, if they are larger than the region size.
	if(length_normalization && 
		(region_end - region_start+1) > l_norm_win)
	{
		total_profile_signal = total_profile_signal / ((double)(region_end - region_start+1) / l_norm_win);
		total_control_signal = total_control_signal / ((double)(region_end - region_start+1) / l_norm_win);

		if(floor(total_control_signal) == 0)
		{
			total_control_signal = 1.0;
		}

		// The total signal hit 0, return 1.
		if(floor(total_profile_signal) == 0)
		{
			//fprintf(stderr, "Hit 0.\n");
			// Return the most insignificant p-value.
			return(xlog(1.0));
		}
	}
	
	int grand_total_int = (int)(total_profile_signal + total_control_signal);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if(_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int+3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(.5);
	double cur_region_log_p_val = xlog(0.0);
	for(int i = 0; i <= total_control_signal; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Free memory if it is allocated in the function.
	if(_log_factorials == NULL)
	{
		delete [] log_factorials;
	}

	return(cur_region_log_p_val);
}

double get_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win)
{
	// Get the p-value for the current region.
	double total_profile_signal = 0;
	double total_control_signal = 0;
	for(int i = region_start; i <= region_end; i++)
	{
		total_profile_signal += floor(signal_track_data[i]);
		total_control_signal += floor(control_track_data[i]);	
	} // i loop.

	//// If the region size is smaller than 2000 base pairs, get the control signal from a larger region.
	//if((region_end - region_start+1) <= 2000)
	//{
	//	int eff_reg_start = (region_start+region_end)/2 - 2500;
	//	int eff_reg_end = (region_start+region_end)/2 + 2500;

	//	for(int i = eff_reg_start; i <= eff_reg_end; i++)
	//	{
	//		total_control_signal += floor(control_track_data[i]);
	//	} // i loop.

	//	// Re-normalize the control signal.
	//	total_control_signal = total_control_signal * (region_end - region_start+1) / 5000;
	//}

	total_profile_signal /= scaling_factor;

	//// Get the control signal: Make sure there is enough reads.
	//double total_control_signal = 0.0;
	//int ext_region_start = region_start;
	//int ext_region_end = region_end;
	//while(total_control_signal < MIN_N_CONTROL_SIGNAL/scaling_factor &&
	//	ext_region_start > 0)
	//{
	//	for(int i = ext_region_start; i <= ext_region_end; i++)
	//	{
	//		total_control_signal += floor(control_track_data[i]);	
	//	} // i loop.
	//} // Extension loop.

	//// Region not significant: There is not enough control reads.
	//if(total_control_signal < MIN_N_CONTROL_SIGNAL/scaling_factor)
	//{
	//	//fprintf(stderr, "Could not find enough control signal, returning p_val=1.\n");
	//	return(1.0);
	//}

	// Scale control signal.
	total_control_signal /= scaling_factor;

	// Normalize both signals with respect to normalization region, if they are larger than the region size.
	if(length_normalization && 
		(region_end - region_start+1) > l_norm_win)
	{
		total_profile_signal = total_profile_signal / ((double)(region_end - region_start+1) / l_norm_win);
		total_control_signal = total_control_signal / ((double)(region_end - region_start+1) / l_norm_win);

		if(floor(total_control_signal) == 0)
		{
			total_control_signal = 1.0;
		}

		// The total signal hit 0, return 1.
		if(floor(total_profile_signal) == 0)
		{
			//fprintf(stderr, "Hit 0.\n");
			// Return the most insignificant p-value.
			return(xlog(1.0));
		}
	}
	
	int grand_total_int = (int)(total_profile_signal + total_control_signal);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if(_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int+3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(.5);
	double cur_region_log_p_val = xlog(0.0);
	for(int i = 0; i <= total_control_signal; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Free memory if it is allocated in the function.
	if(_log_factorials == NULL)
	{
		delete [] log_factorials;
	}

	return(cur_region_log_p_val);
}

double get_modified_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, 
	double flip_prob, 
	int region_start, int region_end, 
	double* _log_factorials, 
	double scaling_factor,
	bool length_normalization, int l_norm_win)
{
	// Get the p-value for the current region.
	double total_profile_signal = 0;
	double total_control_signal = 0;
	for(int i = region_start; i <= region_end; i++)
	{
		total_profile_signal += floor(signal_track_data[i]);
		total_control_signal += floor(control_track_data[i]);	
	} // i loop.

	total_profile_signal /= scaling_factor;

	// Scale control signal.
	total_control_signal /= scaling_factor;

	// Normalize both signals with respect to normalization region, if they are larger than the region size.
	if(length_normalization && 
		(region_end - region_start+1) > l_norm_win)
	{
		total_profile_signal = total_profile_signal / ((double)(region_end - region_start+1) / l_norm_win);
		total_control_signal = total_control_signal / ((double)(region_end - region_start+1) / l_norm_win);

		if(floor(total_control_signal) == 0)
		{
			total_control_signal = 1.0;
		}

		// The total signal hit 0, return 1.
		if(floor(total_profile_signal) == 0)
		{
			//fprintf(stderr, "Hit 0.\n");
			// Return the most insignificant p-value.
			return(xlog(1.0));
		}
	}
	
	int grand_total_int = (int)(total_profile_signal + total_control_signal);

	// Setup the log factorials.
	double* log_factorials = NULL;
	if(_log_factorials == NULL)
	{
		log_factorials = buffer_log_factorials(grand_total_int+3);
	}
	else
	{
		log_factorials = _log_factorials;
	}

	// Compute the p-values.
	double log_flip = log(flip_prob);
	double cur_region_log_p_val = xlog(0.0);
	for(int i = 0; i <= total_control_signal; i++)
	{
		double log_cur_half_pow = log_flip * grand_total_int;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
		cur_region_log_p_val = xlog_sum(cur_region_log_p_val, xlog_mul(log_cur_perm, log_cur_half_pow));
	} // i loop.

	// Free memory if it is allocated in the function.
	if(_log_factorials == NULL)
	{
		delete [] log_factorials;
	}

	return(cur_region_log_p_val);
}

vector<t_annot_region*>* prune_region_ends_per_modified_binomial_p_value_minimization(vector<t_annot_region*>* regions, 
	double flip_prob,
	double* signal_profile, int l_profile, 
	double* control_profile, int l_control,
	int l_fragment)
{
	// List of pruned regions.
	vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();


	double max_grand_total_int = 0.0;
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		// Get the extrema signal values.
		double current_grand_total = 0.0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			if(i < l_profile)
			{
				current_grand_total += signal_profile[i];
			}

			if(i < l_control)
			{
				current_grand_total += control_profile[i];
			}
		} // i loop.

		if(max_grand_total_int < current_grand_total)
		{
			max_grand_total_int = (int)(floor(current_grand_total));
		} 
	} // i_reg loop.

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);

	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Pruning region ends for %ld regions via p-value minimization.\n", regions->size());

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<int>* pruned_starts_per_thresholds = new vector<int>();
	vector<int>* pruned_ends_per_thresholds = new vector<int>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Processing %d.(%ld) region: %d-%d\n", i_reg, regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		}

		// Skip the region if the signal is out of the control signal range.
		if(regions->at(i_reg)->end > l_control)
		{
			fprintf(stderr, "Skipping pruning %d-%d\n", regions->at(i_reg)->start, regions->at(i_reg)->end);
			continue;
		}

		// Get the extrema signal values.
		int max_sig = 0;
		int min_sig = 1000*1000;
		double total_signal = 0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			int cur_floorized_val = (int)(floor(signal_profile[i]));
			if(max_sig < cur_floorized_val)
			{
				max_sig = cur_floorized_val;
			}

			if(min_sig > cur_floorized_val)
			{
				min_sig = cur_floorized_val;
			}

			total_signal += signal_profile[i];
		} // i loop.

		// Limit the maximum at 400.
		max_sig = (max_sig > 500)?(500):(max_sig);

		// Start thresholding over all possible thresholds.			
		for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		{
			// Threshold the ends.
			int pruned_start = regions->at(i_reg)->start;
			int pruned_end = regions->at(i_reg)->end;

			bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					passed_threshold = true;
					pruned_start = i_pr;
					break;
				}
			} // i_pr loop.

			for(int i_pr = pruned_end; i_pr >= pruned_start; i_pr--)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					pruned_end = i_pr;
					break;
				}
			} // i_pr loop.

			// Get the binomial p-value for current region.
			//double cur_pruned_pvalue = get_binomial_pvalue_per_region(signal_profile, control_profile, pruned_start, pruned_end, log_factorials, l_fragment, false, 0);
			double cur_pruned_pvalue = get_modified_binomial_pvalue_per_region(signal_profile, control_profile, flip_prob, pruned_start, pruned_end, log_factorials, l_fragment, false, 0);

			// Add the current p-value.
			p_values_per_threshold->push_back(cur_pruned_pvalue);
			pruned_starts_per_thresholds->push_back(pruned_start);
			pruned_ends_per_thresholds->push_back(pruned_end);

			double total_pruned_sig = 0;
			for(int i_sig = pruned_start; i_sig <= pruned_end; i_sig++)
			{
				total_pruned_sig += signal_profile[i_sig];
			} // i_sig.

			// If we are below 50% of the total probability, stop since most of the signal is lost already.
			if(total_pruned_sig/total_signal <= 0.5)
			{
				break;
			}
		} // sig thresh loop.

		// If the condition is satisfied, we are done.
		if(p_values_per_threshold->size() > 0)
		{
			double min_p_val = p_values_per_threshold->at(0);
			int min_p_val_i = 0;
			for(int i_th = 0; i_th < (int)p_values_per_threshold->size(); i_th++)
			{
				if(min_p_val > p_values_per_threshold->at(i_th))
				{
					min_p_val = p_values_per_threshold->at(i_th);
					min_p_val_i = i_th;
				}
			} // total mass fraction check.

			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_peak->start = pruned_starts_per_thresholds->at(min_p_val_i);
			pruned_peak->end = pruned_ends_per_thresholds->at(min_p_val_i);

			pruned_regions->push_back(pruned_peak);
		}
		else
		{
			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_regions->push_back(pruned_peak);
		}

		// Clear the values.
		p_values_per_threshold->clear();
		pruned_starts_per_thresholds->clear();
		pruned_ends_per_thresholds->clear();
	} // i_reg loop.

	// Dump the bed file.
	//dump_BED(op_bed_fp, pruned_regions);

	return(pruned_regions);
} // prune_region_ends_per_modified_binomial_p_value_minimization


vector<t_annot_region*>* prune_region_ends_per_lost_peak_mass(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control)
{
	return(NULL);
}

vector<t_annot_region*>* prune_region_ends_per_signal_distribution(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control)
{
	return(NULL);
}

vector<t_annot_region*>* prune_region_ends_per_window_average(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control,
	double fdr, int l_win)
{
		// Go over all the windows.
		int cur_win_start = 0;
		int cur_l_win = l_win;

		// Sort the regions.
		sort(regions->begin(), regions->end(), sort_regions);

		// List of pruned regions.
		vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();

		// The windowing schedule matches up with the ms_decomposition.
		while(cur_win_start < l_profile)
		{
			// Set the length of the window 
			if(cur_win_start+l_win > l_profile)
			{
				cur_l_win = l_profile - cur_win_start;
			}
			else
			{
				cur_l_win = l_win;
			}

			//fprintf(stderr, "Processing %d-%d\n", cur_win_start, cur_win_start+cur_l_win);

			// Get the minimum and maximum in the window.
			int max_thresh = 0;
			for(int i = cur_win_start; i < cur_win_start+cur_l_win && i < l_profile; i++)
			{
				if(max_thresh < signal_profile[i])
				{
					max_thresh = (int)(signal_profile[i]);
				}
			} // i loop.


//#define __POISSON_BCKGRND__
//#undef __GEOMETRIC_BCKGRND__
#ifdef __POISSON_BCKGRND__
			double cur_win_thresh = get_poisson_thresholds_per_win(signal_profile,
											l_profile,
											cur_l_win,
											cur_win_start,
											fdr,
											1,
											max_thresh);
#elif defined __GEOMETRIC_BCKGRND__
			double cur_win_thresh = get_geometric_thresholds_per_win(signal_profile,
											l_profile,
											cur_l_win,
											cur_win_start,
											fdr,
											1,
											max_thresh);
#else
			#error "Did not select a background model."
#endif // Bckgrnd selection

			//fprintf(stderr, "Window: %d-%d: Threshold: %lf\n", cur_win_start, cur_win_start + cur_l_win, cur_win_thresh);

			// Go over all the regions in the current window and prune the ends with respect to the threshold.
			for(int i_reg = 0; cur_win_thresh > 0 && i_reg < (int)regions->size(); i_reg++)
			{
				// Process the current region if its start is in the current window.
				if((regions->at(i_reg)->start >= cur_win_start &&
					regions->at(i_reg)->start < cur_win_start + cur_l_win))
				{
					// For the current region, update the ends.
					int pruned_start = regions->at(i_reg)->start;
					int pruned_end = regions->at(i_reg)->end;

					bool passed_threshold = false;
					for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
					{
						if(signal_profile[i_pr] > cur_win_thresh)
						{
							passed_threshold = true;
							pruned_start = i_pr;
							break;
						}
					} // i_pr loop.

					for(int i_pr = pruned_end; i_pr >= pruned_start; i_pr--)
					{
						if(signal_profile[i_pr] > cur_win_thresh)
						{
							pruned_end = i_pr;
							break;
						}
					} // i_pr loop.

					// If this peak passed the threshold, add it to the list of all peaks.
					if(passed_threshold)
					{
						if(pruned_start <= pruned_end)
						{
							t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
							pruned_peak->start = pruned_start;
							pruned_peak->end = pruned_end;

							pruned_regions->push_back(pruned_peak);
						}
						else
						{
							fprintf(stderr, "Pruning error: %d-%d\n", pruned_start, pruned_end);
							exit(0);
						}
					}
				} // Check if this is a region in the current window.
				else
				{
				}
			} // i_reg loop.

			// Update window.
			cur_win_start += l_win;
		} // window start check.

	//	// Dump the bed file.
	//	dump_BED(op_bed_fp, pruned_regions);
	//} // -prune_region_ends_per_window_average

	return(pruned_regions);
}

#define _ADAPTIVE_THRESHOLDING_DIVIDER_ (4)
vector<t_annot_region*>* prune_region_ends_per_adaptive_thresholding(vector<t_annot_region*>* regions, double* signal_track_data, int l_profile_signal, 
	double* control_track_data, int l_control_signal,
	int l_fragment)
{
	fprintf(stderr, "Adaptive thresholding based pruning with threshold %d\n", _ADAPTIVE_THRESHOLDING_DIVIDER_);
	vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();

	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Processing %d.(%ld) region: %d-%d\n", i_reg, regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		}

		// Skip the region if the signal is out of the control signal range.
		if(regions->at(i_reg)->end > l_control_signal)
		{
			fprintf(stderr, "Skipping pruning %d-%d\n", regions->at(i_reg)->start, regions->at(i_reg)->end);
			continue;
		}

		double total_signal = 0;
		double max_sig = 0;
		double min_sig = 10*1000*1000;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			total_signal += signal_track_data[i];
			if(max_sig < signal_track_data[i])
			{
				max_sig = signal_track_data[i];
			}

			if(min_sig > signal_track_data[i])
			{
				min_sig = signal_track_data[i];
			}
		} // i loop.

		// Find the threshold that is above 50% of the mass.
		double cur_thresh = (min_sig + max_sig) / _ADAPTIVE_THRESHOLDING_DIVIDER_;

		int pruned_start = regions->at(i_reg)->start;
		int pruned_end = regions->at(i_reg)->end;
		while(1)
		{
			for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
			{
				if(signal_track_data[i] > cur_thresh)
				{
					pruned_start = i;
					break;
				}
			} // i loop.

			for(int i = regions->at(i_reg)->end; i >= regions->at(i_reg)->start; i--)
			{
				if(signal_track_data[i] > cur_thresh)
				{
					pruned_end = i;
					break;
				}
			} // i loop.

			double total_pruned_sig = 0;
			for(int i = pruned_start; i <= pruned_end; i++)
			{
				total_pruned_sig += signal_track_data[i];
			} // i loop.

			if(total_pruned_sig/total_signal > 0.5)
			{
				break;
			}
			else
			{
				cur_thresh--;
			}
		} // mass checking loop.

		if(pruned_end > pruned_start)
		{
			t_annot_region* pruned_region = duplicate_region(regions->at(i_reg));
			pruned_region->start = pruned_start;
			pruned_region->end = pruned_end;
			pruned_regions->push_back(pruned_region);
		}
	} // i_reg loop.

	return(pruned_regions);
}

vector<t_annot_region*>* prune_region_ends_per_p_value_minimization_per_window(vector<t_annot_region*>* regions, double* signal_track_data, int l_profile_signal, 
	double* control_track_data, int l_control_signal,
	int l_fragment)
{
	vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();

	double max_grand_total_int = 0.0;
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		// Get the extrema signal values.
		double current_grand_total = 0.0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			if(i < l_profile_signal)
			{
				current_grand_total += signal_track_data[i];
			}

			if(i < l_control_signal)
			{
				current_grand_total += control_track_data[i];
			}
		} // i loop.

		if(max_grand_total_int < current_grand_total)
		{
			max_grand_total_int = (int)(floor(current_grand_total));
		} 
	} // i_reg loop.

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);

	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Pruning region ends for %ld regions via p-value minimization.\n", regions->size());

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<int>* pruned_starts_per_thresholds = new vector<int>();
	vector<int>* pruned_ends_per_thresholds = new vector<int>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Processing %d.(%ld) region: %d-%d\n", i_reg, regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		}

		// Skip the region if the signal is out of the control signal range.
		if(regions->at(i_reg)->end > l_control_signal)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Skipping pruning %d-%d\n", regions->at(i_reg)->start, regions->at(i_reg)->end);
			continue;
		}

		double total_signal = 0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			total_signal += signal_track_data[i];
		} // i loop.

		int l_win = (regions->at(i_reg)->end - regions->at(i_reg)->start) / 40;
		if(l_win < 100)
		{
			l_win = 100;
		}

		// Start thresholding over all possible thresholds.			
		//for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		for(int cur_start = regions->at(i_reg)->start; cur_start <= regions->at(i_reg)->end; cur_start += l_win)
		{
			double cur_total_sig_till_end = 0.0;
			for(int i_sig = cur_start; i_sig <= regions->at(i_reg)->end; i_sig++)
			{
				cur_total_sig_till_end += signal_track_data[i_sig];
			} // i_sig.

			// If we do not have enough signal till the end, just break.
			if(cur_total_sig_till_end / total_signal < 0.5)
			{
				break;
			}

			for(int cur_end = cur_start + l_win; cur_end < regions->at(i_reg)->end; cur_end += l_win)
			{
				// Add the current p-value.
				double total_pruned_sig = 0;
				for(int i_sig = cur_start; i_sig <= cur_end; i_sig++)
				{
					total_pruned_sig += signal_track_data[i_sig];
				} // i_sig.

				// If we are below 50% of the total probability, stop since most of the signal is lost already.
				if(total_pruned_sig/total_signal > 0.5)
				{
					// Get the binomial p-value for current region.
					double cur_pruned_pvalue = get_binomial_pvalue_per_region(signal_track_data, control_track_data, cur_start, cur_end, log_factorials, l_fragment, false, 0);

					p_values_per_threshold->push_back(cur_pruned_pvalue);
					pruned_starts_per_thresholds->push_back(cur_start);
					pruned_ends_per_thresholds->push_back(cur_end);
				}
			} // cur_end loop.
		} // sig thresh loop.

		// If the condition is satisfied, we are done.
		if(p_values_per_threshold->size() > 0)
		{
			double min_p_val = p_values_per_threshold->at(0);
			int min_p_val_i = 0;
			for(int i_th = 0; i_th < (int)p_values_per_threshold->size(); i_th++)
			{
				if(min_p_val > p_values_per_threshold->at(i_th))
				{
					min_p_val = p_values_per_threshold->at(i_th);
					min_p_val_i = i_th;
				}
			} // total mass fraction check.

			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_peak->start = pruned_starts_per_thresholds->at(min_p_val_i);
			pruned_peak->end = pruned_ends_per_thresholds->at(min_p_val_i);

			pruned_regions->push_back(pruned_peak);
		}
		else
		{
			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_regions->push_back(pruned_peak);
		}

		// Clear the values.
		p_values_per_threshold->clear();
		pruned_starts_per_thresholds->clear();
		pruned_ends_per_thresholds->clear();
	} // i_reg loop.
	delete [] log_factorials;

	return(pruned_regions);
} // prune_region_ends_per_p_value_minimization_per_window

vector<t_annot_region*>* prune_region_ends_per_p_value_minimization(vector<t_annot_region*>* regions, 
	double* signal_profile, int l_profile, 
	double* control_profile, int l_control,
	int l_fragment)
{
	// List of pruned regions.
	vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();


	double max_grand_total_int = 0.0;
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		// Get the extrema signal values.
		double current_grand_total = 0.0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			if(i < l_profile)
			{
				current_grand_total += signal_profile[i];
			}

			if(i < l_control)
			{
				current_grand_total += control_profile[i];
			}
		} // i loop.

		if(max_grand_total_int < current_grand_total)
		{
			max_grand_total_int = (int)(floor(current_grand_total));
		} 
	} // i_reg loop.

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);
	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Pruning region ends for %ld regions via p-value minimization.\n", regions->size());

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<int>* pruned_starts_per_thresholds = new vector<int>();
	vector<int>* pruned_ends_per_thresholds = new vector<int>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Processing %d.(%ld) region: %d-%d\n", i_reg, regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		}

		// Skip the region if the signal is out of the control signal range.
		if(regions->at(i_reg)->end > l_control)
		{
			fprintf(stderr, "Skipping pruning %d-%d\n", regions->at(i_reg)->start, regions->at(i_reg)->end);
			continue;
		}

		// Get the extrema signal values.
		int max_sig = 0;
		int min_sig = 1000*1000;
		double total_signal = 0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			int cur_floorized_val = (int)(floor(signal_profile[i]));
			if(max_sig < cur_floorized_val)
			{
				max_sig = cur_floorized_val;
			}

			if(min_sig > cur_floorized_val)
			{
				min_sig = cur_floorized_val;
			}

			total_signal += signal_profile[i];
		} // i loop.

		// Limit the maximum at 400.
		max_sig = (max_sig > 500)?(500):(max_sig);

		// Start thresholding over all possible thresholds.			
		for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		{
			// Threshold the ends.
			int pruned_start = regions->at(i_reg)->start;
			int pruned_end = regions->at(i_reg)->end;

			bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					passed_threshold = true;
					pruned_start = i_pr;
					break;
				}
			} // i_pr loop.

			for(int i_pr = pruned_end; i_pr >= pruned_start; i_pr--)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					pruned_end = i_pr;
					break;
				}
			} // i_pr loop.

			// Get the binomial p-value for current region.
			double cur_pruned_pvalue = get_binomial_pvalue_per_region(signal_profile, control_profile, pruned_start, pruned_end, log_factorials, l_fragment, false, 0);

			// Add the current p-value.
			p_values_per_threshold->push_back(cur_pruned_pvalue);
			pruned_starts_per_thresholds->push_back(pruned_start);
			pruned_ends_per_thresholds->push_back(pruned_end);

			double total_pruned_sig = 0;
			for(int i_sig = pruned_start; i_sig <= pruned_end; i_sig++)
			{
				total_pruned_sig += signal_profile[i_sig];
			} // i_sig.

			// If we are below 50% of the total probability, stop since most of the signal is lost already.
			if(total_pruned_sig/total_signal <= 0.5)
			{
				break;
			}
		} // sig thresh loop.

		// If the condition is satisfied, we are done.
		if(p_values_per_threshold->size() > 0)
		{
			double min_p_val = p_values_per_threshold->at(0);
			int min_p_val_i = 0;
			for(int i_th = 0; i_th < (int)p_values_per_threshold->size(); i_th++)
			{
				if(min_p_val > p_values_per_threshold->at(i_th))
				{
					min_p_val = p_values_per_threshold->at(i_th);
					min_p_val_i = i_th;
				}
			} // total mass fraction check.

			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_peak->start = pruned_starts_per_thresholds->at(min_p_val_i);
			pruned_peak->end = pruned_ends_per_thresholds->at(min_p_val_i);

			pruned_regions->push_back(pruned_peak);
		}
		else
		{
			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_regions->push_back(pruned_peak);
		}

		// Clear the values.
		p_values_per_threshold->clear();
		pruned_starts_per_thresholds->clear();
		pruned_ends_per_thresholds->clear();
	} // i_reg loop.

	// Dump the bed file.
	//dump_BED(op_bed_fp, pruned_regions);

	return(pruned_regions);
} // prune_region_ends_per_p_value_minimization

vector<t_annot_region*>* prune_region_ends_per_p_value_minimization_no_control(vector<t_annot_region*>* regions, double* signal_track_data, int l_profile, 
	double scaling_factor,
	double per_win_bckgrnd_win_val)
{
	//return(NULL);
	//char* signal_track_data_fp = argv[2];
	//char* regions_bed_fp = argv[3];
	//int l_bckgrnd_norm_win = atoi(argv[4]);
	//char* op_bed_fp = argv[5];

	// Load the regions.
	//vector<t_annot_region*>* regions = load_BED(regions_bed_fp);
	//fprintf(stderr, "Loaded %d regions.\n", regions->size());

	// List of pruned regions.
	vector<t_annot_region*>* pruned_regions = new vector<t_annot_region*>();

	//// Load the profile signal.
	//fprintf(stderr, "Loading signal.\n");
	//int l_profile_signal = 0;
	//double* _cur_track_data = load_per_nucleotide_binary_profile(signal_track_data_fp, l_profile_signal);
	//double* signal_track_data = get_zero_indexed_per_one_indexed_data(_cur_track_data, l_profile_signal);
	//delete [] _cur_track_data;

	//double per_win_bckgrnd_win_val = get_fdr_based_bckgrnd_window_value(signal_track_data, l_profile_signal,
	//																				0.20, // Target FDR for setting the background window value.
	//																				1000*1000, // Thresholding window length.
	//																				l_bckgrnd_norm_win,
	//																				log(.05),
	//																				200);

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Pruning region ends for %ld regions via p-value minimization without control.\n", regions->size());

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<int>* pruned_starts_per_thresholds = new vector<int>();
	vector<int>* pruned_ends_per_thresholds = new vector<int>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
			fprintf(stderr, "Processing %d.(%ld) region: %d-%d\n", i_reg, regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		}

		// Get the extrema signal values.
		int max_sig = 0;
		int min_sig = 1000*1000;
		double total_signal = 0;
		for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
		{
			int cur_floorized_val = (int)(floor(signal_track_data[i]));
			if(max_sig < cur_floorized_val)
			{
				max_sig = cur_floorized_val;
			}

			if(min_sig > cur_floorized_val)
			{
				min_sig = cur_floorized_val;
			}

			total_signal += signal_track_data[i];
		} // i loop.

		//// We can in principle compute the upper bound on the factorials to be requested.
		//double max_total_signal = (total_signal/200 + cur_region_median_bckgrnd_sig * (double)((regions->at(i_reg)->end - regions->at(i_reg)->start + 1)/l_bckgrnd_norm_win)/200);

		//double* cur_reg_log_factorials = buffer_log_factorials((int)floor(max_total_signal+10));

		// Limit the maximum at 400.
		max_sig = (max_sig > 500)?(500):(max_sig);

		// Start thresholding over all possible thresholds.			
		for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		{
			// Threshold the ends.
			int pruned_start = regions->at(i_reg)->start;
			int pruned_end = regions->at(i_reg)->end;

			bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_track_data[i_pr] > sig_thresh)
				{
					passed_threshold = true;
					pruned_start = i_pr;
					break;
				}
			} // i_pr loop.

			for(int i_pr = pruned_end; i_pr >= pruned_start; i_pr--)
			{
				if(signal_track_data[i_pr] > sig_thresh)
				{
					pruned_end = i_pr;
					break;
				}
			} // i_pr loop.

			// Get the binomial p-value for current region.
			double cur_pruned_pvalue = get_binomial_pvalue_per_region_no_control(signal_track_data, l_profile, 
				pruned_start, pruned_end,
				per_win_bckgrnd_win_val,
				NULL, scaling_factor,
				1,
				true); // Do not normalize.

			// Add the current p-value.
			p_values_per_threshold->push_back(cur_pruned_pvalue);
			pruned_starts_per_thresholds->push_back(pruned_start);
			pruned_ends_per_thresholds->push_back(pruned_end);

			double total_pruned_sig = 0;
			for(int i_sig = pruned_start; i_sig <= pruned_end; i_sig++)
			{
				total_pruned_sig += signal_track_data[i_sig];
			} // i_sig.

			// If we are below 50% of the total probability, stop since most of the signal is lost already.
			if(total_pruned_sig/total_signal <= 0.5)
			{
				break;
			}
		} // sig thresh loop.

		// If the condition is satisfied, we are done.
		if(p_values_per_threshold->size() > 0)
		{
			double min_p_val = p_values_per_threshold->at(0);
			int min_p_val_i = 0;
			for(int i_th = 0; i_th < (int)p_values_per_threshold->size(); i_th++)
			{
				if(min_p_val > p_values_per_threshold->at(i_th))
				{
					min_p_val = p_values_per_threshold->at(i_th);
					min_p_val_i = i_th;
				}
			} // total mass fraction check.

			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_peak->start = pruned_starts_per_thresholds->at(min_p_val_i);
			pruned_peak->end = pruned_ends_per_thresholds->at(min_p_val_i);

			pruned_regions->push_back(pruned_peak);
		}
		else
		{
			t_annot_region* pruned_peak = duplicate_region(regions->at(i_reg));
			pruned_regions->push_back(pruned_peak);
		}

		// Clear the values.
		p_values_per_threshold->clear();
		pruned_starts_per_thresholds->clear();
		pruned_ends_per_thresholds->clear();
	} // i_reg loop.

	// Return the list of regions with pruned ends.
	return(pruned_regions);
}

