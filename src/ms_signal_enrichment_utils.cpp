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
#include "ms_min_max_utils.h"
#include "ms_gsl_fft_filter_utils.h"
#include "ms_mapped_read_tools.h"
#include "ms_profile_normalization.h"
#include "ms_combinatorics.h"
#include "ms_utils.h"
#include "ms_signal_enrichment_utils.h"

bool __DUMP_SIGNAL_ENRICHMENT_MSGS__ = false;

///*
//Here is the contingency table:
//
//			g1	g2
//cond1:		n11	n12
//cond2:		n21	n22
//
//The test asks whether there is a significant difference between two groups (g1 and g2) with respect to the observed conditions, i.e., the significance of the fractions: a/c and b/d .
//cur_peak_profile_signal/200, cur_subpeak_profile_signal/200,cur_peak_control_signal/200, cur_subpeak_control_signal/200
//
//cur_peak_profile_signal/200, cur_subpeak_profile_signal/200,
//cur_peak_control_signal/200, cur_subpeak_control_signal/200
//*/
//double compare_read_numbers_per_peaks(double p1_profile_signal, double p1_control_signal, double p2_profile_signal, double p2_control_signal, double normalization_factor)
//{
//	// Compute the more extreme enrichments:
//	double cur_p1_profile_signal = floor(p1_profile_signal/normalization_factor);
//	double cur_p1_control_signal = floor(p1_control_signal/normalization_factor);
//	double cur_p2_profile_signal = floor(p2_profile_signal/normalization_factor);
//	double cur_p2_control_signal = floor(p2_control_signal/normalization_factor);
//
//	// Following looks at the more extreme profile signal levels.
//	double p_val = xlog(0.0);
//	while(1)
//	{
//		p_val = xlog_sum(p_val, get_fishers_exact_p_val(cur_p1_profile_signal, cur_p2_profile_signal, cur_p1_control_signal, cur_p2_control_signal));
//
//		if(cur_p1_control_signal == 0 || cur_p2_control_signal == 0)
//		{
//			break;
//		}
//
//		cur_p1_profile_signal--;
//		cur_p2_profile_signal--;
//		cur_p1_control_signal++;
//		cur_p2_control_signal++;	
//	}
//
//	return(p_val);
//}
//
//double get_fishers_exact_p_val(double n11, double n12, double n21, double n22)
//{
//	int int_n11 = (int)(floor(n11));
//	int int_n12 = (int)(floor(n12));
//	int int_n21 = (int)(floor(n21));
//	int int_n22 = (int)(floor(n22));
//	int n = int_n11 + int_n12 + int_n21 + int_n22;
//
//	double* log_factorials = buffer_log_factorials(n+5);
//
//	double top1 = xlog_div(log_factorials[int_n11+int_n12], xlog_mul(log_factorials[int_n11], log_factorials[int_n12]));
//
//	double top2 = xlog_div(log_factorials[int_n21+int_n22], xlog_mul(log_factorials[int_n21], log_factorials[int_n22]));
//
//	double bottom = xlog_div(log_factorials[n], xlog_mul(log_factorials[int_n11], log_factorials[int_n21]));
//
//	delete [] log_factorials;
//	return(xlog_div(xlog_mul(top1, top2), bottom));
//}

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

// Following is to used for computing the p-values in punctate/TF peak calls where the p-value is computed from the largest window signal the region, which is not suitable for broad marks.
double get_binomial_pvalue_per_region_neighborhood_window_control_per_max_profile_signal(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win, int l_signal, int l_control, int n_vicinity_wins)
{
	// Get the total profile signal: This is computed by sliding a window over the region.
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

		// Update the normalization window length, this way there is no underpowering for the shorter regions.
		l_norm_win = region_end - region_start;
	}

	// Get the control signal: Get the total profile signal: This is computed by sliding a window over the vicinity extended region (in terms of length of the whole ER).
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

	// If we are past the end of the chromosome, flag this regions as insignificant.
	if(ext_region_start + l_norm_win > l_control)
	{
		return(log(1.0));
	}

	// Must make sure we do not go over the limits.
	double cur_win_total = 0.0;
	for(int i = ext_region_start; i < MIN(l_control, ext_region_start + l_norm_win); i++)
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

		if(cur_win_total > total_control_signal)
		{
			total_control_signal = cur_win_total;
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

vector<t_annot_region*>* filter_regions_per_window_average(vector<t_annot_region*>* regions, 
	double* signal_profile, int l_profile, double* control_profile, int l_control,
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
		vector<double>* cur_region_sigs = new vector<double>();

		// Go over all the regions in the current window and prune the ends with respect to the threshold.
		for(int i_reg = 0; cur_win_thresh > 0 && i_reg < (int)regions->size(); i_reg++)
		{
			// Process the current region if its start is in the current window.
			if((regions->at(i_reg)->start >= cur_win_start &&
				regions->at(i_reg)->start < cur_win_start + cur_l_win))
			{
				// Get the median signal value in the current region.
				for(int i = regions->at(i_reg)->start; i <= regions->at(i_reg)->end; i++)
				{
					cur_region_sigs->push_back(signal_profile[i]);
				} // i loop.

				sort(cur_region_sigs->begin(), cur_region_sigs->end());

				if(cur_region_sigs->at(cur_region_sigs->size() / 2) > cur_win_thresh)
				{
					pruned_regions->push_back(duplicate_region(regions->at(i_reg)));
				}

				cur_region_sigs->clear();
			} // check if the current region is in the current window.
		} // i_reg loop.

		delete cur_region_sigs;

		// Update window.
		cur_win_start += l_win;
	} // window start check.

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
	fprintf(stderr, "Pruned to %d (%d) regions.\n", (int)regions->size(), (int)pruned_regions->size());

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
{
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);
}

	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
{
	fprintf(stderr, "Pruning region ends for %d regions via p-value minimization.\n", (int)regions->size());
}

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<double>* pruned_starts_per_thresholds = new vector<double>();
	vector<double>* pruned_ends_per_thresholds = new vector<double>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		//if(i_reg % 1000 == 0)
		//{
		//	fprintf(stderr, "Processing %d.(%d) region: %d-%d\n", i_reg, (int)regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		//}

		// Skip the region if the signal is out of the control signal range.
		if(regions->at(i_reg)->end > l_control_signal)
		{
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
			pruned_peak->start = (int)pruned_starts_per_thresholds->at(min_p_val_i);
			pruned_peak->end = (int)pruned_ends_per_thresholds->at(min_p_val_i);

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
{
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);
}

	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
{
	fprintf(stderr, "Pruning region ends for %d regions via p-value minimization.\n", (int)regions->size());
}

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<double>* pruned_starts_per_thresholds = new vector<double>();
	vector<double>* pruned_ends_per_thresholds = new vector<double>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		//if(i_reg % 1000 == 0)
		//{
		//	fprintf(stderr, "Processing %d.(%d) region: %d-%d\n", i_reg, (int)regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		//}

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
			pruned_peak->start = (int)(pruned_starts_per_thresholds->at(min_p_val_i));
			pruned_peak->end = (int)(pruned_ends_per_thresholds->at(min_p_val_i));

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
{
	fprintf(stderr, "Setting up the log factorials with maximum grand total of %lf\n", max_grand_total_int);
}

	double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));

if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
{
	fprintf(stderr, "Pruning region ends for %d regions via p-value minimization.\n", (int)regions->size());
}

	// Go over all the regions.
	vector<double>* p_values_per_threshold = new vector<double>();
	vector<double>* pruned_starts_per_thresholds = new vector<double>();
	vector<double>* pruned_ends_per_thresholds = new vector<double>();
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		//if(i_reg % 1000 == 0)
		//{
		//	fprintf(stderr, "Processing %d.(%d) region: %d-%d\n", i_reg, (int)regions->size(), regions->at(i_reg)->start, regions->at(i_reg)->end);
		//}

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
		if((int)p_values_per_threshold->size() > 0)
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
			pruned_peak->start = (int)(pruned_starts_per_thresholds->at(min_p_val_i));
			pruned_peak->end = (int)(pruned_ends_per_thresholds->at(min_p_val_i));

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


