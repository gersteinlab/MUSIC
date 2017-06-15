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
#include "ms_filter_utils.h"
#include "ms_mapped_read_tools.h"
#include "ms_profile_normalization.h"
#include "ms_combinatorics.h"
#include "ms_utils.h"
#include "ms_signal_enrichment_utils.h"

bool __DUMP_SIGNAL_ENRICHMENT_MSGS__ = false;

/*
Start: sensitivity_satisfied=false && target_fpr_satisfied=true
1. sensitivity_satisfied=false && target_fpr_satisfied=false: Wait till sensitivity is satisfied or absolute max fpr is reached.
1. sensitivity_satisfied=true && target_fpr_satisfied=true: Wait will fpr is not satisfied or return if min sensitivity per fpr satisfying sensitivity is observed.
1.5. sensitivity >= 0.99 && target_fpr_satisfied=true::Early stop
1.5. sensitivity_satisfied=true && fpr < 0.00001: Early stop
1.6. fpr >= 0.05: Late stop with not good signal enrichment.
Normal end: sensitivity_satisfied=true && target_fpr_satisfied=false::We are done.
*/
int select_l_p_per_stats_file(const char* l_p_param_stats_fp, 
	double max_target_fpr, double min_target_sensitivity,
	double min_sensitivity_per_satisfying_fpr, // When FPR is satisfied by sensitivity isnt. (Early stopping)
	double max_fpr_per_satisfying_sensitivity, // When sensitivity is satisfied by FPR isnt. (Early stopping)
	double absolute_max_fpr) // This value enables controlling the FPR if sensitivity could not be satisfied for target FPR. (Late stopping)
{
	// This is a sanity check on the early stop conditions: They must be more stringent.
	if(min_sensitivity_per_satisfying_fpr < min_target_sensitivity ||
		max_fpr_per_satisfying_sensitivity > max_target_fpr)
	{
		fprintf(stderr, "Early stop conditions on sensitivity and FPR must be more stringent than target values.\n");
		exit(0);
	}

	FILE* f_l_p_param_stats = open_f("l_p_param_stats.txt", "r"); // This is the file to be read by MUSIC to select the l_p parameter.
	int satisfying_l_win = -1;
	while(1)
	{
/*
		double fc_FNR =			n_insig_pval_sig_FC_wins_per_l_win[l_win] / n_sig_FC_wins_per_l_win[l_win];
		double p_val_FNR =		n_insig_pval_sig_FC_wins_per_l_win[l_win] / n_insig_pval_wins_per_l_win[l_win];
		double fc_FPR =			n_sig_pval_insig_FC_wins_per_l_win[l_win] / n_insig_FC_wins_per_l_win[l_win];
		double p_val_FPR =		n_sig_pval_insig_FC_wins_per_l_win[l_win] / n_sig_pval_wins_per_l_win[l_win];
		double sensitivity =	n_sig_pval_sig_FC_wins_per_l_win[l_win] / n_sig_FC_wins_per_l_win[l_win];			

		fprintf(f_l_p_param_stats, "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", l_win,
			fc_FNR, p_val_FNR, 
			fc_FPR, p_val_FPR,
			sensitivity);
*/
		char cur_line[1000];
		if(fgets(cur_line, 1000, f_l_p_param_stats) == NULL)
		{
			break;
		}

		int l_win;
		double fc_FNR;
		double p_val_FNR;
		double fc_FPR;
		double p_val_FPR;
		double sensitivity;

		if(sscanf(cur_line, "%d %lf %lf %lf %lf %lf", &l_win, &fc_FNR, &p_val_FNR, &fc_FPR, &p_val_FPR, &sensitivity) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", cur_line);
			exit(0);
		}

		// Update window length.
		satisfying_l_win = l_win;

		bool target_fpr_satisfied = (fc_FPR <= max_target_fpr && p_val_FPR <= max_target_fpr);
		bool sensitivity_satisfied = (sensitivity >= min_target_sensitivity);

		// Return if we already passed the max FPR that is allowed no matter what the sensitivity is.
		if(fc_FPR > absolute_max_fpr || 
			p_val_FPR > absolute_max_fpr)
		{
			break;
		}

		if(target_fpr_satisfied || sensitivity_satisfied)
		{
			// Early select 1: Sensitivity reached a very high value before we hit the FPR check.
			if(target_fpr_satisfied && sensitivity_satisfied &&
				sensitivity >= min_sensitivity_per_satisfying_fpr)
			{
				break;
			}

			// Early select 2: Sensitivity is satisfied and FPR is very low, stop increasing window.
			if(sensitivity_satisfied && target_fpr_satisfied && 
				fc_FPR <= max_fpr_per_satisfying_sensitivity && 
				p_val_FPR <= max_fpr_per_satisfying_sensitivity)
			{
				break;
			}

			// This is the normal exit scenario, FPR is not satisfied at this value and sensitivity is satisfied.
			if(sensitivity_satisfied && !target_fpr_satisfied)
			{
				break;
			}
		}
		else 
		{
			// None is satisfied; increase window length.
		}
	} // file reading loop.
	fclose(f_l_p_param_stats);	

	// Return the window length.
	return(satisfying_l_win);
}

int select_l_p_per_stats_file(const char* l_p_param_stats_fp, double target_max_p_val_fc_fpr, double target_min_sensitivity, 
	double min_sensitivity_per_satisfying_FPR)
{
	/******************************************************************
	TODO::Must make sure that this function returns a reasonable 
	window length, always. It is guaranteed to return a window length,
	but not a reasonable one.
	******************************************************************/

	FILE* f_l_p_param_stats = open_f("l_p_param_stats.txt", "r"); // This is the file to be read by MUSIC to select the l_p parameter.
	int satisfying_l_win = -1;
	while(1)
	{
/*
		double fc_FNR =			n_insig_pval_sig_FC_wins_per_l_win[l_win] / n_sig_FC_wins_per_l_win[l_win];
		double p_val_FNR =		n_insig_pval_sig_FC_wins_per_l_win[l_win] / n_insig_pval_wins_per_l_win[l_win];
		double fc_FPR =			n_sig_pval_insig_FC_wins_per_l_win[l_win] / n_insig_FC_wins_per_l_win[l_win];
		double p_val_FPR =		n_sig_pval_insig_FC_wins_per_l_win[l_win] / n_sig_pval_wins_per_l_win[l_win];
		double sensitivity =	n_sig_pval_sig_FC_wins_per_l_win[l_win] / n_sig_FC_wins_per_l_win[l_win];			

		fprintf(f_l_p_param_stats, "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", l_win,
			fc_FNR, p_val_FNR, 
			fc_FPR, p_val_FPR,
			sensitivity);
*/
		char cur_line[1000];
		if(fgets(cur_line, 1000, f_l_p_param_stats) == NULL)
		{
			break;
		}

		int l_win;
		double fc_FNR;
		double p_val_FNR;
		double fc_FPR;
		double p_val_FPR;
		double sensitivity;

		if(sscanf(cur_line, "%d %lf %lf %lf %lf %lf", &l_win, &fc_FNR, &p_val_FNR, &fc_FPR, &p_val_FPR, &sensitivity) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", cur_line);
			exit(0);
		}

		// Moving on increases sensitivity and also FPRs. Following checks whether we should continue to increase l_win at the 
		// expense of increasing FPR.
		satisfying_l_win = l_win;

		// If we are already at very high sensitivity with very low FPR's then return the window.
		if(sensitivity >= min_sensitivity_per_satisfying_FPR && 
			fc_FPR < target_max_p_val_fc_fpr && 
			p_val_FPR < target_max_p_val_fc_fpr)
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
{
			fprintf(stderr, "Selecting l_win = %d @ FPR: %.3f, %.3f; Sensitivity: %.3f\n", l_win, fc_FPR, p_val_FPR, sensitivity);
}

			break;
		}

		// If we are over the maximum FPR level and we have the minimum sensitivity, return. This may have high FPR since we are 
		// enforcing sensitivity requirement. 
		if((fc_FPR > target_max_p_val_fc_fpr || p_val_FPR > target_max_p_val_fc_fpr) && 
			(sensitivity > target_min_sensitivity))
		{
if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
{
			fprintf(stderr, "Selecting l_win = %d @ FPR: %.3f, %.3f; Sensitivity: %.3f\n", l_win, fc_FPR, p_val_FPR, sensitivity);
}

			break;
		}
	} // file reading loop.
	fclose(f_l_p_param_stats);	

	// Return the window length.
	return(satisfying_l_win);
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
        if(i > 0 && i <= l_profile)
        {
			total_sig_in_cur_window += signal_profile[i];
		}
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
		//cur_win_total += floor(control_track_data[i]);
		cur_win_total += control_track_data[i];
	} // i loop.

	total_control_signal = cur_win_total;

	// Following loop computes the per window total control signals.
	int cur_win_start = ext_region_start;
	while(1)
	{
		// Remove the value that just left the window.
		//cur_win_total -= floor(control_track_data[cur_win_start]);
		cur_win_total -= control_track_data[cur_win_start];

		// Add the newly entering value to the window.
		//cur_win_total += floor(control_track_data[cur_win_start+l_norm_win]);
		cur_win_total += control_track_data[cur_win_start+l_norm_win];

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
			double cur_control_value = control_track_data[i];
			//cur_total_control_signal += floor(cur_control_value);
			cur_total_control_signal += (cur_control_value);
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
		total_control_signal += (control_track_data[i]);
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

		// Compute the current permutation.
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

void tune_regions_per_window_average_at_boundaries(vector<t_annot_region*>* regions, double* signal_profile, int l_profile,
	double fdr, int l_win)
{
	// Go over all the windows.
	int cur_win_start = 1;
	int cur_l_win = l_win;

	// Sort the regions.
	sort(regions->begin(), regions->end(), sort_regions);

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

		FILE* f_tuned_feat_ends = NULL;
			
		if(__DUMP_SIGNAL_ENRICHMENT_MSGS__)
		{
			f_tuned_feat_ends = open_f("tuned_feat_ends.bed", "a");
		}

		// Go over all the regions in the current window and prune the ends with respect to the threshold.
		for(int i_reg = 0; cur_win_thresh > 0 && i_reg < (int)regions->size(); i_reg++)
		{
			// Process the current region if its start is in the current window.
			if((regions->at(i_reg)->start >= cur_win_start &&
				regions->at(i_reg)->start < cur_win_start + cur_l_win))
			{
				bool reg_tuned = false;

				// Following is actually a very weak and unstable check, we must evaluate whether there is high signal around the start/end of this region.
				if(signal_profile[regions->at(i_reg)->start] >= cur_win_thresh)
				{
					regions->at(i_reg)->start -= 100;
					regions->at(i_reg)->start = (regions->at(i_reg)->start < 1)?(1):(regions->at(i_reg)->start);					

					reg_tuned = true;
				}
				else if(signal_profile[regions->at(i_reg)->end] >= cur_win_thresh)
				{
					regions->at(i_reg)->end += 100;
					regions->at(i_reg)->end = (regions->at(i_reg)->end > l_profile)?(l_profile):(regions->at(i_reg)->end);

					reg_tuned = true;
				}

				// If the region is tuned, dump the region, if we are dumping.
				if(f_tuned_feat_ends != NULL &&
					reg_tuned)
				{
					fprintf(f_tuned_feat_ends, "%s\t%d\t%d\t%.3f %.3f %.3f\t.\t+\n", 
							regions->at(i_reg)->chrom, regions->at(i_reg)->start, regions->at(i_reg)->end,
							signal_profile[regions->at(i_reg)->start], signal_profile[regions->at(i_reg)->end], cur_win_thresh);
				}
			}
		} // i_reg loop.

		if(f_tuned_feat_ends != NULL)
		{
			fclose(f_tuned_feat_ends);
		}

		// Update window.
		cur_win_start += l_win;
	} // window start check.
}

vector<t_annot_region*>* prune_region_ends_per_window_average(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control,
	double fdr, int l_win)
{
		// Go over all the windows.
		int cur_win_start = 1;
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
						if(i_pr < l_profile &&
							i_pr > 0 &&
							signal_profile[i_pr] > cur_win_thresh)
						{
							passed_threshold = true;
							pruned_start = i_pr;
							break;
						}
					} // i_pr loop.

					for(int i_pr = pruned_end; i_pr >= pruned_start; i_pr--)
					{
						if(i_pr < l_profile &&
							i_pr > 0 &&
							signal_profile[i_pr] > cur_win_thresh)
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
			//fprintf(stderr, "Skipping pruning %d-%d\n", regions->at(i_reg)->start, regions->at(i_reg)->end);
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

		// Limit the maximum at 500.
		max_sig = (max_sig > 500)?(500):(max_sig);

		// Start thresholding over all possible thresholds.
		int prev_pr_s = -1;
		int prev_pr_e = -1;
		for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		{
			// Threshold the ends.
			int pruned_start = regions->at(i_reg)->start;
			int pruned_end = regions->at(i_reg)->end;

			//bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					//passed_threshold = true;
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

            if(prev_pr_s == pruned_start && prev_pr_e == pruned_end)
            {
                    continue;
            }

            prev_pr_s = pruned_start;
            prev_pr_e = pruned_end;

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

			// TODO::Following should actually be removed since it creates issues with overextended ER ends.
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

	delete [] log_factorials;

	// Dump the bed file.
	//dump_BED(op_bed_fp, pruned_regions);

	return(pruned_regions);
} // prune_region_ends_per_p_value_minimization

vector<t_annot_region*>* prune_region_ends_per_p_value_minimization_buffered_binomial_p_vals(vector<t_annot_region*>* regions, 
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

	int buffered_binomial_bin_size = 500;
	double** buffered_binomial_p_vals = buffer_binomial_pvalues((max_grand_total_int / l_fragment) + 100, log_factorials, buffered_binomial_bin_size);

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

		// Limit the maximum at 500.
		max_sig = (max_sig > 500)?(500):(max_sig);

		// Start thresholding over all possible thresholds.			
		for(int sig_thresh = min_sig; sig_thresh <= max_sig; sig_thresh++)
		{
			// Threshold the ends.
			int pruned_start = regions->at(i_reg)->start;
			int pruned_end = regions->at(i_reg)->end;

			//bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					//passed_threshold = true;
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
			//double test_cur_pruned_pvalue = get_binomial_pvalue_per_region(signal_profile, control_profile, pruned_start, pruned_end, log_factorials, l_fragment, false, 0);

			// Get the p-value for the current region.
			double total_profile_signal = 0;
			double total_control_signal = 0;
			for(int i = pruned_start; i <= pruned_end; i++)
			{
				total_profile_signal += floor(signal_profile[i]);
				total_control_signal += (control_profile[i]);
			} // i loop.

			total_profile_signal /= l_fragment;
			total_control_signal /= l_fragment;
	
			int grand_total_int = (int)(total_profile_signal + total_control_signal);

			double log_flip = log(.5);

			int buffered_base_val = ((int)(total_control_signal / buffered_binomial_bin_size)) * buffered_binomial_bin_size;
			double cur_pruned_pvalue = buffered_binomial_p_vals[grand_total_int][buffered_base_val];
			for(int i = buffered_base_val+1; i <= total_control_signal; i++)
			{
				double log_cur_half_pow = log_flip * grand_total_int;
				double log_cur_perm = 0.0; // = xlog(1.0).

				// Compute the current permutation.
				log_cur_perm = xlog_div(log_factorials[grand_total_int], xlog_mul(log_factorials[i], log_factorials[grand_total_int-i]));
				cur_pruned_pvalue = xlog_sum(cur_pruned_pvalue, xlog_mul(log_cur_perm, log_cur_half_pow));
			} // i loop.

			//fprintf(stderr, "p-vals: %lf, %lf\n", test_cur_pruned_pvalue, cur_pruned_pvalue);

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

	// Free buffered value memory.
	delete [] log_factorials;
	delete_buffered_binomial_pvalues(buffered_binomial_p_vals, max_grand_total_int, buffered_binomial_bin_size);

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

			//bool passed_threshold = false;
			for(int i_pr = pruned_start; i_pr <= pruned_end; i_pr++)
			{
				if(signal_profile[i_pr] > sig_thresh)
				{
					//passed_threshold = true;
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

	delete [] log_factorials;

	// Dump the bed file.
	//dump_BED(op_bed_fp, pruned_regions);

	return(pruned_regions);
} // prune_region_ends_per_modified_binomial_p_value_minimization


