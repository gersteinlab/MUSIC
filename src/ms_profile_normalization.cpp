#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <string.h>
#include "ms_profile_normalization.h"

bool __DUMP_INPUT_NORMALIZATION_MSGS__ = false;

//double get_per_window_scaling_factor_for_profile1_per_profile2(int n_meg_wins,
//				vector<t_mapped_fragment*>* chip_seq_fragments,
//				vector<t_mapped_fragment*>* control_fragments,
//				double P_f, 
//				int enrichment_mapped_fragment_length)

void get_per_window_scaling_factors_for_profile1_per_profile2(double* signal_profile1, int l_prof1, double* signal_profile2, int l_prof2, int l_win,
	double& per_win_2DOF_lls_scaling_factor,
	double& per_win_1DOF_lls_scaling_factor,
	double& total_sig_scaling_factor)
{
	/*vector<int>* chip_seq_frag_cnt = new vector<int>();
	vector<int>* control_frag_cnt = new vector<int>();*/

	int l_processable_prof = (l_prof1 < l_prof2)?(l_prof1):(l_prof2);
	int n_wins_2_process = l_processable_prof / l_win;

	double* prof1_win_total_sigs = new double[n_wins_2_process];
	double* prof2_win_total_sigs = new double[n_wins_2_process];

	double total_prof1 = 0;
	double total_prof2 = 0;
	for(int i = 1; i <= l_prof1; i++)
	{
		total_prof1 += signal_profile1[i];
	} // i loop

	for(int i = 1; i <= l_prof2; i++)
	{
		total_prof2 += signal_profile2[i];
	} // i loop

	total_sig_scaling_factor = total_prof2 / total_prof1;

	// Count the numbers per window. 
	//int n_scaling_win = n_meg_wins * MEG_BASE / win_size;

	//int chip_seq_background_cnt = 0;
	//int control_cnt = 0;

	// Exploit the fact that the fragments are sorted with respect to base indices in the lists:
	// Keep track of the current fragment index for both fragment lists.
	//int cur_chr_fragment_list_base_i = 0;
	//int cur_chr_control_fragment_list_base_i = 0;

	int n_processed_wins = 0;
	for(int i_win = 0; i_win < n_wins_2_process; i_win++)
	{
		int cur_win_start = l_win * i_win;
		int cur_win_end = (i_win + 1) * l_win;

		double cur_total_prof1 = 0;
		double cur_total_prof2 = 0;
		for(int i_nuc = cur_win_start; i_nuc < cur_win_end; i_nuc++)
		{
			cur_total_prof1 += signal_profile1[i_nuc];
			cur_total_prof2 += signal_profile2[i_nuc];
		} // i_nuc loop.

		if(cur_total_prof1 > 0 && cur_total_prof2 > 0)
		{
			prof1_win_total_sigs[n_processed_wins] = cur_total_prof1;
			prof2_win_total_sigs[n_processed_wins] = cur_total_prof2;
			n_processed_wins++;
		}
	} // i_win loop.

     // Compute the slopes.
	per_win_2DOF_lls_scaling_factor = PeakSeq_slope(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);
	per_win_1DOF_lls_scaling_factor = slope(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);

	//double pears_corr = pearson_correlation(prof1_win_total_sigs, prof2_win_total_sigs, n_processed_wins);
	delete [] prof1_win_total_sigs;
	delete [] prof2_win_total_sigs;
}

//double PeakSeq_slope(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt, int n_processed_wins)
double PeakSeq_slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
    // Compute the slope from 2 DOF (slope) LLS regression.
    double num = 0.0f;
    double den = 0.0f;

    double sx = 0.0f;
    double sy = 0.0f;
    double sxx = 0.0f;
    double sxy = 0.0f;

    for(int i = 0; i < n_processed_wins; i++)
    {
		sx += prof1_total_sigs_per_win[i];
		sy += prof2_total_sigs_per_win[i];
		sxx += prof1_total_sigs_per_win[i] * prof1_total_sigs_per_win[i];
		sxy += prof1_total_sigs_per_win[i] * prof2_total_sigs_per_win[i];

		//sx += control_frag_cnt->at(i);
		//sy += chip_seq_frag_cnt->at(i);
		//sxx += control_frag_cnt->at(i) * control_frag_cnt->at(i);
		//sxy += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
    } // i loop

    num = n_processed_wins * sxy - sx * sy;
    den = n_processed_wins * sxx - sx * sx;

	return(num / den);
}

//double slope(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt)
double slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
	// Compute the slope from 1 DOF (slope) LLS regression.
	double num = 0.0f;
	double den = 0.0f;
	for(int i = 0; i < n_processed_wins; i++)
	{
		num += prof1_total_sigs_per_win[i] * prof2_total_sigs_per_win[i];
		den += prof1_total_sigs_per_win[i] * prof1_total_sigs_per_win[i];
		//num += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
		//den += control_frag_cnt->at(i) * control_frag_cnt->at(i);
	} // i loop

	return(num / den);
}


// Copied from http://www.vias.org/tmdatanaleng/cc_corr_coeff.html
// Also check http://www.statsdirect.com/help/regression_and_correlation/sreg.htm
//double pearson_correlation(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt)
double pearson_correlation(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins)
{
	double num = 0.0f;
	//double den = 0.0f;

	double prof1_mean = 0.0f;
	double prof2_mean = 0.0f;
	for(int i_win = 0; i_win < n_processed_wins; i_win++)
	{
		//control_mean += control_frag_cnt->at(i_frag);
		prof1_mean += prof1_total_sigs_per_win[i_win];
		prof2_mean += prof2_total_sigs_per_win[i_win];
	} // i_win loop.
	prof1_mean /= n_processed_wins;
	prof2_mean /= n_processed_wins;

    //for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
	for(int i_win = 0; i_win < n_processed_wins; i_win++)
    {
		num += (prof1_total_sigs_per_win[i_win] - prof1_mean) * (prof2_total_sigs_per_win[i_win] - prof2_mean);
	} // i_win loop.
	
	double prof1_var = 0.0f;
    for(int i_win = 0; i_win < n_processed_wins; i_win++)
    {
		prof1_var += (prof1_total_sigs_per_win[i_win] - prof1_mean) * (prof1_total_sigs_per_win[i_win] - prof1_mean);
	}	
	prof1_var = pow(prof1_var, .5);

	double prof2_var = 0.0f;
    for(int i_win = 0; i_win < n_processed_wins; i_win++)
    {
		prof2_var += (prof2_total_sigs_per_win[i_win] - prof2_mean) * (prof2_total_sigs_per_win[i_win] - prof2_mean);
	}
	prof2_var = pow(prof2_var, .5);

	double r = num / (prof1_var * prof2_var);

	return(r);
}



