#ifndef __PROFILE_NORMALIZATION__
#define __PROFILE_NORMALIZATION__

void get_per_window_scaling_factors_for_profile1_per_profile2(double* signal_profile1, int l_prof1, double* signal_profile2, int l_prof2, int l_win,
	double& per_win_2DOF_lls_scaling_factor,
	double& per_win_1DOF_lls_scaling_factor,
	double& total_sig_scaling_factor);

double get_per_window_scaling_factor_per_track2(double* signal_profile1, double* signal_profile2, int l_prof1, int l_prof2);

double PeakSeq_slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);
double slope(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);
double pearson_correlation(double* prof1_total_sigs_per_win, double* prof2_total_sigs_per_win, int n_processed_wins);

#endif // __PROFILE_NORMALIZATION__