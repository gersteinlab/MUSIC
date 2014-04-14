#ifndef __SIGNAL_ENRICHMENT_UTILS__
#define __SIGNAL_ENRICHMENT_UTILS__

#include <vector>
using namespace std;
struct t_annot_region;

#define  __POISSON_BCKGRND__
#undef __GEOMETRIC_BCKGRND__

double* buffer_log_factorials(int n);

vector<t_annot_region*>* prune_region_ends_per_lost_peak_mass(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control);
vector<t_annot_region*>* prune_region_ends_per_signal_distribution(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control);
vector<t_annot_region*>* prune_region_ends_per_window_average(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control,
	double fdr, int l_win);
vector<t_annot_region*>* prune_region_ends_per_p_value_minimization(vector<t_annot_region*>* regions, double* signal_profile, int l_profile, double* control_profile, int l_control,
	int l_fragment);
vector<t_annot_region*>* prune_region_ends_per_p_value_minimization_no_control(vector<t_annot_region*>* regions, double* signal_track_data, int l_profile, 
	double scaling_factor,
	double per_win_bckgrnd_win_val);

vector<t_annot_region*>* prune_region_ends_per_modified_binomial_p_value_minimization(vector<t_annot_region*>* regions, 
	double flip_prob,
	double* signal_profile, int l_profile, 
	double* control_profile, int l_control,
	int l_fragment);

double get_modified_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, 
	double flip_prob, 
	int region_start, int region_end, 
	double* _log_factorials, 
	double scaling_factor,
	bool length_normalization, int l_norm_win);

vector<t_annot_region*>* prune_region_ends_per_p_value_minimization_per_window(vector<t_annot_region*>* regions, 
	double* signal_profile, int l_profile, 
	double* control_profile, int l_control,
	int l_fragment);

vector<t_annot_region*>* prune_region_ends_per_adaptive_thresholding(vector<t_annot_region*>* regions, double* signal_track_data, int l_profile_signal, 
	double* control_track_data, int l_control_signal,
	int l_fragment);

//double get_binomial_p_val_per_enriched_background_values(int enriched_value, int  background_value);
double get_binomial_p_val_per_enriched_background_values(int enriched_value, int background_value, double* _log_factorials);

double compare_read_numbers_per_peaks(double p1_profile_signal, double p1_control_signal, double p2_profile_signal, double p2_control_signal, double normalization_factor);

double get_fishers_exact_p_val(double n11, double n12, double n21, double n22);

double get_per_win_percentile_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, 
	int region_start, int region_end, 
	double* _log_factorials, double scaling_factor,
	int l_norm_win,
	double win_rank_per_region);

double* get_per_win_profile_per_nuc_profile(double* profile, int l_profile, int l_win);

double* get_per_win_median_profile_per_background_window_profile(double* per_bckgrnd_win_signal,
	int n_bckgrnd_wins,
	int l_background_est_win,
	int l_background_norm_win);

double get_binomial_pvalue_per_region_no_control(double* signal_profile, int l_profile, 
	int region_start, int region_end, 
	double per_win_bckgrnd_win_val,
	double* _log_factorials, 
	double scaling_factor,
	int l_norm_win,
	bool normalize);

double get_fdr_based_bckgrnd_window_value(double* signal_profile, int l_profile,
												double target_fdr,
												int thresh_win_length,
												int l_background_est_win,
												int l_pval_norm_win,
												double significant_log_p_val_cutoff,
												double scaling_factor);

double get_binomial_pvalue_per_region_neighborhood_window_control(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win, int l_signal, int l_control, int n_vicinity_windows);

double get_binomial_pvalue_per_region_neighborhood_window_control_per_max_profile_signal(double* signal_track_data, double* control_track_data, int region_start, int region_end, double* _log_factorials, double scaling_factor,
	bool length_normalization, int l_norm_win, int l_signal, int l_control, int n_vicinity_wins);

double get_binomial_pvalue_per_region(double* signal_track_data, double* control_track_data, 
	int region_start, int region_end, 
	double* _log_factorials, 
	double scaling_factor,
	bool length_normalization, int l_norm_win);

double get_geometric_thresholds_per_win(double* signal_profile,
								int l_profile,
								int l_win,
								int win_start,
								double target_fdr,
								int min_thresh,
								int max_thresh);

double get_poisson_thresholds_per_win(double* signal_profile,
								int l_profile,
								int l_win,
								int win_start,
								double target_fdr,
								int min_thresh,
								int max_thresh);


#endif // __SIGNAL_ENRICHMENT_UTILS__