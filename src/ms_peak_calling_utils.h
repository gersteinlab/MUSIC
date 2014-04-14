#ifndef __PEAK_CALLING_UTILS__
#define __PEAK_CALLING_UTILS__

struct t_annot_region;

vector<t_annot_region*>* filter_strand_uneven_peaks(vector<t_annot_region*>* peak_regions,
		char* chrom,
		char* chip_reads_dir,
		int l_fragment,
		double min_forward_reverse_strand_signal_fraction);

enum {P_VAL_PER_WIN_MEAN_SIGNAL, P_VAL_PER_WIN_MAX_SIGNAL};

void get_peaks(char* chip_reads_dir,
		char* control_reads_dir,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_p_val_norm_win,
		int l_mapability_filtering_win,
		double max_normalized_mapability_signal,
		double filtered_vs_non_filtered_max_scaler,
		double signal_profile_scaling_factor,
		double rd_end_pruning_p_val,
		double min_per_strand_evennes_fraction,
		int p_value_computation_type,
		bool do_replace_profiles_w_smoothed_profiles,
		bool do_BJ_correction_on_minima,
		bool do_filter_minima_per_length_min_base_scale_l_win,
		bool do_post_peak_p_value_minimization,
		double p_val_min_pruning_flip_prob);

void get_peaks_no_control(char* chip_reads_dir,
	int l_fragment,
	char* mapability_signal_dir,
	int l_read_mapability_signal,
	double base_scale_l_win,
	double end_scale_l_win, 
	double log_step,
	int l_p_val_norm_win,
	int l_mapability_filtering_win,
	double max_normalized_mapability_signal,
	double filtered_vs_non_filtered_max_scaler);

#endif // __PEAK_CALLING_UTILS__