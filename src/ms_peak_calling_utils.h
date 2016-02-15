#ifndef __PEAK_CALLING_UTILS__
#define __PEAK_CALLING_UTILS__

struct t_annot_region;

struct t_ER_info
{
	double rev_strand_mass_of_center;
	double fore_strand_mass_of_center;

	int fore_climax_posn;
	int rev_climax_posn;

	double fore_max_val;
	double rev_max_val;

	double total_fore_mass;
	double total_rev_mass;

	double total_chip_mass;
	double total_control_mass;

	int climax_posn;

	double max_chip_mass;
	double max_control_mass;

	int mappable_trough;
};

void write_decomposition_WAV(char* chip_reads_dir,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_mapability_filtering_win,
		double max_normalized_mapability_signal);

void write_logR_profiles(char* chip_reads_dir,
							char* control_reads_dir, 
							int l_fragment,
							char* mapability_signal_dir,
							int l_read_mapability_signal,
							double base_scale_l_win,
							double end_scale_l_win,
							double log_step,
							int l_mapability_filtering_win,
							double max_normalized_mapability_signal);

void write_decomposition_bedGraphs(char* chip_reads_dir,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_mapability_filtering_win,
		double max_normalized_mapability_signal);

double* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile);

int get_trough_posn_per_ER(double* signal_profile, int l_profile, 
							double* multi_mapp_signal, int l_multi_map_signal, double max_multi_mapp_val,
							int ER_start, int ER_end, 
							int l_trough_win);

void dump_ERs_per_q_value_and_summit(vector<t_annot_region*>* regions, char* op_fp);
void dump_broadPeak_formatted_op(vector<t_annot_region*>* regions, char* op_fp);

vector<t_annot_region*>* filter_strand_uneven_peaks(vector<t_annot_region*>* peak_regions,
		char* chrom,
		char* chip_reads_dir,
		int l_fragment,
		double min_forward_reverse_strand_signal_fraction);

vector<t_annot_region*>* set_per_strand_info_per_peaks(vector<t_annot_region*>* peak_regions,
		char* chrom,
		char* chip_reads_dir,
		int l_fragment);

double* load_normalized_multimappability_profile(char* mapability_signal_profile_fp, int& l_mapability_profile);

void dump_ERs_per_q_value_and_summit(vector<t_annot_region*>* regions, char* op_fp);

enum {P_VAL_PER_WIN_MEAN_SIGNAL, P_VAL_PER_WIN_MAX_SIGNAL};

void get_peaks(char* chip_reads_dir,
		char* control_reads_dir,
		vector<char*>* chr_ids_2_use,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_bin,
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
		double p_val_min_pruning_flip_prob,
		bool find_mappable_troughs,
		double trimmed_SSER_p_value_cutoff,
		double q_value_cutoff);

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