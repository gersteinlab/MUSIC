#ifndef __GSL_FFT_FILTER_UTILS__
#define __GSL_FFT_FILTER_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;

int sort_doubles_descending(const void* p1, const void* p2);

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
	const char* op_file_prefix);

vector<double*>* multiscale_avg_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix);

double* mean_filter_data(double* track_data,
						int l_track_data, 
						int l_averaging_win,
						int skip_value); // This is the value to be skipped from the window values.

double* median_filter_data(double* track_data,
	int l_track_data, 
	int l_averaging_win,
	int skip_value); // This is the value to be skipped from the window values.

double* mapability_aware_median_filter(double* signal_profile, int l_profile,
	double* scaled_mapability_profile, int l_mapability_profile,
	double max_mapable_signal_2_use_in_filter,
	int l_mapability_filtering_win);

void get_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale);

void get_mapability_aware_median_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale,
	double* normalized_mapability_signal,
	int l_mapability_profile,
	double max_mapability_signal);

bool get_origin_of_symmetry(double* filtered_data, int l_signal, int search_start_i, int& origin, int l_mirror_check);

double* sum_tracks(vector<double*>* multitrack_data, int l_signal);

void get_next_2_exp(int val, int& larger_exp_val, int& expon);

#endif // __GSL_FFT_FILTER_UTILS__