#ifndef __GSL_FFT_FILTER_UTILS__
#define __GSL_FFT_FILTER_UTILS__

#include <vector>
using namespace std;

struct t_plateau
{
	int start;
	int end;
};

enum{EXTREMA_MAX, EXTREMA_MIN};
struct t_extrema_node
{
	int extrema_posn;
	int extrema_type;

	// This is the scale at which this extrema is. Each extrema is found at a scale. 
	int scale;

	// The peak to which this extrema belongs. (For maxima only).
	//t_ms_peak* peak;

	// These are the nodes that pro/preceed this node in scale space.
	t_extrema_node* lower_node;
	t_extrema_node* higher_node;

	// Height at the extrema value.
	double height_at_extrema;

	// This is the flag that identifies whether this node has a multipath issue.
	bool node_ok;
};

void delete_extrema_nodes(vector<t_extrema_node*>* extrema_nodes);

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale,
	double zero_deriv);

struct t_annot_region;

int sort_doubles_descending(const void* p1, const void* p2);

void get_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale);

double* median_filter_data(double* track_data,
	int l_track_data, 
	int l_averaging_win,
	int skip_value);

double* mapability_aware_median_filter(double* signal_profile, int l_profile,
	double* scaled_mapability_profile, int l_mapability_profile,
	double max_mapable_signal_2_use_in_filter,
	int l_mapability_filtering_win);

vector<double*>* multiscale_conv_filter_data(double* cur_real_track_data, int i_t, int l_track_data, double scale_start, double scale_end, double log_scale_step, vector<double>* scales_per_i_scale);
double* filter_data_per_odd_length_symmetric_filter(double* cur_real_track_data, int l_track_data, double* cur_filter_array, int l_filter, int n_iters);

double* sum_tracks(vector<double*>* multitrack_data, int l_signal);

void get_next_2_exp(int val, int& larger_exp_val, int& expon);

#endif // __GSL_FFT_FILTER_UTILS__