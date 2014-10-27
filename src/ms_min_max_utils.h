#ifndef __MIN_MAX_UTILS__
#define __MIN_MAX_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;

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

	// These are the nodes that pro/preceed this node in scale space.
	t_extrema_node* lower_node;
	t_extrema_node* higher_node;

	// Height at the extrema value.
	double height_at_extrema;

	// This is the flag that identifies whether this node has a multipath issue.
	bool node_ok;
};

void delete_extrema_nodes(vector<t_extrema_node*>* extrema_nodes);

// Given a set of peaks, assign id's to the lowest levels peaks, based on the set of maxima.
void set_matching_maxima_per_minima_regions(vector<t_annot_region*>* min_regions, vector<t_annot_region*>* max_regions);

bool sort_extremas_per_posn(t_extrema_node* node1, t_extrema_node* node2);
bool sort_extremas_per_height(t_extrema_node* node1, t_extrema_node* node2);
bool sort_extremas_per_decreasing_height(t_extrema_node* node1, t_extrema_node* node2);

//double* load_data(char* data_fp, int l_signal);
void trace_local_extrema_per_decomp(vector<double*>* decomps_per_scale,
									int l_decomp_signal,
									vector<vector<t_extrema_node*>*>* maxes_per_scale, 
									vector<vector<t_extrema_node*>*>* mins_per_scale);

void get_matching_extrema_regions_per_extrema_nodes(vector<t_extrema_node*>* minima_nodes, vector<t_extrema_node*>* maxima_nodes,
													vector<t_annot_region*>* matching_mins, vector<t_annot_region*>* matching_maxes);

int extrema_posn_accessor(void* extrema_ptr);

void check_multitraced_extrema(int i_scale, 
	vector<t_extrema_node*>* cur_scale_extrema, vector<t_extrema_node*>* higher_scale_extrema, 
	bool* traced_path_indicator, 
	int l_decomp,
	bool& multipath_extrema_exists, bool& no_path_extrema_exists);

void prefilter_data_per_peak_side(double* signal, int l_signal, int l_smallest_side);

void prune_extrema_derivative_map(vector<t_extrema_node*>* mins,
	vector<t_extrema_node*>* maxes,
	signed int* derivative_map, 
	int l_smallest_side,
	int l_signal, 
	vector<t_extrema_node*>* pruned_mins,
	vector<t_extrema_node*>* pruned_maxes,
	signed int* pruned_derivative_map);

void trace_local_minima(vector<t_extrema_node*>* higher_scale_minima, 
	vector<t_extrema_node*>* cur_scale_minima, 	
	bool* cur_scale_min_indicator,
	bool* cur_scale_traced_min_indicator,
	int* cur_scale_derivative_map,	
	int l_signal,
	int i_scale);

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale);

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale,
	double zero_deriv);

#endif // __MIN_MAX_UTILS__
