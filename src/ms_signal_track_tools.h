#ifndef __SIGNAL_TRACK_FILE_INTERFACE__
#define __SIGNAL_TRACK_FILE_INTERFACE__

#include <vector>

using namespace std;

struct t_annot_region;
class t_rng;

double* load_per_nucleotide_binary_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_binary_profile_per_bedgraph(char* bgr_fp, bool dump_binary, char* op_dir);
void dump_per_nucleotide_binary_profile(double* profile, int l_profile, char* op_fp);
void dump_bedGraph_per_per_nucleotide_binary_profile(double* signal_profile_buffer, int l_profile, char* chrom, char* op_fp);

double* get_zero_indexed_per_one_indexed_data(double* one_indexed_data, int l_profile);
double* get_one_indexed_per_zero_indexed_data(double* zero_indexed_data, int l_profile);

double* quantize_per_nucleotide_profiles(double* profile, int l_profile, vector<double>* thresholds, vector<double>* quantized_vals);

double** joint_prune_multiple_profiles_per_thresholds(double** profiles, int n_profiles, int l_profile, double* threshold_per_profile, int& l_pruned_profile);

unsigned char* load_per_nucleotide_binary_uchar_profile(char* binary_per_nucleotide_profile_fp, int& l_profile);
void dump_per_nucleotide_uchar_binary_profile(unsigned char* signal_profile, int l_profile, char* op_fp);

double* get_block_permute_profile(double* profile_data, t_rng* rng, int l_profile, int& l_permuted, int l_block);
double* get_random_circular_shifted_profile(double* profile_data, t_rng* rng, int l_profile);

void get_profile_extrema(double* profile, int l_profile, double& prof_min, double& prof_max);

double* extract_one_indexed_profile_per_profile(double* signal_profile_buffer, int l_profile, int start, int end, int& l_extracted_profile);
double* copy_profile(double* signal_profile, int l_profile);

vector<t_annot_region*>** get_peaks_per_per_nucleotide_signal_profile(double* signal_profile, char* chrom, int l_data, double min_thresh, double max_thresh);
vector<t_annot_region*>** get_valleys_per_per_nucleotide_signal_profile(double* signal_profile, char* chrom, int l_data, double min_thresh, double max_thresh);

void exclude_regions_from_signal_profiles(double* signal_profile, int l_profile, vector<t_annot_region*>* regions_2_exclude, double* pruned_signal_profile, int& l_pruned_profile);

void floorize_profile(double* signal_profile, int l_profile);
void get_log_plus_one_profile(double* signal_profile, double base, int l_profile);
double* get_zero_profile(int l_profile);

#endif // __SIGNAL_TRACK_FILE_INTERFACE__

