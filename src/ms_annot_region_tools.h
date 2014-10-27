#ifndef __ANNOT_REGION_TOOLS__
#define  __ANNOT_REGION_TOOLS__

#include <vector>

using namespace std;

// Stores the sorting information for the regions.
struct t_sorting_info
{
	// Following are the cumulative ends for the region: Basically, for a current read, sets the smallest and largest start and end at this read position.
	// This becomes useful when doing binary searches over the regions that overlap.
	int cumulative_sorted_start;
	int cumulative_sorted_end;	
};

struct t_significance_info
{
	double log_p_val;
	double log_q_val;
};

// This is just a region on a chromosome.
struct t_annot_region
{
	char* chrom;
	char strand;
	int start; // If strand is '-', this coordinate is with respect to reverse complement on the forward strand.
	unsigned int score;
	char* name;
	int end;

	// Load the remaining information if it is available.
	int n_exons;
	int thick_start;
	int thick_end;
	vector<t_annot_region*>* intervals; // These are the exons in the region. Some bed files include this information.

	t_sorting_info* sort_info;
	t_significance_info* significance_info;
	void* data; // Extra data about the region, this is a pointer to a data structure that holds specific data about the region depending on the loaded data type: narrowPeak, GFF, ... Can also be used for other purposes for storing data about a region.
};

struct t_pooled_peak_association_info
{
	vector<t_annot_region*>* agg_peak_oling_peaks;
	int* oling_peaks_counts;
};

vector<t_annot_region*>* get_top_region_per_total_coverage(vector<t_annot_region*>* regions, double max_covg);

vector<t_annot_region*>* extract_region_mids(vector<t_annot_region*>* regions, int l_mid_reg);

vector<t_annot_region*>* extract_region_ends_by_fraction_extension(vector<t_annot_region*>* regions, double ext_fraction,
	bool extract_5p, 
	bool extract_3p);

vector<t_annot_region*>* extract_region_ends_by_l_extension(vector<t_annot_region*>* regions, int l_ext,
	bool extract_5p, 
	bool extract_3p);

vector<vector<t_annot_region*>*>** get_shuffled_merged_regions_per_profile(vector<vector<t_annot_region*>*>* regions_per_profile, int n_shufflings);

double** get_pw_association_counts_per_regions(vector<char*>* ids_per_profile, 
	vector<t_annot_region*>* aggregated_signal_peaks,
	vector<vector<t_annot_region*>*>* per_profile_peak_regions);

struct t_intersect_info
{
	t_annot_region* src_reg;
	t_annot_region* dest_reg;
	int l_overlap;
};

void sort_set_sorting_info(vector<t_annot_region*>* regions, bool (sort_regions_callback)(t_annot_region*, t_annot_region*));

/* 
This is a re-structured representation of a list of regions which holds three different 
lists that do not re-allocate but divides the regions, these lists can be used directly in operations like overlapping,
merging, ... because they are directly comparable.
*/
struct t_sorted_annot_region_lists
{
	vector<t_annot_region*>** regions_per_chrom;
	vector<t_annot_region*>** pos_strand_regions_per_chrom;
	vector<t_annot_region*>** neg_strand_regions_per_chrom;
	vector<char*>* chr_ids; // These are the ids that have smae ordering in the region lists.
};

// Following is the extra information that is stored in a narrowPeak file. This information is stored for each narrowPeak entry.
#define NARROWPEAK_INFO (0x1234)
struct t_narrowPeak_info
{
	int info_type;
	int score;
	char strand;
	double signal_value;
	double p_value;
	double q_value;
	int peak_pos; // This is usually important. This is where the peak is.
};

void get_benjamini_hochberg_corrected_p_values(vector<t_annot_region*>* regions);

void dump_GFF_per_region_intervals(vector<t_annot_region*>* regions, char* gff_fp, char* source_name, char* interval_parent_feature_str, char* interval_child_feature_str);

vector<t_annot_region*>* merge_regions_per_super_regions(vector<t_annot_region*>* regions, vector<t_annot_region*>* super_regions);

t_sorted_annot_region_lists* restructure_annot_regions(vector<t_annot_region*>* regions); // Re-structure the list using the chromosome ids retreived from itself.
t_sorted_annot_region_lists* restructure_annot_regions(vector<t_annot_region*>* regions, vector<char*>* chr_ids); // Re-structure the list for a given list of chromosome ids.
void delete_restructured_annot_regions(t_sorted_annot_region_lists* restructured_region_lists);

vector<t_annot_region*>* load_BED_w_p_values(char* regions_bed_fp);
void dump_BED_w_q_values(vector<t_annot_region*>* regions, char* op_fp);
void dump_BED_w_p_values(vector<t_annot_region*>* regions, char* op_fp);
vector<t_annot_region*>* load_BED(char* bed_fp);
vector<t_annot_region*>* load_BED12(char* bed_fp);
vector<t_annot_region*>* load_BED_with_line_information(char* bed_fp);
void dump_BED(char* bed_fp, vector<t_annot_region*>* annot_regions);
vector<t_annot_region*>* load_Interval(char* interval_fp);
vector<t_annot_region*>* get_all_intervals_per_composite_intervals(vector<t_annot_region*>* composite_intervals);
void dump_Interval(char* interval_fp, vector<t_annot_region*>* annot_regions);
vector<t_annot_region*>* load_narrowPeak(char* narrowPeak_fp); // This function also loads the broadPeak files.
void extend_BED(vector<t_annot_region*>* regions, int l_extend, bool extend_3p, bool extend_5p);

bool validate_region_coords(vector<t_annot_region*>* regions);

//t_annot_region* copy_region(t_annot_region* region);

void load_chromosome_lengths_per_tabbed_file(char* chr_lengths_fp, vector<char*>* chr_ids, vector<int>* chr_lengths);

t_annot_region* get_empty_region();

// Divide into strands first, then divide into chromosomes, then sort.
vector<t_annot_region*>* sort_regions_per_chromosome_per_strand(vector<t_annot_region*>* annot_region_list);

void delete_annot_regions(vector<t_annot_region*>* region_list);
void delete_annot_regions(t_annot_region* region);
void delete_annot_regions_with_line_information(vector<t_annot_region*>* region_list);

// Get the list of interregion distances for a list of regions. 
vector<int>* inter_region_distances(vector<t_annot_region*>* annot_regions);

bool check_line_skip(char* cur_line);

vector<char*>* get_chr_ids(vector<t_annot_region*>* annot_regions);

vector<t_annot_region*>* get_regions_per_chromosome(vector<t_annot_region*>* annot_regions,
													char* chr_id);

vector<t_annot_region*>* get_regions_per_strand(vector<t_annot_region*>* annot_regions,
													char strand);

vector<t_annot_region*>* get_regions_per_min_run(vector<t_annot_region*>* annot_regions, int l_min_run);
vector<t_annot_region*>* get_regions_per_max_run(vector<t_annot_region*>* annot_regions, int l_max_run);

vector<t_annot_region*>* get_regions_per_end_max(vector<t_annot_region*>* annot_regions, int end_max);

// Get left/right flanking regions (with respect to forward strand).
vector<t_annot_region*>* get_left_flanking_regions(vector<t_annot_region*>* regions, int l_flank);
vector<t_annot_region*>* get_right_flanking_regions(vector<t_annot_region*>* regions, int l_flank);

double coverage(vector<t_annot_region*>* annot_regions);

/*
Region intersection
*/
vector<t_annot_region*>* intersect_annot_regions(vector<t_annot_region*>* annot_regions1,
														vector<t_annot_region*>* annot_regions2,
														bool match_strands,
														bool find_all_overlaps);

vector<t_annot_region*>* intersect_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2,
													bool find_all_overlaps);

void delete_intersect_info(vector<t_annot_region*>* regions);

/*
Region merging
*/
vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* annot_regions1,
												vector<t_annot_region*>* annot_regions2,
												int max_gap, 
												bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<vector<t_annot_region*>*>* annot_regions_list,
												int max_gap, 
												bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* total_annot_regions,
											int max_gap, 
											bool match_strands);

vector<t_annot_region*>* merge_annot_regions(vector<t_annot_region*>* total_annot_regions,
											int max_gap);

/*
Region exclusion.
*/
vector<t_annot_region*>* exclude_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2, 
													bool strand_specific);

// Backend functions.
vector<t_annot_region*>* exclude_annot_regions(vector<t_annot_region*>* annot_regions1,
													vector<t_annot_region*>* annot_regions2);

/*
Region subtraction.
*/
vector<t_annot_region*>* subtract_annot_regions(t_annot_region* region1, t_annot_region* region2, bool strand_specific);

// Dump the distances between the regions. This gives an overall statistic about how dispersed the data is.
vector<int>* inter_region_distances(vector<t_annot_region*>* annot_regions);

void dump_inter_region_distances(vector<t_annot_region*>* annot_regions1,
									vector<t_annot_region*>* annot_regions2);


//// Load the regions, find the transcripts with the same gene id's.
//vector<t_annot_region*>* get_gene_regions_per_transcripts(vector<t_annot_region*>* mixed_regions, 
//	const char* gene_feature_id, 
//	char* (*get_gene_id_per_group_str(char*)));
//
//// Associate the exons for the 
//vector<t_annot_region*>* get_transcript_regions_per_exons(vector<t_annot_region*>* mixed_regions, 
//	const char* transcript_feature_id, 
//	char* (*get_transcript_id_per_group_str(char*)));

// Necessary for fast comparison of the regions.
bool sort_regions(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_increasing_length(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_decreasing_length(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_ends(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_score(t_annot_region* region1, t_annot_region* region2);
bool sort_regions_per_increasing_p_value(t_annot_region* reg1, t_annot_region* reg2);
bool sort_regions_per_increasing_q_value(t_annot_region* reg1, t_annot_region* reg2);

// Duplicate a region.
t_annot_region* duplicate_region(t_annot_region* region_2_dup);

//vector<t_annot_region*>* intersect_regions_binary_search(vector<t_annot_region*>* src_regions, vector<t_annot_region*>* dest_regions, bool match_strands);

vector<int>* get_inter_region_distance_distribution(vector<t_annot_region*>* annot_regions1, 
	vector<t_annot_region*>* annot_regions2);

vector<t_annot_region*>* subtract_annot_regions(vector<t_annot_region*>* regions1, vector<t_annot_region*>* regions2);

int get_reg2reg_distance(t_annot_region* reg1, t_annot_region* reg2);

// Locate the index of the region (out of a list of regions) whose start is just to the left of start_posn.
int locate_posn_per_sorted_obj_list(void* prop_ptr, vector<void*>* obj_list, int i, int j, int (*comparor)(void*, void*));
int locate_posn_per_sorted_obj_list(int posn, vector<void*>* obj_list, int i, int j, int (*prop_accession)(void*));
int locate_posn_region_per_region_starts(int start_posn, vector<t_annot_region*>* region_list, int i, int j);
int locate_posn_per_sorted_posn_list(int posn, vector<int>* sorted_posn_list, int i, int j);

void get_overlap_statistics(vector<t_annot_region*>* regs1, vector<t_annot_region*>* regs2, 
							vector<double>* bases, vector<int>* stat_counts, double sim_step);

int region_3p_accessor(void* void_ptr);
int region_5p_accessor(void* void_ptr);

#endif // __ANNOT_REGION_TOOLS__
