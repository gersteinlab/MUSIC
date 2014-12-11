#ifndef __MAPPED_READ_FILE_INTERFACE__
#define __MAPPED_READ_FILE_INTERFACE__

#include <vector>
using namespace std;

//struct t_profile_site;
struct t_mapped_fragment;
struct t_frag_end;

const int MEG_BASE = 1000 * 1000;
const int K_BASE = 1000;

class t_string;
struct t_annot_region;


// MApping information for a mapped read.
struct t_mapping_info
{
	char* chrom;
	int base_posn;
	int sequenced_length;
	char strand;
	char* mapping_quality_str; // CIGAR string.
};

// This is a sequenced tag w/o any mapping information: Does not have any information about mapping, initially. In case there is
// mapping information, it is stored under the member named mapping_info.
struct t_sequenced_read
{
	char* id;
	t_mapping_info* mapping_info;
	char* quality_str;
	char* nucs;
};

// This is a mapped fragment: It can be a part of a read.
struct t_mapped_fragment
{
	int base_index;
	int sequenced_fragment_length; // Necessary when building enrichment profile. This is the fragment length coming from ChIP-Seq data.
	char strand_char;
};

struct t_mapped_read
{
	char* mapping_str;

	// The start of the read: 5' position.
	int base_index;

	// Left side posn.
	//int left_base_posn;

	// This is the total span of the read.
	int span; 
	char strand;
};

#define MAX_N_PAIRS (10)

bool check_genome_index_update_per_CIGAR_entry(char entry_char);
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char);

// Load sequenced reads directly from the fastq file.
// Loads and pools the reads.
void load_sequenced_reads_per_fastq(char* fastq_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void dump_phred_quality_distribution(vector<t_sequenced_read*>* sequenced_reads, char* op_fp);

// Subsample the reads to a desired size.
void subsample_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads, 
	vector<t_sequenced_read*>* subsampled_sequenced_reads, 
	int n_subsample_size);

void subsample_mapped_reads(vector<t_mapped_read*>* mapped_reads, 
	vector<t_mapped_read*>* subsampled_mapped_reads, 
	int n_subsample_size);

void dump_fastq(vector<t_sequenced_read*>* sequenced_reads, char* op_fastq_fp);
void delete_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads);

// Load the mapped read files and preprocess them, from which the sequenced fragments can be loaded.
void load_mapped_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_ELAND(char* eland_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_tagAlign(char* tagalign_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_bowtie(char* bowtie_fp, vector<t_sequenced_read*>* sequenced_reads);

void count_mapped_reads_per_file(char* mrf_fp, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	double n_total_reads);

// Following are the preprocessing functions for preprocessing the mapped reads file 
void preprocess_mapped_reads_file(char* mrf_fp, char* parsed_reads_op_dir, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_ids);
void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);
void preprocess_BED5_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str);

int get_l_signal_per_reads(char* reads_fp, int l_ext_tag);

//void buffer_per_nucleotide_profile_no_buffer(char* sorted_read_fp, int l_extended_tag, double* signal_profile_buffer, int l_buffer, int& l_data);
void buffer_per_nucleotide_profile_no_buffer(char* sorted_read_fp, const int l_extended_tag, 
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int l_buffer, int& l_data);

void buffer_per_nucleotide_preprocessed_read_profile_no_buffer(char* sorted_read_fp,
	double* signal_profile_buffer, double* forward_strand_signal, double* reverse_strand_signal, 
	int max_l_read, // This is the length of the longest spliced read to be processed.
	int l_buffer, int& l_data);

double get_n_mapped_nucs(vector<t_mapped_fragment*>* fragments);

struct t_read_line_sorting_info
{
	int start;
	char* read_line;
};

bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2);

vector<char*>* sort_bucket_read_lines(char* bucket_fp);














void get_per_strand_read_stats(vector<t_annot_region*>* regions,
	char* preprocessed_reads_dir);




// BINARY INTERFACE IS OBSOLETE
void load_mapped_reads_file(char* read_file_path, vector<char*>* chr_ids, 
														vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
														void (preprocess_mapped_read_line)(char* cur_line, 
														char* chrom, 
														int& chr_index, int& sequenced_length, 
														char& strand_char, 
														char* mapping_quality_str));

// Load reads reads/fragments.
void load_fragments_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
	vector<vector<t_mapped_fragment*>*>* fore_strand_fragments_per_chr, vector<vector<t_mapped_fragment*>*>* rev_strand_fragments_per_chr, 
	int max_n_pcr_amplified_reads);

void load_fragments(char* mapped_reads_fp, 
	vector<t_mapped_fragment*>* fore_strand_fragments, vector<t_mapped_fragment*>* rev_strand_frag, 
	int max_n_pcr_amplified_reads = 10000);

void load_reads_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
	vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
	int max_n_pcr_amplified_reads);

void load_reads(char* mapped_reads_fp, 
	vector<t_mapped_read*>* fore_strand_reads, vector<t_mapped_read*>* rev_strand_reads, 
	int max_n_pcr_amplified_reads);

void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments);
void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments);

// Following are for doing binary searches over the fragments.
bool sort_mapped_reads_per_3p(t_mapped_read* read1, t_mapped_read* read2);
bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2);
int read_5p_accessor(void* obj_ptr);
//int read_3p_accessor(void* obj_ptr);

void delete_mapped_reads(vector<t_mapped_read*>* mapped_read);
void delete_mapped_read(t_mapped_read* mapped_read);

int fragment_5p_accessor(void* obj_ptr);
int fragment_3p_accessor(void* obj_ptr);

void get_read_statistics_per_region(vector<t_mapped_fragment*>* fragments, 
									vector<int>* all_fragment_3p_posns, 
									t_annot_region* region, 
									double& n_mapped_reads, 
									double& n_mapped_nucs);

bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced);
void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										int& l_cur_entry,
										char& entry_type_char);

void delete_fragments(vector<t_mapped_fragment*>* fragment_list);
void delete_fragments(t_mapped_fragment** fragment_list);

bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str);
bool sort_mapped_fragments(t_mapped_fragment* frag1, t_mapped_fragment* frag2);
bool sort_mapped_fragments_per_3p(t_mapped_fragment* frag1, t_mapped_fragment* frag2);

void preprocessed_read_file_iterator(char* mapped_reads_fp, 
	void (per_read_callback)(char*, char, int, void*), 
	void (per_fragment_callback)(char*, char, int, void*),
	void* per_read_callback_param,
	void* per_fragment_callback_param);

// Prune fragments/reads
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
				vector<t_mapped_read*>* pruned_forward_reads, 
				vector<t_mapped_read*>* pruned_reverse_reads);
//vector<t_mapped_fragment*>* prune_fragment_list_by_strand(vector<t_mapped_fragment*>* chr_fragments, char strand_2_prune, int n_max_reps);

/*
Following function is the only function that does tag extension with enrichment_mapped_fragment_length parameter. This should not be handled anywhere else.
*/
vector<t_mapped_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_mapped_fragment*>* fore_frag_list, vector<t_mapped_fragment*>* rev_frag_list, int enrichment_mapped_fragment_length);

int* get_n_reads_per_window(int n_wins, vector<t_mapped_fragment*>* frag_list);

void exclude_reads_per_regions(vector<char*>* read_chr_ids, 
	vector<vector<t_mapped_read*>*>* reads_per_chrs, 
	vector<vector<t_mapped_read*>*>* no_overlap_reads_per_chrs,
	vector<t_annot_region*>* regions_2_exclude);

void exclude_reads_per_regions_per_chr(char* read_chr_id, 
	vector<t_mapped_read*>* cur_chr_reads,
	vector<t_mapped_read*>* no_overlap_reads,
	vector<t_annot_region*>* regions_2_exclude);

#endif // __MAPPED_READ_FILE_INTERFACE__

