#ifndef __GENOME_SEQUENCE_TOOLS__
#define __GENOME_SEQUENCE_TOOLS__

#include <vector>

using namespace std;

struct t_annot_region;

bool get_seq_i_per_genome_i(char strand, 
	int genome_start, 
	int genome_end, 
	int genome_i, int& seq_i);

bool get_genome_i_per_seq_i(char strand, 
	int genome_start, 
	int genome_end, 
	int& genome_i, int seq_i);

// Separate the chromosome data from genome files, dump it in binary format into output directory.
void binarize_genome_per_plain(char* fp, char* op_dir);
void binarize_fasta_file(char* fasta_fp, char* bin_dir);

char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq);

//char* numerize_sequence_signal(char* seq_signal, int l_seq);

void get_numerized_sequence_signals(char* sequence, int l_signal,
	double* walk_signal,
	double* nucleotide_signal,
	double* gc_content_signal);

bool check_kmer_composition(int* num_seq, int left_win_pos, int k);

void get_kmer_signal_per_numerized_sequence_signal(double* numerized_sequence_signal, int l_signal, double* kmer_signal, int k);

// Extract sequence from binary data.
char* extract_sequence(char* chr_data_dir, char* chr_name, int start, int end, char strand, bool all_upper_case = false);
char* load_whole_chromosome_sequence_per_bin_file(char* chr_fp, int& n_nucs);

void batch_extract_region_sequences(char* chr_data_dir, vector<t_annot_region*>* regions);
void batch_extract_interval_sequences(char* chr_data_dir, vector<t_annot_region*>* regions);

vector<t_annot_region*>* load_BED_with_sequences(char* bed_w_seq_fp);
vector<t_annot_region*>* load_BED_with_ancestral_derived_seqs(char* bed_w_seq_fp);

// Specifies the matching position of a k-mer on a list of sequences.
struct t_k_mer_match_posn_info
{
	int seq_i;
	int left_kmer_posn;
	int strand;
};

struct t_kmer_freq_info
{
	char* kmer;
	int count;
	int bckgrnd_count;
	vector<t_k_mer_match_posn_info*>* k_mer_match_posns;
};

struct t_kmer_count_tree_entry
{
	int k;
	t_kmer_count_tree_entry** next_nuc_entries;
	int count;
	int bckgrnd_count; // This is the number of times this kmer is seen in the background.
};

double get_max_match_per_kmers(char* kmer1, char* kmer2, int l_kmer1, int l_kmer2);

void delete_kmer_freq_tree(t_kmer_count_tree_entry* cur_tree_entry_2_delete);
int get_total_kmer_count(t_kmer_count_tree_entry* cur_tree_entry, int k_2_count);
void get_bckgrnd_counts(t_kmer_count_tree_entry* cur_tree_entry, t_kmer_count_tree_entry* cur_neg_tree_entry);
void get_kmer_count_array(vector<t_kmer_freq_info*>* cur_kmer_list, t_kmer_count_tree_entry* cur_entry, char* cur_kmer);
t_kmer_count_tree_entry* build_kmer_frequency_tree(vector<int*>* positive_num_seqs, vector<int>* positive_num_seq_lengths, int k);
void get_kmer_frequency_enrichments(vector<char*>* positive_sequences, vector<char*>* negative_sequences, 
								vector<int*>* positive_sequence_per_nuc_depths,
								int min_k, int max_k);

void update_kmer_frequency_tree(t_kmer_count_tree_entry* main_entry, 
	int* new_seq, 
	int l_cur_seq, 
	int k);

t_kmer_count_tree_entry* init_kmer_frequency_tree(int k);

void dump_seqs_freqs(FILE* f_op, t_kmer_count_tree_entry* cur_entry, char* cur_kmer);

char* get_forward_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max);
char* get_reverse_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max);

char* get_reverse_complement(char* cur_aln_line, int i_nuc_min, int i_nuc_max);

void get_shuffled_motif_matching_statistics(char* sequence, double** kmer_emission_probs, int k, double& mean, double& std_dev, int n_randomizations);
void scan_motif_in_sequence(char* sequence, double** kmer_emission_probs, int k, double* postv_match_score, double* negtv_match_score);

// Following are k-mer enrichment related functions.
void get_kmer_probability_profiles(vector<char*>* positive_sequences, vector<char*>* negative_sequences,	// These are the list of sequences to scan.
									vector<int*>* positive_sequence_per_nuc_depths,							// This is the per nucleotide score for the positive sequences.
									int min_k, int max_k, int n_inits, int n_iters_per_init,
									int max_matches_per_sequence,  // This is the number of matches to be traced per sequence.
									char* op_dir); // Output directory for output files.

void update_motif_counts_per_sequence_list(vector<int*>* numerized_sequences,
	vector<int>* numerized_seq_lengths,
	vector<int*>* per_nuc_depths,
	int min_k, int max_k, 
	double*** kmer_emission_log_probs,
	double*** kmer_observation_log_counts,
	bool positive_seq);

void dump_embedded_sequences_w_kmer(char* kmer, int l_seqs, int n_seqs, double mutation_prob, char* op_fp);
void dump_mutated_embedded_sequences_w_multiple_kmers(vector<char*>* kmers, vector<double>* kmer_fractions, int l_seqs, int n_seqs, double mutation_prob, char* op_fp);

#endif // __GENOME_SEQUENCE_TOOLS__

