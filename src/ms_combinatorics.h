#ifndef __PERMUTATIONS__
#define __PERMUTATIONS__

#include <vector>
using namespace std;

class t_rng;

bool is_equal(vector<int>* combinatorial_indices1, vector<int>* combinatorial_indices2);

bool is_gt(vector<int>* combinatorial_indices1, vector<int>* combinatorial_indices2);

bool get_next_combination_indices(vector<int>* combination_indices, 
								vector<int>* current_counting_indices,
								int n_all_elements_per_posn);

bool get_next_counting_indices(vector<int>* counting_indices, 
								vector<int>* n_all_elements_per_digit);

bool get_random_counting_indices(vector<vector<int>*>* random_counting_indices_list,
								int k,
								vector<int>* n_all_elements_per_digit,
								int n_indices_2_generate,
								t_rng* rng);

void dump_binomial_p_vals(int max_grand_total, double* log_factorials);
double** buffer_binomial_pvalues(int max_grand_total, double* log_factorials, int BINOMIAL_P_VAL_BIN_SIZE);
void delete_buffered_binomial_pvalues(double** binomial_p_vals, int max_grand_total, int BINOMIAL_P_VAL_BIN_SIZE);
double* buffer_log_factorials(int n);

#endif // __PERMUTATIONS__