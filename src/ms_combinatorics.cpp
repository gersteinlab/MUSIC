#include <stdio.h>
#include <stdlib.h>
#include "ms_combinatorics.h"
#include <algorithm>
#include "ms_xlog_math.h"
#include <algorithm>

double factorial(int n)
{
	double cur_fact = 1.0;
	for(int i = 1; i <= n; i++)
	{
		cur_fact *= i;
	} 

	return(cur_fact);
}

bool is_gt(vector<int>* combinatorial_indices1, vector<int>* combinatorial_indices2)
{
	for(int i = 0; i < (int)combinatorial_indices1->size(); i++)
	{
		if(combinatorial_indices1->at(i) > combinatorial_indices2->at(i))
		{
			return(true);
		}
		else if(combinatorial_indices1->at(i) == combinatorial_indices2->at(i))
		{
			
		}
	} // i loop.

	return(false);
}

bool is_equal(vector<int>* combinatorial_indices1, vector<int>* combinatorial_indices2)
{
	for(int i = 0; i < (int)combinatorial_indices1->size(); i++)
	{
		if(combinatorial_indices1->at(i) != combinatorial_indices2->at(i))
		{
			return(false);
		}
	} // i loop.

	return(true);
}

bool is_valid_indexing(vector<int>* sorted_combinatorial_indices)
{
	for(int i = 0; i < (int)sorted_combinatorial_indices->size()-1; i++)
	{
		if(sorted_combinatorial_indices->at(i) == sorted_combinatorial_indices->at(i+1))
		{
			return(false);
		}
	} // i loop.

	return(true);
}

bool is_increasing(vector<int>* sorted_combinatorial_indices)
{
	for(int i = 0; i < (int)sorted_combinatorial_indices->size()-1; i++)
	{
		if(sorted_combinatorial_indices->at(i) >= sorted_combinatorial_indices->at(i+1))
		{
			return(false);
		}
	} // i loop.

	return(true);
}

bool get_next_counting_indices(vector<int>* counting_indices, 
								vector<int>* n_all_elements_per_digit)
{
	bool updated_comb_indices = false;
	for(int i = (int)counting_indices->size()-1; i >= 0; i--)
	{
		if(counting_indices->at(i)+1 < n_all_elements_per_digit->at(i))
		{
			counting_indices->at(i) = counting_indices->at(i)+1;
			updated_comb_indices = true;
			break;
		}
		else
		{
			counting_indices->at(i) = 0;
		}
	} // i loop.

	// No more updates available.
	if(!updated_comb_indices)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

double* buffer_log_factorials(int n)
{
	double* factorials = new double[n+2];
	factorials[0] = xlog(1.0);
	for(int i = 1; i <= n; i++)
	{
		factorials[i] = xlog_mul(factorials[i-1], xlog(i));
	} // i loop.

	return(factorials);
}

bool get_next_combination_indices(vector<int>* combination_indices, 
								vector<int>* current_counting_indices,
								int n_all_elements)
{
	while(1)
	{
		// Update the last index in the sorted list.
		bool updated_comb_indices = false;
		for(int i = (int)current_counting_indices->size()-1; i >= 0; i--)
		{
			if(current_counting_indices->at(i)+1 < n_all_elements)
			{
				current_counting_indices->at(i) = current_counting_indices->at(i)+1;
				updated_comb_indices = true;
				break;
			}
			else
			{
				current_counting_indices->at(i) = 0;
			}
		} // i loop.

		// No more updates available.
		if(!updated_comb_indices)
		{
			return(false);
		}

		// Sort the current indices.
		if(is_increasing(current_counting_indices))
		{	
			break;
		}
	} // Get a new permutation.

	// Copy the indices.
	for(int i = 0; i < (int)current_counting_indices->size(); i++)
	{
		combination_indices->at(i) = current_counting_indices->at(i);
	} // i loop.

	return(true);
}

