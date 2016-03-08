#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "ms_min_max_utils.h"
#include "ms_filter_utils.h"
#include "ms_utils.h"
#include "ms_combinatorics.h"
#include "ms_xlog_math.h"
#include "ms_ansi_string.h"
#include "ms_annot_region_tools.h"
#include <list>

bool __DUMP_MIN_MAX_UTILS_MSGS__ = false;

void prune_extrema_derivative_map(vector<t_extrema_node*>* mins,
	vector<t_extrema_node*>* maxes,
	signed int* derivative_map, 
	int l_smallest_side,
	int l_signal, 
	vector<t_extrema_node*>* pruned_mins,
	vector<t_extrema_node*>* pruned_maxes,
	signed int* pruned_derivative_map)
{
	vector<t_extrema_node*>* extrema = new vector<t_extrema_node*>();
	extrema->insert(extrema->end(), mins->begin(), mins->end());
	extrema->insert(extrema->end(), maxes->begin(), maxes->end());

	sort(extrema->begin(), extrema->end(), sort_extremas_per_posn);

	// Copy pruned derivative map.
	for(int i = 0; i < l_signal; i++)
	{
		pruned_derivative_map[i] = derivative_map[i];
	}	
	
	int i_ext = 0;
	while(i_ext < (int)extrema->size()-1)
	{
		int l_cur_side = extrema->at(i_ext+1)->extrema_posn - extrema->at(i_ext)->extrema_posn;

		if(l_cur_side < l_smallest_side)
		{
			if(extrema->at(i_ext)->extrema_type == EXTREMA_MAX)
			{
				// This is from MAX to MIN: Fix the derivative map, do not add these extrema to the pruned list.
if(__DUMP_MIN_MAX_UTILS_MSGS__)
				fprintf(stderr, "Pruning MAX-MIN transition at %d, %d\n",extrema->at(i_ext)->extrema_posn, extrema->at(i_ext+1)->extrema_posn);
				for(int i = extrema->at(i_ext)->extrema_posn; i <= extrema->at(i_ext+1)->extrema_posn; i++)
				{
					if(derivative_map[i] > 0)
					{
						fprintf(stderr, "The derivative is set to negative from max to min transition.\n");
						exit(0);
					}
					else
					{
						pruned_derivative_map[i] = 1;
					}
				} // i loop.
			}
			else if(extrema->at(i_ext)->extrema_type == EXTREMA_MIN)
			{
if(__DUMP_MIN_MAX_UTILS_MSGS__)
				fprintf(stderr, "Pruning MIN-MAX transition at %d, %d\n",extrema->at(i_ext)->extrema_posn, extrema->at(i_ext+1)->extrema_posn);
				// This is from MIN to MAX: Fix the derivative map, do not add these extrema to the pruned list.
				for(int i = extrema->at(i_ext)->extrema_posn; i <= extrema->at(i_ext+1)->extrema_posn; i++)
				{
					if(derivative_map[i] < 0)
					{
						fprintf(stderr, "The derivative is set to negative from max to min transition.\n");
						exit(0);
					}
					else
					{
						pruned_derivative_map[i] = -1;
					}
				} // i loop.
			}

			// Remove the extrema at i_ext and i_ext+1.
			extrema->erase(extrema->begin()+i_ext);
			extrema->erase(extrema->begin()+i_ext);

			// At this point, we are basically at the next element.
		}
		else
		{
			// This side's length is ok. Move to the next side.
			i_ext++;
		}
	} // i_ext loop.

	// Must check whether the last extrema nodes are removed, then make sure those are added for correct tracing.
	if(pruned_derivative_map[0] != 0)
	{
		fprintf(stderr, "The first extrema is removed, adding it back.\n");
		//getc(stdin);

		t_extrema_node* begin_ext_node = new t_extrema_node();
		begin_ext_node->extrema_posn = 0;
		begin_ext_node->extrema_type = EXTREMA_MIN;
		begin_ext_node->higher_node = NULL;
		begin_ext_node->lower_node = NULL;
		begin_ext_node->scale = extrema->at(0)->scale;

		if(pruned_derivative_map[0] > 0)
		{
			// Add a minima to the beginning.
			begin_ext_node->extrema_type = EXTREMA_MIN;
		}
		else if(pruned_derivative_map[0] < 0)
		{
			// Add a maxima to the beginning.
			begin_ext_node->extrema_type = EXTREMA_MAX;
		}

		pruned_derivative_map[0] = 0;
		extrema->push_back(begin_ext_node);
	}

	if(pruned_derivative_map[l_signal-1] != 0)
	{
		fprintf(stderr, "The last extrema is removed, adding it back.\n");
		//getc(stdin);

		t_extrema_node* end_ext_node = new t_extrema_node();
		end_ext_node->extrema_posn = l_signal-1;
		end_ext_node->extrema_type = EXTREMA_MIN;
		end_ext_node->higher_node = NULL;
		end_ext_node->lower_node = NULL;
		end_ext_node->scale = extrema->back()->scale;

		if(pruned_derivative_map[l_signal-1] > 0)
		{
			end_ext_node->extrema_type = EXTREMA_MAX;
		}
		else if(pruned_derivative_map[l_signal-1] < 0)
		{
			end_ext_node->extrema_type = EXTREMA_MIN;
		}

		pruned_derivative_map[l_signal-1] = 0;
		extrema->push_back(end_ext_node);
	}

	// Sort the extrema.
	sort(extrema->begin(), extrema->end(), sort_extremas_per_posn);

	// Separate the total list into two lists.
	for(i_ext = 0; i_ext < (int)extrema->size(); i_ext++)
	{
		if(extrema->at(i_ext)->extrema_type == EXTREMA_MIN)
		{
			pruned_mins->push_back(extrema->at(i_ext));
		}
		else if(extrema->at(i_ext)->extrema_type == EXTREMA_MAX)
		{
			pruned_maxes->push_back(extrema->at(i_ext));
		}
	} // i_ext loop.

	//vector<t_extrema_node*>* all_pruned_nodes = new vector<t_extrema_node*>();
	//all_pruned_nodes->insert(all_pruned_nodes->begin(), pruned_mins->begin(), pruned_mins->end());
	//all_pruned_nodes->insert(all_pruned_nodes->begin(), pruned_maxes->begin(), pruned_maxes->end());

	//sort(all_pruned_nodes->begin(), all_pruned_nodes->end(), sort_extremas);

	// Set the derivative map.
	//int i_cur_nuc = 0;
	for(int i_ext = 0; i_ext < (int)extrema->size(); i_ext++)
	{
		if(i_ext > 0 &&
			extrema->at(i_ext-1)->extrema_type == extrema->at(i_ext)->extrema_type)
		{
			fprintf(stderr, "Consecutive same extrema: %d, %d: %d\n", 
				extrema->at(i_ext-1)->extrema_posn, 
				extrema->at(i_ext)->extrema_posn, 
				extrema->at(i_ext)->extrema_type);
			exit(0);
		}
	} // i_ext loop.
}

int extrema_posn_accessor(void* extrema_ptr)
{
	t_extrema_node* ptr = (t_extrema_node*)extrema_ptr;
	return(ptr->extrema_posn);
}


void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale)
{
	double zero_deriv = pow(10.0, -10.0);

	get_extrema_per_plateaus(data, l_signal,
		maxes, 
		mins,
		derivative_sign_map,
		i_scale,
		zero_deriv);
}

void get_extrema_per_plateaus(double* data, int l_signal,
	vector<t_extrema_node*>* maxes, 
	vector<t_extrema_node*>* mins,
	int* derivative_sign_map,
	int i_scale,
	double zero_deriv)
{
	vector<t_extrema_node*>* extrema_nodes = new vector<t_extrema_node*>();

	double* deriv = new double[l_signal + 2];
	//vector<int>* plateau_locs = new vector<int>();
	//// Following can be implemented with 
	//deriv[0] = 0;
	//for(int i = 1; i < l_signal; i++)
	//{
	//	deriv[i] = data[i] - data[i-1];

	//	if(fabs(deriv[i]) < zero_deriv)
	//	{
	//		//fprintf(stderr, "Pleateau loc: %d\n", i);
	//		plateau_locs->push_back(i);
	//	}
	//	else
	//	{
	//	}
	//} // i loop.

	//// Form the plateau signal.
	//vector<t_plateau*>* plateaus = new vector<t_plateau*>();

	//int* plateau_signal = new int[l_signal];
	//memset(plateau_signal, 0, sizeof(int) * l_signal);
	//int n_pl_pts = (int)plateau_locs->size(); 
	//int i = 0; 
	//while(i < n_pl_pts)
	//{
	//	int cur_start = plateau_locs->at(i);
	//	int cur_end = plateau_locs->at(i);

	//	while(i+1 < n_pl_pts &&
	//		plateau_locs->at(i+1) == plateau_locs->at(i)+1)
	//	{			
	//		cur_end = plateau_locs->at(i+1);

	//		i++;
	//	} // plateau finder loop.

	//	i++;

	//	t_plateau* new_plateau = new t_plateau();
	//	new_plateau->start = cur_start;
	//	new_plateau->end = cur_end;
	//	//fprintf(stderr, "Pleateau: %d-%d\n", cur_start, cur_end);

	//	plateaus->push_back(new_plateau);
	//} // i loop.	

	deriv[0] = 0;
	for(int i = 1; i <= l_signal; i++)
	{
			deriv[i] = data[i] - data[i-1];
	} // i loop.

	vector<t_plateau*>* plateaus = new vector<t_plateau*>();

	int i = 1;
	int cur_start = -1;
	int cur_end = -1;
	while(i <= l_signal)
	{
		// Are we starting a new plateau here?
		if(fabs(deriv[i]) < zero_deriv)
		{
			// Plateau: Get the beginning and the end.
			cur_start = i;
			cur_end = i;
			while(i <= l_signal &&
					fabs(deriv[i]) < zero_deriv)
			{
				cur_end = i;
				i++;
			} // i loop.

			t_plateau* new_plateau = new t_plateau();
			new_plateau->start = cur_start;
			new_plateau->end = cur_end;
			plateaus->push_back(new_plateau);
		}
		else
		{
			i++;
		}
	} // i loop.

	// Go over all the plateaus and set the max/min ver plateau.
	int l_win = 1;
	for(int i_pl = 0; i_pl < (int)plateaus->size(); i_pl++)
	{
		if(plateaus->at(i_pl)->start > l_win &&
			plateaus->at(i_pl)->end+l_win < l_signal &&
			deriv[plateaus->at(i_pl)->start-l_win] > 0 &&
			deriv[plateaus->at(i_pl)->end+l_win] < 0)
		{
			t_extrema_node* new_ext_node = new t_extrema_node();
			new_ext_node->extrema_posn = (plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2;
			new_ext_node->extrema_type = EXTREMA_MAX;
			new_ext_node->higher_node = NULL;
			new_ext_node->lower_node = NULL;
			new_ext_node->scale = i_scale;
			new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
			extrema_nodes->push_back(new_ext_node);
			//fprintf(stderr, "Adding max\n");
			//p2n->push_back((plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2);
		}

		if(plateaus->at(i_pl)->start > l_win &&
			plateaus->at(i_pl)->end+l_win < l_signal &&
			deriv[plateaus->at(i_pl)->start-l_win] < 0 &&
			deriv[plateaus->at(i_pl)->end+l_win] > 0)
		{
			t_extrema_node* new_ext_node = new t_extrema_node();
			new_ext_node->extrema_posn = (plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2;
			new_ext_node->extrema_type = EXTREMA_MIN;
			new_ext_node->higher_node = NULL;
			new_ext_node->lower_node = NULL;
			new_ext_node->scale = i_scale;
			new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
			extrema_nodes->push_back(new_ext_node);
			//fprintf(stderr, "Adding max\n");
			//n2p->push_back((plateaus->at(i_pl)->start + plateaus->at(i_pl)->end) / 2);
		}
	} // i_pl loop.

	// Add the sudden changes.
	for(int i = 1; i < l_signal; i++)
	{
		if(i+1 < l_signal)
		{
			if(deriv[i] > zero_deriv &&
				deriv[i+1] < -1 * zero_deriv)
			{
				t_extrema_node* new_ext_node = new t_extrema_node();
				new_ext_node->extrema_posn = i;
				new_ext_node->extrema_type = EXTREMA_MAX;
				new_ext_node->higher_node = NULL;
				new_ext_node->lower_node = NULL;
				new_ext_node->scale = i_scale;
				new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
				extrema_nodes->push_back(new_ext_node);
				//p2n->push_back(i);
			}

			if(deriv[i] < -1 * zero_deriv &&
				deriv[i+1] > zero_deriv)
			{
				t_extrema_node* new_ext_node = new t_extrema_node();
				new_ext_node->extrema_posn = i;
				new_ext_node->extrema_type = EXTREMA_MIN;
				new_ext_node->higher_node = NULL;
				new_ext_node->lower_node = NULL;
				new_ext_node->scale = i_scale;
				new_ext_node->height_at_extrema = data[new_ext_node->extrema_posn];
				extrema_nodes->push_back(new_ext_node);
				//n2p->push_back(i);
			}
		}
	} // i loop.

	if((int)extrema_nodes->size() == 0 ||
		(int)extrema_nodes->size() == 1)
	{
if(__DUMP_MIN_MAX_UTILS_MSGS__)
{
		fprintf(stderr, "Found %d extrema.\n", (int)extrema_nodes->size());
}
		return;
	}

	sort(extrema_nodes->begin(), extrema_nodes->end(), sort_extremas_per_posn);

	t_extrema_node* begin_ext_node = new t_extrema_node();
	begin_ext_node->extrema_posn = 0;
	begin_ext_node->extrema_type = EXTREMA_MIN;
	begin_ext_node->higher_node = NULL;
	begin_ext_node->lower_node = NULL;
	begin_ext_node->scale = i_scale;
	begin_ext_node->height_at_extrema = data[begin_ext_node->extrema_posn];

	t_extrema_node* end_ext_node = new t_extrema_node();
	end_ext_node->extrema_posn = l_signal-1;
	end_ext_node->extrema_type = EXTREMA_MIN;
	end_ext_node->higher_node = NULL;
	end_ext_node->lower_node = NULL;
	end_ext_node->scale = i_scale;
	end_ext_node->height_at_extrema = data[end_ext_node->extrema_posn];

	if(extrema_nodes->at(0)->extrema_type == EXTREMA_MAX)
	{
		begin_ext_node->extrema_type = EXTREMA_MIN;
	}
	else
	{
		begin_ext_node->extrema_type = EXTREMA_MAX;
	}

	if(extrema_nodes->back()->extrema_type == EXTREMA_MAX)
	{
		end_ext_node->extrema_type = EXTREMA_MIN;
	}
	else
	{
		end_ext_node->extrema_type = EXTREMA_MAX;
	}

	extrema_nodes->push_back(begin_ext_node);
	extrema_nodes->push_back(end_ext_node);

	sort(extrema_nodes->begin(), extrema_nodes->end(), sort_extremas_per_posn);

	// Set the derivative map.
	int i_cur_nuc = 0;
	for(int i_ext = 0; i_ext < (int)extrema_nodes->size(); i_ext++)
	{
		if(i_ext > 0 &&
			extrema_nodes->at(i_ext-1)->extrema_type == extrema_nodes->at(i_ext)->extrema_type)
		{
			fprintf(stderr, "Consecutive same extrema: %d, %d: %d\n", 
				extrema_nodes->at(i_ext-1)->extrema_posn, 
				extrema_nodes->at(i_ext)->extrema_posn, 
				extrema_nodes->at(i_ext)->extrema_type);
			exit(0);
		}

		int der_val = 0;
		if(extrema_nodes->at(i_ext)->extrema_type == EXTREMA_MAX)
		{
			der_val = 1;
		}
		else
		{
			der_val = -1;
		}

		for(int i_nuc = i_cur_nuc; i_nuc <= extrema_nodes->at(i_ext)->extrema_posn-1; i_nuc++)
		{
			derivative_sign_map[i_nuc] = der_val;
		} // i_nuc loop.

		derivative_sign_map[extrema_nodes->at(i_ext)->extrema_posn] = 0;
		i_cur_nuc = extrema_nodes->at(i_ext)->extrema_posn+1;
		//cur_type = extrema_nodes->at(i_ext)->extrema_type;
	} // i_ext loop.

if(__DUMP_MIN_MAX_UTILS_MSGS__)
{
	fprintf(stderr, "Found %d extrema.\n", (int)extrema_nodes->size());
}
	//for(int i = 0; i < 20; i++)
	//{
	//	fprintf(stderr, "%d: %d\n", extrema_nodes->at(i)->extrema_posn, extrema_nodes->at(i)->extrema_type);
	//}

	for(int i = 0; i < (int)extrema_nodes->size(); i++)
	{
		if(extrema_nodes->at(i)->extrema_type == EXTREMA_MAX)
		{
			maxes->push_back(extrema_nodes->at(i));
		}
		else
		{
			mins->push_back(extrema_nodes->at(i));
		}
	} // i loop.

	// Free memory. These are pretty large chunks of memory.
	delete [] deriv;
	//delete [] plateau_signal;
	for(int i_p = 0; i_p < (int)plateaus->size(); i_p++)
	{
		delete(plateaus->at(i_p));
	} // i_p loop.
	//delete plateau_locs;
}


bool sort_extremas_per_decreasing_height(t_extrema_node* node1, t_extrema_node* node2)
{
	return(node1->height_at_extrema > node2->height_at_extrema);
}

bool sort_extremas_per_height(t_extrema_node* node1, t_extrema_node* node2)
{
	return(node1->height_at_extrema < node2->height_at_extrema);
}

bool sort_extremas_per_posn(t_extrema_node* node1, t_extrema_node* node2)
{
	return(node1->extrema_posn < node2->extrema_posn);
}

void delete_extrema_nodes(vector<t_extrema_node*>* extrema_nodes)
{
	for(int i_n = 0; i_n < (int)extrema_nodes->size(); i_n++)
	{
		delete extrema_nodes->at(i_n);
	}

	delete extrema_nodes;
}

