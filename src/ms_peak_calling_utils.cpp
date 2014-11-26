#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "ms_signal_track_tools.h"
#include "ms_ansi_string.h"
#include "ms_utils.h"
#include "ms_annot_region_tools.h"
#include "ms_genomics_coords.h"
#include "ms_nomenclature.h"
#include "ms_xlog_math.h"
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include "ms_min_max_utils.h"
#include "ms_gsl_fft_filter_utils.h"
#include "ms_mapped_read_tools.h"
#include "ms_profile_normalization.h"
#include "ms_combinatorics.h"
#include "ms_utils.h"
#include "ms_signal_enrichment_utils.h"
#include "ms_peak_calling_utils.h"

bool __DUMP_PEAK_CALLING_UTILS_MSGS__ = false;

void write_decomposition_bedGraphs(char* chip_reads_dir,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_mapability_filtering_win,
		double max_normalized_mapability_signal)
{
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", chip_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's for the chip data.\n");
		exit(0);
	}

	// Process each chromosome.
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Writing %s..", chr_ids->at(i_chr));

		// Skip chromosome M for now.
		if(strcmp(chr_ids->at(i_chr), "M") == 0)
		{
			continue;
		}

		// Generate the current signal profile.
		int l_buffer = 300*1000*1000;
		int l_profile = 0;
		//int l_read = l_read_mapability_signal;
		double* buffered_signal_profile = new double[l_buffer + 2];	
		char cur_chr_chip_reads_fp[1000];
		sprintf(cur_chr_chip_reads_fp, "%s/%s_mapped_reads.txt", chip_reads_dir, chr_ids->at(i_chr));
		buffer_per_nucleotide_profile_no_buffer(cur_chr_chip_reads_fp, l_fragment, 
			buffered_signal_profile, NULL, NULL,
			l_buffer, l_profile);

		double* signal_profile = new double[l_profile + 2];	
		for(int i = 1; i <= l_profile; i++)
		{
			signal_profile[i] = buffered_signal_profile[i];
		} // i loop.
		delete [] buffered_signal_profile;

		//// Generate the current control profile.
		//l_buffer = 300*1000*1000;
		//int l_control = 0;
		//double* buffered_control_profile = new double[l_buffer + 2];
		//char cur_chr_control_reads_fp[1000];
		//sprintf(cur_chr_control_reads_fp, "%s/%s_mapped_reads.txt", control_reads_dir, chr_ids->at(i_chr));
		//buffer_per_nucleotide_profile_no_buffer(cur_chr_control_reads_fp, l_fragment, 
		//	buffered_control_profile, NULL, NULL, 
		//	l_buffer, l_control);


		//double* control_profile = new double[l_control + 2];
		//for(int i = 1; i <= l_control; i++)
		//{
		//	control_profile[i] = buffered_control_profile[i];
		//} // i loop.
		//delete [] buffered_control_profile;
		
		// Load and process the mapability profile.
		char mapability_signal_profile_fp[1000];
		sprintf(mapability_signal_profile_fp, "%s/%s.bin", mapability_signal_dir, chr_ids->at(i_chr));

		double* mapability_aware_smoothed_signal_profile = NULL;
		int l_mapability_profile = 0;
		if(check_file(mapability_signal_profile_fp))
		{
			double* mapability_signal_profile = NULL;

#define __UCHAR_MAPPABILITY__
#ifdef __DOUBLE_MAPPABILITY__
			// Load the mapability map signal profile, do filtering.
			mapability_signal_profile = load_per_nucleotide_binary_profile(mapability_signal_profile_fp, l_mapability_profile);

			// Do mapability aware median filtering on the current signal profile.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
			fprintf(stderr, "Scaling the mapability map with %d.\n", l_read_mapability_signal * 2);

			int mapability_scaling = l_read * 2;
			for(int i = 1; i <= l_mapability_profile; i++)
			{
				mapability_signal_profile[i] /= mapability_scaling;
			} // i loop.
#elif defined(__UCHAR_MAPPABILITY__)
			// Following loads the mappability signal profile from the char version of the multi-mappability profile.
			// Load the mapability map signal profile, do filtering.
			unsigned char* mapability_signal_char_profile = load_per_nucleotide_binary_uchar_profile(mapability_signal_profile_fp, l_mapability_profile);
			mapability_signal_profile = new double[l_mapability_profile + 2];
			for(int i = 1; i <= l_mapability_profile; i++)
			{
				unsigned char unsigned_char_val = (unsigned char)(mapability_signal_char_profile[i]);
				mapability_signal_profile[i] = (double)(unsigned_char_val);
				mapability_signal_profile[i] /= 100;

				if(mapability_signal_profile[i] < 0)
				{
					fprintf(stderr, "Sanity check failed.\n");
					exit(0);
				}
			} // i loop.
			delete [] mapability_signal_char_profile;
#else
			#error "Must define the type of mappability."
#endif

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
			fprintf(stderr, "Smoothing the signal profile.\n");
			mapability_aware_smoothed_signal_profile = mapability_aware_median_filter(signal_profile, l_profile,
																					mapability_signal_profile, l_mapability_profile,
																					max_normalized_mapability_signal,
																					l_mapability_filtering_win);

			delete [] mapability_signal_profile;
		} // mapability based filtering check.
		else
		{
			fprintf(stderr, "Could not find the mapability map profile, skipping it.");
			mapability_aware_smoothed_signal_profile = new double[l_profile + 2];
			for(int i = 1; i <= l_profile; i++)
			{
				mapability_aware_smoothed_signal_profile[i] = signal_profile[i];
			} // i loop.
		}

//		// Scale the control with respect to signal.
//		double per_win_2DOF_lls_scaling_factor = 0;
//		double per_win_1DOF_lls_scaling_factor = 0;
//		double total_sig_scaling_factor = 0;
//		get_per_window_scaling_factors_for_profile1_per_profile2(control_profile, l_control, signal_profile, l_profile, 10000,
//																	per_win_2DOF_lls_scaling_factor,
//																	per_win_1DOF_lls_scaling_factor,
//																	total_sig_scaling_factor);
//
//		double scaling_factor = per_win_1DOF_lls_scaling_factor;
//
//if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
//		fprintf(stderr, "Scaling control with factor of %lf\n", scaling_factor);
//
//		// Go over the whole control signal.
//		for(int i = 0; i < l_control; i++)
//		{
//			control_profile[i] *= scaling_factor;
//		} // i loop.

		// Now, do the multiscale decomposition, 
		vector<double>* scales_per_i_scale = new vector<double>();
		vector<double*>* filtered_tracks = multiscale_median_filter_data(mapability_aware_smoothed_signal_profile, l_profile, 
			base_scale_l_win, end_scale_l_win, log_step, 
			scales_per_i_scale, 
			true,
			false, 
			false, 
			false, 
			NULL, 
			NULL,
			"temp");

		// Dump the tracks.
		for(int i_scale = 0; i_scale < (int)scales_per_i_scale->size(); i_scale++)
		{
			char cur_decomp_fp[1000];
			sprintf(cur_decomp_fp, "%s_scale_%d.bin", "temp", i_scale);

			if(check_file(cur_decomp_fp))
			{
				int l_cur_decomp;
				double* cur_filtered_track = load_per_nucleotide_binary_profile(cur_decomp_fp, l_cur_decomp);

				// Delete the temp file.
				remove(cur_decomp_fp);

				char cur_bgr_fp[1000];
				sprintf(cur_bgr_fp, "chr%s_scale_%d.bgr", chr_ids->at(i_chr), (int)(floor(scales_per_i_scale->at(i_scale))));
				dump_bedGraph_per_per_nucleotide_binary_profile(cur_filtered_track, l_cur_decomp, chr_ids->at(i_chr), cur_bgr_fp);
				delete [] cur_filtered_track;
			}
		} // i_scale loop.
		delete filtered_tracks;

		delete scales_per_i_scale;
		delete [] signal_profile;
		//delete [] control_profile;
		delete [] mapability_aware_smoothed_signal_profile;
	} // i_chr loop.
}

void get_peaks(char* chip_reads_dir,
		char* control_reads_dir,
		int l_fragment,
		char* mapability_signal_dir,
		int l_read_mapability_signal,
		double base_scale_l_win,
		double end_scale_l_win,
		double log_step,
		int l_p_val_norm_win,
		int l_mapability_filtering_win,
		double max_normalized_mapability_signal,
		double filtered_vs_non_filtered_max_scaler,
		double signal_profile_scaling_factor,
		double rd_end_pruning_p_val,
		double min_per_strand_evennes_fraction,
		int p_value_computation_type,
		bool do_replace_profiles_w_smoothed_profiles,
		bool do_BJ_correction_on_minima,
		bool do_filter_minima_per_length_min_base_scale_l_win,
		bool do_post_peak_p_value_minimization,
		double p_val_min_pruning_flip_prob,
		bool find_mappable_troughs,
		double q_value_cutoff)
{
	char test_op_bed_fp[1000];
	sprintf(test_op_bed_fp, "ERs_%.1f_%.1f_%.2f_%d_%.1f.bed", base_scale_l_win, end_scale_l_win, log_step, l_p_val_norm_win, filtered_vs_non_filtered_max_scaler);
	if(check_file(test_op_bed_fp))
	{
		fprintf(stderr, "Found the output file %s, will not overwrite, exiting.\n", test_op_bed_fp);
		exit(0);
	}

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", chip_reads_dir);
	vector<char*>* chr_ids = buffer_file(chr_ids_fp);
	if(chr_ids == NULL)
	{
		fprintf(stderr, "Could not load the chromosome id's for the chip data.\n");
		exit(0);
	}

	// Process each chromosome.
	vector<t_annot_region*>* all_peaks = new vector<t_annot_region*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "%s..", chr_ids->at(i_chr));

		// Skip chromosome M for now.
		if(strcmp(chr_ids->at(i_chr), "M") == 0)
		{
			continue;
		}

		// Generate the current signal profile.
		int l_buffer = 300*1000*1000;
		int l_profile = 0;
		int l_read = l_read_mapability_signal;
		double* buffered_signal_profile = new double[l_buffer + 2];	
		char cur_chr_chip_reads_fp[1000];
		sprintf(cur_chr_chip_reads_fp, "%s/%s_mapped_reads.txt", chip_reads_dir, chr_ids->at(i_chr));

		// Check IP file.
		if(!check_file(cur_chr_chip_reads_fp))
		{
			fprintf(stderr, "Processed control reads for %s does not exist @ %s. Skipping.\n", chr_ids->at(i_chr), cur_chr_chip_reads_fp);
			continue;
		}

		buffer_per_nucleotide_profile_no_buffer(cur_chr_chip_reads_fp, l_fragment, 
			buffered_signal_profile, NULL, NULL,
			l_buffer, l_profile);

		double* signal_profile = new double[l_profile + 2];	
		for(int i = 1; i <= l_profile; i++)
		{
			signal_profile[i] = buffered_signal_profile[i];
		} // i loop.
		delete [] buffered_signal_profile;

		// Generate the current control profile.
		l_buffer = 300*1000*1000;
		int l_control = 0;
		double* buffered_control_profile = new double[l_buffer + 2];
		char cur_chr_control_reads_fp[1000];
		sprintf(cur_chr_control_reads_fp, "%s/%s_mapped_reads.txt", control_reads_dir, chr_ids->at(i_chr));
		buffer_per_nucleotide_profile_no_buffer(cur_chr_control_reads_fp, l_fragment, 
			buffered_control_profile, NULL, NULL, 
			l_buffer, l_control);

		// Check control file.
		if(!check_file(cur_chr_control_reads_fp))
		{
			fprintf(stderr, "Processed control reads for %s does not exist @ %s. Skipping.\n", chr_ids->at(i_chr), cur_chr_control_reads_fp);
			continue;
		}

		double* control_profile = new double[l_control + 2];
		for(int i = 1; i <= l_control; i++)
		{
			control_profile[i] = buffered_control_profile[i];
		} // i loop.
		delete [] buffered_control_profile;
		
		// Load and process the mapability profile.
		char mapability_signal_profile_fp[1000];
		sprintf(mapability_signal_profile_fp, "%s/%s.bin", mapability_signal_dir, chr_ids->at(i_chr));

		double* mapability_aware_smoothed_signal_profile = NULL;
		int l_mapability_profile = 0;
		if(check_file(mapability_signal_profile_fp))
		{
			double* mapability_signal_profile = NULL;

#define __UCHAR_MAPPABILITY__
#ifdef __DOUBLE_MAPPABILITY__
			// Load the mapability map signal profile, do filtering.
			mapability_signal_profile = load_per_nucleotide_binary_profile(mapability_signal_profile_fp, l_mapability_profile);

			// Do mapability aware median filtering on the current signal profile.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
			fprintf(stderr, "Scaling the mapability map with %d.\n", l_read * 2);

			int mapability_scaling = l_read * 2;
			for(int i = 1; i <= l_mapability_profile; i++)
			{
				mapability_signal_profile[i] /= mapability_scaling;
			} // i loop.
#elif defined(__UCHAR_MAPPABILITY__)
			// Following loads the mappability signal profile from the char version of the multi-mappability profile.
			// Load the mapability map signal profile, do filtering.
			unsigned char* mapability_signal_char_profile = load_per_nucleotide_binary_uchar_profile(mapability_signal_profile_fp, l_mapability_profile);
			mapability_signal_profile = new double[l_mapability_profile + 2];
			for(int i = 1; i <= l_mapability_profile; i++)
			{
				unsigned char unsigned_char_val = (unsigned char)(mapability_signal_char_profile[i]);
				mapability_signal_profile[i] = (double)(unsigned_char_val);
				mapability_signal_profile[i] /= 100;

				if(mapability_signal_profile[i] < 0)
				{
					fprintf(stderr, "Sanity check failed.\n");
					exit(0);
				}
			} // i loop.
			delete [] mapability_signal_char_profile;
#else
			#error "Must define the type of mappability."
#endif

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
			fprintf(stderr, "Smoothing the signal profile.\n");

			mapability_aware_smoothed_signal_profile = mapability_aware_median_filter(signal_profile, l_profile,
																					mapability_signal_profile, l_mapability_profile,
																					max_normalized_mapability_signal,
																					l_mapability_filtering_win);

			delete [] mapability_signal_profile;
		} // mapability based filtering check.
		else
		{
			fprintf(stderr, "Could not find the mapability map profile, skipping it.");
			mapability_aware_smoothed_signal_profile = new double[l_profile + 2];
			for(int i = 1; i <= l_profile; i++)
			{
				mapability_aware_smoothed_signal_profile[i] = signal_profile[i];
			} // i loop.
		}

		// Signal profile scaling loop.
		for(int i = 1; i <= l_profile; i++)
		{
			signal_profile[i] *= signal_profile_scaling_factor;
		} // i loop.

		// Scale the control with respect to signal.
		double per_win_2DOF_lls_scaling_factor = 0;
		double per_win_1DOF_lls_scaling_factor = 0;
		double total_sig_scaling_factor = 0;
		get_per_window_scaling_factors_for_profile1_per_profile2(control_profile, l_control, signal_profile, l_profile, 10000,
																	per_win_2DOF_lls_scaling_factor,
																	per_win_1DOF_lls_scaling_factor,
																	total_sig_scaling_factor);

		double scaling_factor = per_win_1DOF_lls_scaling_factor;

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
		fprintf(stderr, "Scaling control with factor of %lf\n", scaling_factor);

		// Go over the whole control signal.
		for(int i = 0; i < l_control; i++)
		{
			control_profile[i] *= scaling_factor;
		} // i loop.

		// Dump the signal and control profiles after fixing.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{			
		char mapable_aware_filtered_control_profile_fp[] = "mapability_aware_profile.bin";
		dump_per_nucleotide_binary_profile(mapability_aware_smoothed_signal_profile, l_profile, mapable_aware_filtered_control_profile_fp);
}

		// Do multiscale decomposition and get the minima.
		vector<double>* scales_per_i_scale = new vector<double>();
		vector<vector<t_annot_region*>*>* per_scale_minima_regions = new vector<vector<t_annot_region*>*>();

		get_filtered_maxima_regions_multiscale_filtered_data(mapability_aware_smoothed_signal_profile, 
			l_profile, 
			base_scale_l_win, end_scale_l_win, log_step,
			scales_per_i_scale,
			per_scale_minima_regions,
			filtered_vs_non_filtered_max_scaler);

		delete [] mapability_aware_smoothed_signal_profile;

if(__DUMP_PEAK_CALLING_UTILS_MSGS__ )
{
			FILE* f_scales = open_f("per_scale_l_win.txt", "w");
			for(int i_scale = 0; i_scale < (int)scales_per_i_scale->size(); i_scale++)
			{
				fprintf(f_scales, "%lf\n", scales_per_i_scale->at(i_scale));
			}
			fclose(f_scales);
}

		// Filter the minima per p-value and fdr based end pruning.
		vector<t_annot_region*>* all_filtered_pruned_minima_regions = new vector<t_annot_region*>();
		for(int i_scale = 0; i_scale < (int)per_scale_minima_regions->size(); i_scale++)
		{
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{
			fprintf(stderr, "Filtering minima at %d. scale.\n", i_scale);
}

			if(scales_per_i_scale->at(i_scale) < base_scale_l_win)
			{
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
				fprintf(stderr, "Skipping filtering %d. scale.\n", i_scale);
			}
			else
			{
if(__DUMP_PEAK_CALLING_UTILS_MSGS__ )
{
				// Dump the filtered minima.
				char cur_scaled_filtered_minima_fp[1000];
				sprintf(cur_scaled_filtered_minima_fp, "all_minima_scale_%d.bed", i_scale);
				dump_BED(cur_scaled_filtered_minima_fp, per_scale_minima_regions->at(i_scale));
}

				vector<t_annot_region*>* rd_pruned_minima = prune_region_ends_per_window_average(per_scale_minima_regions->at(i_scale), signal_profile, l_profile, 
					control_profile, l_control, 
					rd_end_pruning_p_val, 1000*1000);

				//fprintf(stderr, "P-value minimization based end pruning.\n");
				//vector<t_annot_region*>* p_value_end_pruned_minima = prune_region_ends_per_p_value_minimization(rd_pruned_minima, signal_profile, l_profile, control_profile, l_control, l_fragment);
				//delete_annot_regions(rd_pruned_minima);

				//// Rename this.
				//rd_pruned_minima = p_value_end_pruned_minima;

				// Get the grand total to fasten up the p-value based filtering of the minimas.
				double max_grand_total_int = 0.0;
				for(int i_reg = 0; i_reg < (int)rd_pruned_minima->size(); i_reg++)
				{
					// Get the extrema signal values.
					double current_grand_total = 0.0;
					for(int i = rd_pruned_minima->at(i_reg)->start; i <= rd_pruned_minima->at(i_reg)->end; i++)
					{
						if(i < l_profile)
						{
							current_grand_total += signal_profile[i];
						}

						if(i < l_control)
						{
							current_grand_total += control_profile[i];
						}
					} // i loop.

					if(max_grand_total_int < current_grand_total)
					{
						max_grand_total_int = (int)(floor(current_grand_total));
					} 
				} // i_reg loop.

				double* log_factorials = buffer_log_factorials((int)(max_grand_total_int+20));
				
				// Get the p-values and filter the peaks.
				for(int i_reg = 0; i_reg < (int)rd_pruned_minima->size(); i_reg++)
				{
					// Get the p-value for the current region.
					double total_profile_signal = 0;
					double total_control_signal = 0;
					double cur_region_log_p_val = 0.0;
				
					if(rd_pruned_minima->at(i_reg)->end < l_control && 
						rd_pruned_minima->at(i_reg)->end < l_profile)
					{
						for(int i = rd_pruned_minima->at(i_reg)->start; i <= rd_pruned_minima->at(i_reg)->end; i++)
						{
							total_profile_signal += floor(signal_profile[i]);
							total_control_signal += floor(control_profile[i]);	
						} // i loop.

						cur_region_log_p_val = get_binomial_pvalue_per_region_neighborhood_window_control(signal_profile, control_profile, 
																			rd_pruned_minima->at(i_reg)->start, rd_pruned_minima->at(i_reg)->end,
																			log_factorials,
																			l_fragment, 
																			true, 
																			l_p_val_norm_win, 
																			l_profile, 
																			l_control,
																			1);

						rd_pruned_minima->at(i_reg)->significance_info = new t_significance_info();
						rd_pruned_minima->at(i_reg)->significance_info->log_p_val = cur_region_log_p_val;
					}
					else
					{
						rd_pruned_minima->at(i_reg)->significance_info = new t_significance_info();
						rd_pruned_minima->at(i_reg)->significance_info->log_p_val = xlog(1.0);
					}
				} // i_reg loop.

				// Delete the buffered factorials.
				delete [] log_factorials;

				// Do the benj. hochberg correction on the p-values.
				vector<t_annot_region*>* cur_scale_filtered_minima = new vector<t_annot_region*>();

				// Should we do BJ correction on the minima p-values?
				if(do_BJ_correction_on_minima)
				{
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
					fprintf(stderr, "Doing Benjamini-Hochberg correction on the filtered minima regions.\n");

					get_benjamini_hochberg_corrected_p_values(rd_pruned_minima);
				} // do_BJ_correction_on_minima

				// Do the p-value filtering on the minima.
				for(int i_reg = 0; i_reg < (int)rd_pruned_minima->size(); i_reg++)
				{
					// Check if this minima is significant.
					bool is_minima_significant = (rd_pruned_minima->at(i_reg)->significance_info->log_p_val < xlog(0.05));;

					// If BJ correction on the minima is requested, look at the q-values for filtering.
					if(do_BJ_correction_on_minima)
					{
						is_minima_significant = (rd_pruned_minima->at(i_reg)->significance_info->log_q_val < xlog(0.05));
					}

					// Use only the regions that are significant.
					if(is_minima_significant &&
						(!do_filter_minima_per_length_min_base_scale_l_win || (rd_pruned_minima->at(i_reg)->end - rd_pruned_minima->at(i_reg)->start > base_scale_l_win)))
					{
						cur_scale_filtered_minima->push_back(rd_pruned_minima->at(i_reg));
					}
					else
					{
						//delete rd_pruned_minima->at(i_reg)->significance_info;
						//delete_annot_regions(rd_pruned_minima->at(i_reg));
					}
				} // i_reg loop.					

				// Do p-value minimization based end pruning.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
				fprintf(stderr, "P-value minimization based end pruning.\n");

				vector<t_annot_region*>* p_value_end_pruned_minima = prune_region_ends_per_p_value_minimization(cur_scale_filtered_minima, signal_profile, l_profile, control_profile, l_control, l_fragment);

				// Delete the rd pruned minima regions: cur_scale_filtered_minima is a subset of rd_pruned_minima and is not copied, do not delete it.
				for(int i_reg = 0; i_reg < (int)rd_pruned_minima->size(); i_reg++)
				{
					delete rd_pruned_minima->at(i_reg)->significance_info;
				} // i_reg loop.
				delete_annot_regions(rd_pruned_minima);

				all_filtered_pruned_minima_regions->insert(all_filtered_pruned_minima_regions->end(), 
					p_value_end_pruned_minima->begin(), p_value_end_pruned_minima->end());

				// Replace the chromosome ids for each region.
				for(int i_reg = 0; i_reg < (int)p_value_end_pruned_minima->size(); i_reg++)
				{
					delete [] p_value_end_pruned_minima->at(i_reg)->chrom;
					p_value_end_pruned_minima->at(i_reg)->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
				} // i_reg loop.

				char cur_scaled_filtered_minima_fp[1000];
				sprintf(cur_scaled_filtered_minima_fp, "SSERs_%s_scale_%d.bed", chr_ids->at(i_chr), (int)scales_per_i_scale->at(i_scale));
				dump_BED(cur_scaled_filtered_minima_fp, p_value_end_pruned_minima);
			} // base_scale check.
		} // i_scale loop.

		// Free the memory for the per scale minima regions.
		for(int i_scale = 0; i_scale < (int)per_scale_minima_regions->size(); i_scale++)
		{
			if(per_scale_minima_regions->at(i_scale) != NULL)
			{
				delete_annot_regions(per_scale_minima_regions->at(i_scale));
			}
		} // i_scale loop.

		// Now merge the filtered/pruned minima regions.
		vector<t_annot_region*>* filtered_pruned_merged_minima_regions = merge_annot_regions(all_filtered_pruned_minima_regions, 1);

		// Do RD based and p-value based end pruning to the final set of peaks.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
		fprintf(stderr, "RD based end pruning.\n");

		vector<t_annot_region*>* rd_end_pruned_minima = prune_region_ends_per_window_average(filtered_pruned_merged_minima_regions, signal_profile, l_profile, 
																								control_profile, l_control,
																								rd_end_pruning_p_val, 1000*1000);

		//dump_BED("rd_end_pruned_minima.bed", rd_end_pruned_minima);
		delete_annot_regions(filtered_pruned_merged_minima_regions);

		// Do p-value minimization based end pruning.
		vector<t_annot_region*>* p_value_end_pruned_minima = rd_end_pruned_minima;
		if(do_post_peak_p_value_minimization)
		{
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
			fprintf(stderr, "P-value minimization based end pruning.\n");

			//p_value_end_pruned_minima = prune_region_ends_per_p_value_minimization(rd_end_pruned_minima, signal_profile, l_profile, control_profile, l_control, l_fragment);
			p_value_end_pruned_minima = prune_region_ends_per_modified_binomial_p_value_minimization(rd_end_pruned_minima, 
																									p_val_min_pruning_flip_prob,
																									signal_profile, l_profile, control_profile, l_control, l_fragment);
			delete_annot_regions(rd_end_pruned_minima);
		}

		// Replace the chromosome id's; necessary before strand signal based filtering.
		for(int i_reg = 0; i_reg < (int)p_value_end_pruned_minima->size(); i_reg++)
		{
			delete [] p_value_end_pruned_minima->at(i_reg)->chrom;
			p_value_end_pruned_minima->at(i_reg)->chrom = t_string::copy_me_str(chr_ids->at(i_chr));

			// Set the peak info for the current region.
			t_ER_info* cur_peak_info = new t_ER_info();
			cur_peak_info->fore_strand_mass_of_center = 0;
			cur_peak_info->rev_strand_mass_of_center = 0;
			cur_peak_info->total_fore_mass = 0;
			cur_peak_info->total_rev_mass = 0;
			p_value_end_pruned_minima->at(i_reg)->data = cur_peak_info;
		} // i_reg loop.

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
		fprintf(stderr, "Filtering peaks with respect to strand signal evenness.\n");

		// Set the per filter information.
		set_per_strand_info_per_peaks(p_value_end_pruned_minima,
				chr_ids->at(i_chr),
				chip_reads_dir,
				l_fragment);

		// Filter the per strand information.
		vector<t_annot_region*>* strand_filtered_peaks = new vector<t_annot_region*>();
		vector<t_annot_region*>* removed_peaks = new vector<t_annot_region*>();
		for(int i_reg = 0; i_reg < (int)p_value_end_pruned_minima->size(); i_reg++)
		{
			t_ER_info* cur_peak_info = (t_ER_info*)(p_value_end_pruned_minima->at(i_reg)->data);

			if(cur_peak_info->total_fore_mass > 0 &&
				cur_peak_info->total_rev_mass > 0 &&
				cur_peak_info->total_fore_mass / cur_peak_info->total_rev_mass > min_per_strand_evennes_fraction &&
				cur_peak_info->total_rev_mass / cur_peak_info->total_fore_mass > min_per_strand_evennes_fraction)
			{
				strand_filtered_peaks->push_back(p_value_end_pruned_minima->at(i_reg));
			}
			else
			{
				// Delete the peak info.
				delete cur_peak_info;
				removed_peaks->push_back(p_value_end_pruned_minima->at(i_reg));
			}
		} // i_reg loop.
		fprintf(stderr, "%s: %d peaks passed. (%d)\n", chr_ids->at(i_chr), (int)strand_filtered_peaks->size(), (int)p_value_end_pruned_minima->size());

		// Free unused peak memory.
		delete_annot_regions(removed_peaks);
		delete(p_value_end_pruned_minima);

		// Compute the p-values for each region.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
		fprintf(stderr, "Computing p-values.\n");

		for(int i_reg = 0; i_reg < (int)strand_filtered_peaks->size(); i_reg++)
		{
			t_ER_info* cur_peak_info = (t_ER_info*)(strand_filtered_peaks->at(i_reg)->data);

			// Get the p-value for the current region.
			double total_profile_signal = 0;
			double total_control_signal = 0;
			double cur_region_log_p_val = 0.0;
				
			if(strand_filtered_peaks->at(i_reg)->end < l_control && 
				strand_filtered_peaks->at(i_reg)->end < l_profile)
			{
				cur_peak_info->max_chip_mass = 0;
				cur_peak_info->climax_posn = 0;
				cur_peak_info->max_control_mass = 0;
				int cur_plateau_start = 0;
				//int cur_plateau_end = 0;
				//bool signal_increasing = false;
				bool in_plateau = false;

				// Following updates the ER information.
				for(int i = strand_filtered_peaks->at(i_reg)->start; i <= strand_filtered_peaks->at(i_reg)->end; i++)
				{
					total_profile_signal += floor(signal_profile[i]);
					total_control_signal += floor(control_profile[i]);

					int cur_signal_val = (int)signal_profile[i];
					int prev_signal_val = (int)signal_profile[i-1];

					if(!in_plateau &&
						cur_signal_val == prev_signal_val)
					{
						in_plateau = true;
						cur_plateau_start = i-1;
					}

					if(cur_signal_val > prev_signal_val)
					{
						cur_plateau_start = i;
					}

					// Did we start decreasing?
					if(cur_signal_val < prev_signal_val &&
						in_plateau)
					{
						// Is this a maximum within the peak? Then the climax position is the mid point of the plateau.
						if(cur_peak_info->max_chip_mass < floor(signal_profile[i-1]))
						{
							cur_peak_info->climax_posn = (int)(i-1 + cur_plateau_start) / 2;
							cur_peak_info->max_chip_mass = floor(signal_profile[i-1]);
						}

						cur_plateau_start = i;
						in_plateau = true;
					}
					
					/*if(cur_peak_info->max_chip_mass < floor(signal_profile[i]))
					{
						cur_peak_info->max_chip_mass = floor(signal_profile[i]);
					}*/

					if(cur_peak_info->max_control_mass < floor(control_profile[i]))
					{
						cur_peak_info->max_control_mass = floor(control_profile[i]);
					}
				} // i loop.

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{
				fprintf(stderr, "%s:%d-%d: %lf\n", chr_ids->at(i_chr), 
					strand_filtered_peaks->at(i_reg)->start, strand_filtered_peaks->at(i_reg)->end, cur_peak_info->climax_posn);
				getc(stdin);
}

				// Update the chip and control signal masses.
				cur_peak_info->total_chip_mass = total_profile_signal;
				cur_peak_info->total_control_mass = total_control_signal;

				// Compute the p-value with the option chosen.
				if(p_value_computation_type == P_VAL_PER_WIN_MEAN_SIGNAL)
				{
						cur_region_log_p_val = get_binomial_pvalue_per_region_neighborhood_window_control(signal_profile, control_profile, 
																			strand_filtered_peaks->at(i_reg)->start, strand_filtered_peaks->at(i_reg)->end, 
																			NULL,
																			l_fragment, 
																			true, 
																			l_p_val_norm_win,
																			l_profile,
																			l_control, 
																			1);
				}
				else if(p_value_computation_type == P_VAL_PER_WIN_MAX_SIGNAL)
				{
					cur_region_log_p_val = get_binomial_pvalue_per_region_neighborhood_window_control_per_max_profile_signal(signal_profile, control_profile, 
																		strand_filtered_peaks->at(i_reg)->start, strand_filtered_peaks->at(i_reg)->end, 
																		NULL,
																		l_fragment, 
																		true, 
																		l_p_val_norm_win,
																		l_profile,
																		l_control,
																		3);
				}
				else
				{
					fprintf(stderr, "Error: Did not understand the p-value computation type.\n");
					exit(0);
				}

				strand_filtered_peaks->at(i_reg)->significance_info = new t_significance_info();
				strand_filtered_peaks->at(i_reg)->significance_info->log_p_val = cur_region_log_p_val;
			} // profile length checks for ER start and end.
		} // i_reg loop.

		// Find the mappable troughs in the peaks.
		if(find_mappable_troughs)
		{
			// Load the mappability signal.
			char mapability_signal_profile_fp[1000];
			sprintf(mapability_signal_profile_fp, "%s/%s.bin", mapability_signal_dir, chr_ids->at(i_chr));

			if(check_file(mapability_signal_profile_fp))
			{
				// Load the mapability map signal profile, do filtering.
				int l_mapability_profile = 0;
				double* mapability_signal_profile = load_per_nucleotide_binary_profile(mapability_signal_profile_fp, l_mapability_profile);
				int mapability_scaling = l_read * 2;
				for(int i = 1; i <= l_mapability_profile; i++)
				{
					mapability_signal_profile[i] /= mapability_scaling;
				} // i loop.

				for(int i_reg = 0; i_reg < (int)strand_filtered_peaks->size(); i_reg++)
				{
					t_ER_info* cur_ER_info = (t_ER_info*)(strand_filtered_peaks->at(i_reg)->data);

					int cur_ER_trough_posn = get_trough_posn_per_ER(signal_profile, l_profile, 
																	mapability_signal_profile, l_mapability_profile, max_normalized_mapability_signal,
																	strand_filtered_peaks->at(i_reg)->start, strand_filtered_peaks->at(i_reg)->end, 
																	0);

					cur_ER_info->mappable_trough = cur_ER_trough_posn;
				} // i_reg loop.

				// Free mappability memory.
				delete [] mapability_signal_profile;
			}
			else
			{
				fprintf(stderr, "Could not find the multi-mappabibility signal, skipping mappable trough detection.\n");
			}
		} // mappable trough identification check.

		// Insert the pruned peaks to the end of peaks list.
		all_peaks->insert(all_peaks->end(), strand_filtered_peaks->begin(), strand_filtered_peaks->end());

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{
		// Dump the filtered peaks.
		char cur_chr_regions_fp[1000];
		sprintf(cur_chr_regions_fp, "scored_peaks_%s_%.1f_%.1f_%.2f.bed", chr_ids->at(i_chr), base_scale_l_win, end_scale_l_win, log_step);
		FILE* f_cur_chr_regions = open_f(cur_chr_regions_fp, "w");
		for(int i_reg = 0; i_reg < (int)strand_filtered_peaks->size(); i_reg++)
		{
			t_ER_info* cur_ER_info = (t_ER_info*)(strand_filtered_peaks->at(i_reg)->data);
			fprintf(f_cur_chr_regions, "%s\t%d\t%d\t.\t%lf\t+\t%d\t%d\t%d\t%d\t%lf\t%lf\n", chr_ids->at(i_chr), 
				strand_filtered_peaks->at(i_reg)->start, strand_filtered_peaks->at(i_reg)->end, strand_filtered_peaks->at(i_reg)->significance_info->log_p_val,
				(int)((cur_ER_info->fore_strand_mass_of_center + cur_ER_info->rev_strand_mass_of_center) / 2),
				(int)((cur_ER_info->fore_climax_posn + cur_ER_info->rev_climax_posn) / 2), 
				(int)(cur_ER_info->climax_posn), 
				(int)(cur_ER_info->mappable_trough), 
				cur_ER_info->total_chip_mass / cur_ER_info->total_control_mass,
				cur_ER_info->max_chip_mass / cur_ER_info->max_control_mass);
		} // i_reg loop.
		fclose(f_cur_chr_regions);
}

		// Free signal memory.
		delete [] signal_profile;
		delete [] control_profile;
	} // i_chr loop.

	// Do BH corrections, sort with respect to q-values.
if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
	fprintf(stderr, "Doing BH correction on %d peaks.\n", (int)all_peaks->size());

	get_benjamini_hochberg_corrected_p_values(all_peaks);

	// Filter the peaks wrt q-value.
	//double q_val_cutoff = log(0.05);
	double log_q_val_cutoff = log(q_value_cutoff);
	vector<t_annot_region*>* significant_peaks = new vector<t_annot_region*>();
	for(int i_reg = 0; i_reg < (int)all_peaks->size(); i_reg++)
	{
		// check if the q-value is below the cutoff.
		if(all_peaks->at(i_reg)->significance_info->log_q_val < log_q_val_cutoff)
		{
			significant_peaks->push_back(all_peaks->at(i_reg));
		}
	} // i_reg loop.

	// Dump the final set of peaks.
	char op_bed_fp[1000];
	sprintf(op_bed_fp, "ERs_%.1f_%.1f_%.2f_%d_%.1f.bed", base_scale_l_win, end_scale_l_win, log_step, l_p_val_norm_win, filtered_vs_non_filtered_max_scaler);
	dump_ERs_per_q_value_and_summit(significant_peaks, op_bed_fp);
} // get_peaks function.

/*
Find the maxima, sort them, find the smallest positions in between the top maxima regions.
*/
int get_trough_posn_per_ER(double* signal_profile, int l_profile, 
							double* multi_mapp_signal, int l_multi_map_signal, double max_multi_mapp_val,
							int ER_start, int ER_end, 
							int l_trough_win)
{
	//vector<int>* maxima_posns = new vector<int>();
	//double min_max_val = 0.0;

	vector<t_extrema_node*>* maxima_nodes = new vector<t_extrema_node*>();
	vector<t_extrema_node*>* minima_nodes = new vector<t_extrema_node*>();
	int* derivative_map = new int[ER_end - ER_start + 2];
	memset(derivative_map, 0, sizeof(int) * (ER_end - ER_start + 2));
	get_extrema_per_plateaus(&signal_profile[ER_start], ER_end-ER_start+1, maxima_nodes, minima_nodes, derivative_map, 0);

	if((int)maxima_nodes->size() == 0 ||
		(int)minima_nodes->size() == 0)
	{
		return(0);
	}

	sort(maxima_nodes->begin(), maxima_nodes->end(), sort_extremas_per_decreasing_height);

	// Sort the minima with increasing height then take the top minima.
	sort(minima_nodes->begin(), minima_nodes->end(), sort_extremas_per_height);

	// Start from the top minimum and check if it is around one of the top maxima.
	vector<int>* rank_sum_per_minimum = new vector<int>();
	for(int i_min = 0; i_min < (int)minima_nodes->size(); i_min++)
	{
		t_extrema_node* prev_max = NULL;
		t_extrema_node* next_max = NULL;
		int prev_max_rank = 0;
		int next_max_rank = 0;

		for(int i_max = 0; i_max < (int)maxima_nodes->size(); i_max++)
		{
			if(prev_max == NULL &&
				maxima_nodes->at(i_max)->extrema_posn < minima_nodes->at(i_min)->extrema_posn)
			{
				prev_max_rank = i_max;
				prev_max = maxima_nodes->at(i_max);
			}

			if(next_max == NULL &&
				maxima_nodes->at(i_max)->extrema_posn > minima_nodes->at(i_min)->extrema_posn)
			{
				next_max_rank = i_max;
				next_max = maxima_nodes->at(i_max);
			}
		} // i_max loop.

		if(prev_max == NULL ||
			next_max == NULL)
		{
			rank_sum_per_minimum->push_back((int)maxima_nodes->size() + (int)maxima_nodes->size() + (int)minima_nodes->size());
			continue;
		}

		double total_multi_mapp_signal = 0;
		for(int i = prev_max->extrema_posn; i <= next_max->extrema_posn; i++)
		{
			total_multi_mapp_signal += multi_mapp_signal[i];
		} // i loop.

		// Make sure that the region in this trough is mappable.
		if(total_multi_mapp_signal / (next_max->extrema_posn - prev_max->extrema_posn+1) < max_multi_mapp_val)
		{
			rank_sum_per_minimum->push_back(i_min + prev_max_rank + next_max_rank);
		}
		else
		{
			// Set the maximum possible rank to this minimum.
			rank_sum_per_minimum->push_back((int)maxima_nodes->size() + (int)maxima_nodes->size() + (int)minima_nodes->size());
		}
	} // i_min loop.

	// Find the minimum with smallest rank-sum statistic.
	int smallest_rank_sum_stat = 1000;
	int smallest_rank_sum_min_i = 0;
	for(int i_min = 0; i_min < (int)rank_sum_per_minimum->size(); i_min++)
	{
		if(smallest_rank_sum_stat > rank_sum_per_minimum->at(i_min))
		{
			smallest_rank_sum_stat = rank_sum_per_minimum->at(i_min);
			smallest_rank_sum_min_i = i_min;
		}
	} // i_min loop.

	int trough_posn = minima_nodes->at(smallest_rank_sum_min_i)->extrema_posn;

	delete_extrema_nodes(minima_nodes);
	delete_extrema_nodes(maxima_nodes);

	// Free the current decomposition.
	delete [] derivative_map;

	return(trough_posn);
}

void dump_ERs_per_q_value_and_summit(vector<t_annot_region*>* regions, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "w");

	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		t_ER_info* cur_reg_info = (t_ER_info*)(regions->at(i_reg)->data);

		double FC = cur_reg_info->total_chip_mass / cur_reg_info->total_control_mass;

		//t_significance_info* sig_info = (t_significance_info*)(regions->at(i_reg)->significance_info);
		//fprintf(f_op, "%s\t%d\t%d\t.\t%lf\t%c\t%d\t%d\t%d\t%d\t%lf\t%lf\n", regions->at(i_reg)->chrom, regions->at(i_reg)->start, regions->at(i_reg)->end, 
		//	sig_info->log_q_val, 
		//	regions->at(i_reg)->strand,
		//	(int)((cur_reg_info->rev_strand_mass_of_center + cur_reg_info->fore_strand_mass_of_center) / 2), 
		//	(int)((cur_reg_info->fore_climax_posn + cur_reg_info->rev_climax_posn) / 2), 
		//	(int)(cur_reg_info->climax_posn), 
		//	(int)(cur_reg_info->mappable_trough), 
		//	FC,
		//	cur_reg_info->max_chip_mass / cur_reg_info->max_control_mass);
		t_significance_info* sig_info = (t_significance_info*)(regions->at(i_reg)->significance_info);
		fprintf(f_op, "%s\t%d\t%d\t.\t%lf\t%c\t%d\t%d\t%lf\n", regions->at(i_reg)->chrom, regions->at(i_reg)->start, regions->at(i_reg)->end, 
			sig_info->log_q_val, 
			regions->at(i_reg)->strand,
			(int)(cur_reg_info->climax_posn), 
			(int)(cur_reg_info->mappable_trough), 
			FC);
	} // i_reg loop.

	fclose(f_op);
}

vector<t_annot_region*>* filter_strand_uneven_peaks(vector<t_annot_region*>* peak_regions,
		char* chrom,
		char* chip_reads_dir,
		int l_fragment,
		double min_forward_reverse_strand_signal_fraction)
{
	vector<char*>* chr_ids = NULL;

	if(chrom == NULL)
	{
		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", chip_reads_dir);
		chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load the chromosome id's for the chip data.\n");
			exit(0);
		}
	}
	else
	{
		chr_ids = new vector<char*>();
		chr_ids->push_back(t_string::copy_me_str(chrom));
	}

	// Process each chromosome.
	vector<t_annot_region*>* filtered_peaks = new vector<t_annot_region*>();
	//vector<t_annot_region*>* removed_peaks = new vector<t_annot_region*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "%s..", chr_ids->at(i_chr));

		// Skip chromosome M for now.
		if(strcmp(chr_ids->at(i_chr), "M") == 0)
		{
			continue;
		}

		// Generate the current signal profile.
		int l_buffer = 300*1000*1000;
		int l_profile = 0;
		double* forward_signal = new double[l_buffer + 2];	
		double* reverse_signal = new double[l_buffer + 2];	
		char cur_chr_chip_reads_fp[1000];
		sprintf(cur_chr_chip_reads_fp, "%s/%s_mapped_reads.txt", chip_reads_dir, chr_ids->at(i_chr));
		buffer_per_nucleotide_profile_no_buffer(cur_chr_chip_reads_fp, l_fragment, 
			NULL, forward_signal, reverse_signal,
			l_buffer, l_profile);

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{
		char fore_strand_signal_fp[1000];
		sprintf(fore_strand_signal_fp, "%s_fore_strand.bin", chr_ids->at(i_chr));
		dump_per_nucleotide_binary_profile(forward_signal, l_profile, fore_strand_signal_fp);
		char rev_strand_signal_fp[1000];
		sprintf(rev_strand_signal_fp, "%s_rev_strand.bin", chr_ids->at(i_chr));
		dump_per_nucleotide_binary_profile(reverse_signal, l_profile, rev_strand_signal_fp);
}
			
		vector<t_annot_region*>* cur_chr_regs = get_regions_per_chromosome(peak_regions, chr_ids->at(i_chr));
		fprintf(stderr, "Filtering %d regions on chromosome %s\n", (int)cur_chr_regs->size(), chr_ids->at(i_chr));
			
		int n_passing_peaks = 0;
		for(int i_reg = 0; i_reg < (int)cur_chr_regs->size(); i_reg++)
		{
			double total_forward_signal = 0.0;
			double total_reverse_signal = 0.0;
			for(int i = cur_chr_regs->at(i_reg)->start; i <= cur_chr_regs->at(i_reg)->end; i++)
			{
				if(i < l_profile)
				{
					total_forward_signal += forward_signal[i];
					total_reverse_signal += reverse_signal[i];
				}
			} // i loop.

			// Both strand signals much be larger than 0.
			if(total_forward_signal > 0 &&
				total_reverse_signal > 0 &&
				total_reverse_signal / total_forward_signal > min_forward_reverse_strand_signal_fraction &&
				total_forward_signal / total_reverse_signal > min_forward_reverse_strand_signal_fraction)
			{
				t_annot_region* new_region = duplicate_region(cur_chr_regs->at(i_reg));
				filtered_peaks->push_back(new_region);
				n_passing_peaks++;
			}
			else
			{
				//removed_peaks->push_back(cur_chr_regs->at(i_reg));
			}
		} // i_reg loop.

		fprintf(stderr, "%s: %d peaks passed. (%d)\n", chr_ids->at(i_chr), n_passing_peaks, (int)cur_chr_regs->size());

		// For each peak, compute the fraction of read signal that is in mapable regions.
		delete [] forward_signal;
		delete [] reverse_signal;
		delete cur_chr_regs;
	} // chromosome loop.

	//dump_BED("filtered_peaks.bed", filtered_peaks);
	//dump_BED("removed_peaks.bed", removed_peaks);
	//delete removed_peaks;

	t_string::clean_string_list(chr_ids);

	return(filtered_peaks);
}

vector<t_annot_region*>* set_per_strand_info_per_peaks(vector<t_annot_region*>* peak_regions,
		char* chrom,
		char* chip_reads_dir,
		int l_fragment)
{
	vector<char*>* chr_ids = NULL;

	if(chrom == NULL)
	{
		char chr_ids_fp[1000];
		sprintf(chr_ids_fp, "%s/chr_ids.txt", chip_reads_dir);
		chr_ids = buffer_file(chr_ids_fp);
		if(chr_ids == NULL)
		{
			fprintf(stderr, "Could not load the chromosome id's for the chip data.\n");
			exit(0);
		}
	}
	else
	{
		chr_ids = new vector<char*>();
		chr_ids->push_back(t_string::copy_me_str(chrom));
	}

	// Process each chromosome.
	vector<t_annot_region*>* filtered_peaks = new vector<t_annot_region*>();
	//vector<t_annot_region*>* removed_peaks = new vector<t_annot_region*>();
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "%s..", chr_ids->at(i_chr));

		// Skip chromosome M for now.
		if(strcmp(chr_ids->at(i_chr), "M") == 0)
		{
			continue;
		}

		// Generate the current signal profile.
		int l_buffer = 300*1000*1000;
		int l_profile = 0;
		double* forward_signal = new double[l_buffer + 2];	
		double* reverse_signal = new double[l_buffer + 2];	
		char cur_chr_chip_reads_fp[1000];
		sprintf(cur_chr_chip_reads_fp, "%s/%s_mapped_reads.txt", chip_reads_dir, chr_ids->at(i_chr));
		buffer_per_nucleotide_profile_no_buffer(cur_chr_chip_reads_fp, l_fragment, 
			NULL, forward_signal, reverse_signal,
			l_buffer, l_profile);

if(__DUMP_PEAK_CALLING_UTILS_MSGS__)
{
		char fore_strand_signal_fp[1000];
		sprintf(fore_strand_signal_fp, "%s_fore_strand.bin", chr_ids->at(i_chr));
		dump_per_nucleotide_binary_profile(forward_signal, l_profile, fore_strand_signal_fp);
		char rev_strand_signal_fp[1000];
		sprintf(rev_strand_signal_fp, "%s_rev_strand.bin", chr_ids->at(i_chr));
		dump_per_nucleotide_binary_profile(reverse_signal, l_profile, rev_strand_signal_fp);
}
			
		vector<t_annot_region*>* cur_chr_regs = get_regions_per_chromosome(peak_regions, chr_ids->at(i_chr));
		fprintf(stderr, "Filtering %d regions on chromosome %s\n", (int)cur_chr_regs->size(), chr_ids->at(i_chr));
			
		//int n_passing_peaks = 0;
		for(int i_reg = 0; i_reg < (int)cur_chr_regs->size(); i_reg++)
		{
			t_ER_info* cur_peak_info = (t_ER_info*)(cur_chr_regs->at(i_reg)->data);

			double total_forward_signal = 0.0;
			double total_reverse_signal = 0.0;
			double total_fore_moment = 0;
			double total_rev_moment = 0;
			cur_peak_info->fore_max_val = 0;
			cur_peak_info->fore_climax_posn = 0;
			cur_peak_info->rev_max_val = 0;
			cur_peak_info->rev_climax_posn = 0;

			// Process forward strand info.
			for(int i = cur_chr_regs->at(i_reg)->start; i <= cur_chr_regs->at(i_reg)->end; i++)
			{
				if(i < l_profile)
				{
					total_forward_signal += forward_signal[i];
					total_fore_moment += (i - cur_chr_regs->at(i_reg)->start + 1) * forward_signal[i];

					if(cur_peak_info->fore_max_val < forward_signal[i])
					{
						cur_peak_info->fore_max_val = forward_signal[i];
						cur_peak_info->fore_climax_posn = i;
					}
				}
			} // i loop.

			// Process reverse strand info.
			for(int i = cur_chr_regs->at(i_reg)->end; i >= cur_chr_regs->at(i_reg)->start; i--)
			{
				if(i < l_profile)
				{
					total_reverse_signal += reverse_signal[i];
					total_rev_moment += (i - cur_chr_regs->at(i_reg)->start + 1) * reverse_signal[i];

					if(cur_peak_info->rev_max_val < reverse_signal[i])
					{
						cur_peak_info->rev_max_val = reverse_signal[i];
						cur_peak_info->rev_climax_posn = i;
					}
				}
			} // i loop.

			cur_peak_info->total_fore_mass = total_forward_signal;
			cur_peak_info->total_rev_mass = total_reverse_signal;

			cur_peak_info->fore_strand_mass_of_center = (int)(total_fore_moment / total_forward_signal) + cur_chr_regs->at(i_reg)->start;
			cur_peak_info->rev_strand_mass_of_center = (int)(total_rev_moment / total_reverse_signal) + cur_chr_regs->at(i_reg)->start;
		} // i_reg loop.

		//fprintf(stderr, "%s: %d peaks passed. (%d)\n", chr_ids->at(i_chr), n_passing_peaks, (int)cur_chr_regs->size());

		// For each peak, compute the fraction of read signal that is in mapable regions.
		delete [] forward_signal;
		delete [] reverse_signal;
		delete cur_chr_regs;
	} // chromosome loop.

	//dump_BED("filtered_peaks.bed", filtered_peaks);
	//dump_BED("removed_peaks.bed", removed_peaks);
	//delete removed_peaks;

	t_string::clean_string_list(chr_ids);

	return(filtered_peaks);
}