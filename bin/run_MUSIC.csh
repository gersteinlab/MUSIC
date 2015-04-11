#!/bin/bash 

# Check out the readme at the end.

# Check the number of arguments.
if [ $# -lt 3 ]
then 
	echo -e "USAGE: run_MUSIC.csh [Options] [Arguments]
-preprocess [Reads file path] [Output dir]
-remove_duplicates [Preprocessed reads dir] [Sorted reads dir] [Pruned reads dir]
-get_[relaxed/optimal]_[punctate/broad]_ERs [chip preprocessed dir] [input preprocessed dir] [Multi-mappability profile directory]\n";

	exit;
fi

if [ "$1" == "-preprocess" ]
then
	read_file_name=`basename $2`;
	output_dir=$3;
	file_extension=${read_file_name#*.};

	if [[ "$read_file_name" == *.tagAlign.gz ]]
	then
		gzip -cd $2 | MUSIC -preprocess tagAlign stdin $output_dir
	elif [[ "$read_file_name" == *.tagAlign ]]
	then
		cat $2 | MUSIC -preprocess tagAlign stdin $output_dir
	elif [[ "$read_file_name" == *.bed.gz ]]
	then
		gzip -cd $2 | MUSIC -preprocess BED6 stdin $output_dir
	elif [[ "$read_file_name" == *.bed ]]
	then
		cat $2 | MUSIC -preprocess BED stdin $output_dir
	elif [[ "$read_file_name" == *.bowtie ]]
	then
		cat $2 | MUSIC -preprocess bowtie stdin $output_dir
	elif [[ "$read_file_name" == *.bowtie.gz ]]
	then
		gzip -cd $2 | MUSIC -preprocess bowtie stdin $output_dir
	elif [[ "$read_file_name" == *.bam ]]
	then
		samtools view $2 | MUSIC -preprocess SAM stdin $output_dir
	fi
elif [ "$1" == "-remove_duplicates" ]
then
	preprocessed_reads_dir=$2;
	sorted_reads_dir=$3;
	pruned_reads_dir=$4;
	MUSIC -sort_reads $preprocessed_reads_dir $sorted_reads_dir
	MUSIC -remove_duplicates $sorted_reads_dir 2 $pruned_reads_dir
elif [ "$1" == "-get_optimal_broad_ERs" ]
then
	chip_processed_dir=$2;
	input_processed_dir=$3;
	mappability_dir=$4;

	# Clean directory.
	rm -f *.bed

	# Run MUSIC: -l_p selection.
	rm -f l_p_param_stats.txt
	MUSIC -get_per_win_p_vals_vs_FC -chip ${chip_processed_dir} -control ${input_processed_dir} -l_win_step 50 -l_win_min 200 -l_win_max 5000

	# Finally, get peaks; select l_p from the file.
	MUSIC -get_multiscale_broad_ERs -chip ${chip_processed_dir} -control ${input_processed_dir} -l_p 0 -mapp ${mappability_dir} -l_mapp 50 -q_val 0.05
elif [ "$1" == "-get_relaxed_broad_ERs" ]
then
	chip_processed_dir=$2;
	input_processed_dir=$3;
	mappability_dir=$4;

	# Clean directory.
	rm -f *.bed

	# Run MUSIC: -l_p selection.
	rm -f l_p_param_stats.txt
	MUSIC -get_per_win_p_vals_vs_FC -chip ${chip_processed_dir} -control ${input_processed_dir} -l_win_step 50 -l_win_min 200 -l_win_max 5000

	# Finally, get peaks; select l_p from the file.
	MUSIC -get_multiscale_broad_ERs -chip ${chip_processed_dir} -control ${input_processed_dir} -l_p 0 -mapp ${mappability_dir} -l_mapp 50 -q_val 0.1
elif [ "$1" == "-get_relaxed_punctate_ERs" ]
then
	chip_processed_dir=$2;
	input_processed_dir=$3;
	mappability_dir=$4;

	# Clean directory.
	rm -f *.bed

	# Do peak calling directly, no parameter selection for calling punctate ERs.
	MUSIC -get_multiscale_punctate_ERs -chip ${chip_processed_dir} -control ${input_processed_dir} -l_p 1500 -mapp ${mappability_dir} -l_mapp 50 -q_val 0.1
elif [ "$1" == "-get_optimal_punctate_ERs" ]
then
	chip_processed_dir=$2;
	input_processed_dir=$3;
	mappability_dir=$4;

	# Clean directory.
	rm -f *.bed

	# Do peak calling directly, no parameter selection for calling punctate ERs.
	MUSIC -get_multiscale_punctate_ERs -chip ${chip_processed_dir} -control ${input_processed_dir} -l_p 1500 -mapp ${mappability_dir} -l_mapp 50 -q_val 0.05
fi

# Done.
exit;

#
# This script can be used to run MUSIC:
# Punctate ER calling options: -get_optimal_punctate_ERs (Optimal ERs), -get_relaxed_punctate_ERs (Relaxed ER calls)
# Use with: H3K4me3, H3K27ac, H3K4me1, H3K4me3, H3K9ac, H2az.  
# 
# Broad ER calling options: -get_optimal_broad_ERs (Optimal ERs), -get_relaxed_broad_ERs (Relaxed ER calls)
# Use With: H3K9me3, H3K36me3, H3K27me3, H3K79me2, H4K20me1. 
#
# Following commands give an example run for replicates and pooled run. 
# Note that the input (control) needs to be preprocessed and duplicate removed before ER calling. This is used in ER calling.
#
# In order to call punctate ERs (see above for which marks are considered punctate), just change "broad" to "punctate" in ER calling command line.
# 
# Multi-mappability profiles for the read length is very useful for running especially the broad peaks. 
# You can find some of these at http://archive.gersteinlab.org/proj/MUSIC/multimap_profiles/ 
# Email to me (arif.harmanci@yale.edu) for generating a new one.
#

# The input (control) should be preprocessed first.
input_fp="../../wgEncodeBroadHistoneGm12878ControlStdAlnRep1.bam";
mkdir input;mkdir input/preprocessed input/sorted input/pruned
run_MUSIC.csh -preprocess ${input_fp} input/preprocessed 
run_MUSIC.csh -remove_duplicates input/preprocessed input/sorted input/pruned

# Rep run.
rep_fp="../../wgEncodeBroadHistoneGm12878H3k27me3StdAlnRep1.tagAlign.gz";
mappability_map="../../../mappability/36bp";
input_processed_dir="../../input/pruned";
rm -f -r preprocessed sorted pruned
mkdir preprocessed sorted pruned
run_MUSIC.csh -preprocess ${rep_fp} preprocessed
run_MUSIC.csh -remove_duplicates preprocessed sorted pruned
run_MUSIC.csh -get_optimal_broad_ERs pruned ${input_processed_dir} ${mappability_map}

# Pooled processing
rep1_fp="../../wgEncodeBroadHistoneGm12878H3k36me3StdAlnRep1.tagAlign.gz";
rep2_fp="../../wgEncodeBroadHistoneGm12878H3k36me3StdAlnRep2.tagAlign.gz";
mappability_map="../../../mappability/36bp";
input_processed_dir="../../input/pruned";
rm -f -r preprocessed sorted pruned
mkdir preprocessed sorted pruned
run_MUSIC.csh -preprocess ${rep1_fp} preprocessed
run_MUSIC.csh -preprocess ${rep2_fp} preprocessed
run_MUSIC.csh -remove_duplicates preprocessed sorted pruned
run_MUSIC.csh -get_optimal_broad_ERs pruned ${input_processed_dir} ${mappability_map}