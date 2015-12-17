<html>
<font face="arial">
<!------------------- <title>MUSIC</title> ---->
<div style="text-align:center"><img src="music_logo.png" alt="Could not load logo." width="1000" align="center"></div>
<!-------------------<h1>MUSIC: MUltiScale enrIchment Calling</h1>---->
<!---------- MUSIC is an algorithm for identification of enriched regions at multiple scales.<br>
It takes as input:(1) Mapped ChIP and control reads, (2) Smoothing scale window lengths (in base pairs), <br>
and outputs: (1) Enriched regions at multiple scales, (2) Significantly enriched regions from all the scales.<br>
<br><br><br><br><br>---------->
<title>MUSIC</title>
<br><br>
MUSIC is an algorithm for identification of enriched regions at multiple scales in the read depth signals from ChIP-Seq experiments.<br>It takes as input: <br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- Mapped ChIP and control reads (optional),<br>
- Smoothing scale window lengths (in base pairs),
</div><br>
and outputs:<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- Enriched regions at multiple scales,<br>
- Significantly enriched regions from all the scales.
</div>
<br>
Unlike other ER identification methods, MUSIC allows analyzing the scale length spectrum of ChIP-Seq datasets and also selecting a user specific slice in the 
scale length spectrum with custom granularity and generates the enriched regions at each length scale.

MUSIC does not strictly require a control experiment (for example mock IP, input DNA) to be performed. It is, however, strongly advised to generate control datasets matching 
the sample of interest (See Landt et al 2012).

<h2>Download and Installation</h2>
You can download MUSIC C++ code <a href="https://github.com/gersteinlab/MUSIC/archive/master.zip">here</a>. There are no dependencies for building MUSIC. After download, type:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
unzip MUSIC.zip<br>
cd MUSIC<br>
make clean<br>
make
</font></i>
</div><br>
to build MUSIC. The executable is located under directory <font face="courier">bin/</font>. It may be useful to install <a href="http://samtools.sourceforge.net/">samtools</a> for processing BAM files.

To get help on which options are available, use:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
MUSIC -help
</font>
</div>

<h2>Usage</h2>
MUSIC run starts with preprocessing the reads for ChIP and control samples (Note that we use samtools for converting BAM file to SAM files.):
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
mkdir chip;mkdir input<br>
samtools view chip.bam | MUSIC -preprocess SAM stdin chip/ <br>
samtools view input.bam | MUSIC -preprocess SAM stdin input/ <br>
</font>
</div><br>
If there are multiple replicates to be pooled, they can be done at once or separately. If done separately, MUSIC pools the reads automatically. Then it is 
necessary to sort and remove duplicate reads in control and ChIP samples:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
mkdir chip/sorted;mkdir chip/dedup;mkdir input/sorted;mkdir input/dedup<br>
MUSIC -sort_reads chip chip/sorted <br>
MUSIC -sort_reads input input/sorted <br>
MUSIC -remove_duplicates chip/sorted 2 chip/dedup <br>
MUSIC -remove_duplicates input/sorted 2 input/dedup <br>
</font>
</div><br>
We do enriched region identification:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
MUSIC -get_multiscale_broad_ERs -chip chip/dedup -control input/dedup -mapp Mappability_36bp -l_mapp 36 -begin_l 1000 -end_l 16000 -step 1.5
</font>
</div><br>
This code tells MUSIC to identify the enriched regions starting from 1kb smoothing window length upto 16kb with multiplicative factor of 1.5 using the default
parameters for the remaining parameters. The ERs for each scale are dumped. 

In case there is no control, skip the option that specifies preprocessed control reads directory (i.e., '-control input/dedup')

There are 3 different ER identification modes by default. 
<h3>-get_TF_peaks</h3>
Identifies point binding events. This uses very small scale level to identify the small transcription factor binding peaks. Use this mode for TF's like CTCF. MUSIC aims
at trimming the reported regions and identifying the peaks with most dense signal in it.

<h3>-get_multiscale_punctate_ERs</h3>
Identifies punctate enriched regions. Uses 100 base pairs to 2 kbps scale levels. This mode is useful for punctate histone marks like H3K4me3, H3K27ac, and marks that behave
in a mixed manner with a dominating spectrum at punctate scales like H3K4me1.

<h3>-get_multiscale_broad_ERs</h3>
Identifies the ERs at broad scales that can span megabases. Use this option for marks like H3K9me3, H3K27me3, H3K36me3, H3K79me2, H4K20me1... Make sure you select p-value 
normalization parameter to balance power and false positive rate.

<br><br>
MUSIC can save the smoothed tracks in <a href="http://genome.ucsc.edu/goldenPath/help/bedgraph.html">bedGraph</a> format that can be viewed locally:
</div><br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
MUSIC -write_MS_decomposition -chip chip/dedup -control input/dedup -mapp Mappability_36bp -l_mapp 36 
</font>
</div><br>
The smoothed bedGraph files are usually very small in size and can easily be stored/transferred.

<h2>Output format</h2>
MUSIC output a large number of files that contain the SSERs at different scales (named SSER_....bed). <br>

The final set of ERs are reported in two files: One is broadPeak formatted (http://genome.ucsc.edu/FAQ/FAQformat.html#format13)<br>

Other file is in an extended BED format and has 9 columns:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
[chromosome]	[start]	[end]	["."]	[log_10 Q-value]	[Strand ("+")]	[Summit position]	[Mappable Trough Position]	[Fold Change]
</font>
</div><br>
The entries are sorted with respect to increasing Q-values.

Note that MUSIC reports log_10(Q-values) and not -log_10(Q-value).

<h2>Parameter Selection Guideline for the Studied ChIP-Seq Datasets</h2>
MUSIC has a set of default parameter sets for broad marks (like H3K36me3, H3K27me3), punctate marks (H3K4me3), and point binding (like transcription factors) which work well in most 
cases. When one is not sure about the ER scale spectrum for a signal profile at hand, one can perform a scale spectrum analysis using a large spectrum with dense sampling and get an 
idea about the dominant ER length scale (if there is one) for the signal profile and match the parameters used in the manuscript.

<h2>Parametrization for New ChIP-Seq Datasets</h2>
When a new ChIP-Seq dataset is going to be processed, it is necessary to choose begin and end length scales and p-value normalization window length. The length
scales enables one to concentrate on the correct scale spectrum and p-value normalization window length compensates for variation in sequencing depth and allows
controlling the estimated false positive rates. <br>

There are two steps to parameter selection:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
1. Select a stringent p-value normalization window length: -get_per_win_p_vals_FC option. 
</font>
</div><br>
This option estimates the false positive and negative rates using a large selection of p-value window lengths. The output is a file where each row look like this:
<br>
...<br>
l_win: 1700	FNR: (FC:0.001) (p-val:0.001)	FPR: (FC:0.010) (p-val:0.005)	Sentitivity: 0.999<br>
...<br>
<br>
This option evaluates several windows lengths and estimates the false positive rate and false negative rates. We recommends using the maximum window length (l_win) where false
positive rate (FPR) for FC and p-val values are below 1%.<br>
<br>
Note that you can skip going through the file manually if you specify p-value normalization window length as 0 ('-l_p 0') in the peak calling step; which tells MUSIC to select p-value normalization window length 
from the above file automatically. For this to work, make sure that you do not delete any files after running -get_per_win_p_vals_FC option, otherwise MUSIC will complain that it cannot find 
the file. See below for complete automation (of p-value normalization window length selection) with default parameters. <br>
<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
2. Using the p-value normalization window length in step 1, generate the scale specific ER scale spectrum: -get_scale_spectrum option.<br>

</font><br>
</div><br>
MUSIC generates the spectrum (using the scale lengths 100 base pairs to 1 megabase. The output is reported in a text file where each 
row corresponds to a scale length. In each row, the coverage of the SSERs is given. For example: 
<br>
...<br>
27	56815.13	111	11675769<br>
...<br>
<br>
where 2nd column is the scale length and 4th column is the coverage of the ERs that are specific to that scale. It is best to plot the spectrum, i.e., the scale lengths versus 
the fraction of coverage of the SSERs (4th column in the file), then match the spectrum with the studied HMs in the 
manuscript. If the spectrum is very different from all of the parametrized ChIP-Seq datasets, it is useful to generate the statistics on ER length distribution 
and distribution of ER-ER distances and use the procedure described in the manuscript to select the scale levels.<br>
<br>
The other parameters (namely, gamma and sigma) do not depend on experimental variables and are optimized for minimizing overmerging and maximizing sensitivity.
We suggest using the values specified in the manuscript, i.e., gamma=4, sigma=1.5.<br>

<h2>Important Note on Punctate ERs</h2>
After this analysis, if the scale spectrum turns out to be punctate, i.e., the dominating scale is smaller than 10kb; there is a possibility that the p-value normalization window
length parameter (l_p) yields low sensitivity. To compensate for this, we recommend running with default l_p parameter of MUSIC.

<h2>Running MUSIC with Default Parameters and Automatic Selection of l_p Parameter: </h2>
We just added a new script run_MUSIC.csh. This script automates the parameter selection for -l_p option. This script automates running MUSIC with default parameters. It simply calls MUSIC 
with above parameters in order. Here is how this script can be used: 
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
read_fp="wgEncodeBroadHistoneGm12878H3k27me3StdAlnRep1.bam";<br>
mappability_map="mappability/36bp"; <br>
input_processed_dir="input/pruned";<br>
mkdir preprocessed sorted pruned<br>
run_MUSIC.csh -preprocess ${read_fp} preprocessed<br>
run_MUSIC.csh -remove_duplicates preprocessed sorted pruned<br>
run_MUSIC.csh -get_optimal_broad_ERs pruned ${input_processed_dir} ${mappability_map}<br>
</font>
</div><br>
Make sure that run_MUSIC.csh is included in PATH. Currently only BAM files are supported in preprocessing. We will add more file types, soon. Note that it is necessary to install samtools for making 
sure this runs smoothly. The above code calls the ER's using the optimal l_p selection with the default parameters.

<h2>Multi-Mappability Signals</h2>
Using Mappability correction increases the accuracy of MUSIC. You can download the multi-mappability signals for several common read lengths <a href="http://archive.gersteinlab.org/proj/MUSIC/multimap_profiles/">here</a>.

<h2>Multi-Mappability Profile Generation</h2>
To generate the multi-mappability profile, MUSIC depends on a short read aligner. By default, <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">bowtie2</a> generated SAM alignments are supported. Following are necessary for multi-Mappability profile generation:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
- The FASTA file for the genomic sequence with all the chromosomes<br>
- The read length<br>
- bowtie2 installation and the genome indices for bowtie2<br>
- Read length<br>
</font>
</div><br>
Use this script generate_multimappability_signal.csh under bin/ directory to generate the multi-mappability profile. For example, to generate the 50 bp multi-Mappability profile for human genome, 
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
cd bin<br>
./bin/generate_multimappability_signal.csh hg19.fa 50 hg19_indexes/hg19
</font>
</div><br>
Note that the bowtie indices must be supplied in the command line. For the above command, they are under hg19_indexes/. This processes the FASTA file and 
writes two temporary scripts (temp_map_reads.csh, temp_process_mapping.csh) that fragments the genome, maps the fragments, and then generates the profile.
Next, it is necessary to run these two scripts, in order:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
./temp_map_reads.csh<br>
./temp_process_mapping.csh
</font>
</div><br>
After the scripts are run, multi-mappability profile for each chromosome should be created with '.bin' extension. 'temp_map_reads.csh' is the most time consuming script that maps the fragments to the genome. 
Each line in the script is a command that maps the reads to a chromosome that can be run in parallel on a cluster. It is important to make sure 'temp_map_reads.csh' 
finishes before running 'temp_process_mapping.csh'.

Your can also email me (arif.harmanci@yale.edu) to generate multimappability profiles for new species.

<h2>Datasets</h2>
The ENCODE datasets can be downloaded from <a href="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/">UCSC Genome Browser</a>. <br>
The H3K36me3 datasets for K562 and GM12878 cell lines can be downloaded from <a href="http://archive.gersteinlab.org/proj/MUSIC/h3k36me3.tar.bz2">here</a>.
<!--------------- The enriched regions identified by MUSIC for the HMs and Polymerase ChIP-Seq datasets for several cell lines can be downloaded from <a href="multiMappability_signals/">here</a>. ---->
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<b><i>Arif Harmanci (arif.harmanci@yale.edu), Joel Rozowsky, Mark Gerstein, 2014</i></b>
</div>
</font>
</html>