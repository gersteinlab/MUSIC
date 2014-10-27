ll: multiscale_peak_analyzer 

CC = g++
comp_flags = -c -Wall -O3
gsl_flags = -lgsl -lgslcblas
exec_name = ../../bin/multiscale_peak_analyzer
LIB_DIR = ../../lib

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = \
main.o \
peak_calling_utils.o \
${LIB_DIR}/genomics_utils/mapped_read/mapped_read_tools.o \
${LIB_DIR}/nomenclature/nomenclature.o \
${LIB_DIR}/genomics_utils/signal_processing/gsl_fft_filter_utils.o \
${LIB_DIR}/genomics_utils/signal_processing/pattern_analysis_utils.o \
${LIB_DIR}/utils/stats/histogram.o \
${LIB_DIR}/utils/stats/full_histogram.o \
${LIB_DIR}/genomics_utils/annotation/annot_region_tools.o \
${LIB_DIR}/genomics_utils/annotation/gff_utils.o \
${LIB_DIR}/genomics_utils/signal_track/signal_track_tools.o \
${LIB_DIR}/utils/combinatorics/combinatorics.o \
${LIB_DIR}/utils/geometry/multidimensional_geometry.o \
${LIB_DIR}/genomics_utils/signal_processing/signal_enrichment_utils.o \
${LIB_DIR}/genomics_utils/signal_processing/multiscale_peak_utils.o \
${LIB_DIR}/genomics_utils/signal_processing/min_max_utils.o \
${LIB_DIR}/genomics_utils/signal_processing/profile_normalization.o \
${LIB_DIR}/utils/rng/rng.o \
${LIB_DIR}/utils/rng/seed_manager.o \
${LIB_DIR}/utils/xmath/log/xlog_math.o \
${LIB_DIR}/utils/xmath/linear/linear_math.o \
${LIB_DIR}/utils/ansi_cli/config.o \
${LIB_DIR}/utils/file/utils.o \
${LIB_DIR}/utils/ansi_string/ansi_string.o \

multiscale_peak_analyzer: ${objs}
	${CC} -O3 ${gsl_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
