all: MUSIC

CC = g++
comp_flags = -c -Wall -O3
exec_name = bin/MUSIC
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/main.o \
${LIB_DIR}/ms_peak_calling_utils.o \
${LIB_DIR}/ms_mapped_read_tools.o \
${LIB_DIR}/ms_nomenclature.o \
${LIB_DIR}/ms_gsl_fft_filter_utils.o \
${LIB_DIR}/ms_annot_region_tools.o \
${LIB_DIR}/ms_signal_track_tools.o \
${LIB_DIR}/ms_signal_enrichment_utils.o \
${LIB_DIR}/ms_genome_sequence_tools.o \
${LIB_DIR}/ms_profile_normalization.o \
${LIB_DIR}/ms_combinatorics.o \
${LIB_DIR}/ms_nucleotide.o \
${LIB_DIR}/ms_min_max_utils.o \
${LIB_DIR}/ms_xlog_math.o \
${LIB_DIR}/ms_utils.o \
${LIB_DIR}/ms_ansi_cli.o \
${LIB_DIR}/ms_ansi_string.o \

MUSIC: ${objs}
	@echo "Building..."
	@${CC} -O3 -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 
