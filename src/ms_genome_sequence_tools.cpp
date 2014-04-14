#include <stdio.h>
#include <stdlib.h>
#include "ms_genome_sequence_tools.h"
#include "ms_utils.h"
#include "ms_ansi_string.h"
#include "ms_genomics_coords.h"
#include "ms_annot_region_tools.h"
#include "ms_xlog_math.h"
#include "ms_nomenclature.h"
#include "ms_nucleotide.h"
#include <vector>
#include <algorithm>
#include <ctype.h>
#include <string.h>

using namespace std;

bool _DUMP_GENOME_SEQUENCE_MESSAGES_ = false;

int map_nuc_to_sig(char nuc)
{
	switch(nuc)
	{
		case 'A':
			return(0);
			break;
		case 'C':
			return(1);
			break;
		case 'G':
			return(2);
			break;
		case 'U':
			return(0);
			break;
		case 'T':
			return(4);
			break;
		default:
			return(0);
			break;
	}
}

void get_numerized_sequence_signals(char* sequence, int l_signal,
	double* walk_signal,
	double* nucleotide_signal,
	double* gc_content_signal)
{
	fprintf(stderr, "Generating the sequence signal for %d nucleotides.\n", l_signal);

	double gc_count = 0;
	int cur_cumul_walk_sig = 0;
	double cur_gc_content = 0;
	double l_sequenced = 0;
	for(int i_nuc = 0; i_nuc < l_signal; i_nuc++)
	{
		cur_gc_content = 0;
		if(toupper(sequence[i_nuc]) != 'N')
		{
			l_sequenced++;

			if(toupper(sequence[i_nuc]) == 'C' || 
						toupper(sequence[i_nuc]) == 'G')
			{
				cur_gc_content = 1;
				gc_count++;
			}
		}

		if(toupper(sequence[i_nuc]) == 'A' || toupper(sequence[i_nuc]) == 'G' || 
		toupper(sequence[i_nuc]) == 'C' || toupper(sequence[i_nuc]) == 'T' ||
		toupper(sequence[i_nuc]) == 'U')
		{
			if(toupper(sequence[i_nuc]) == 'A' || toupper(sequence[i_nuc]) == 'G')
			{
				cur_cumul_walk_sig--;
			}
			else
			{
				cur_cumul_walk_sig++;
			}
		}

		gc_content_signal[i_nuc+1] = cur_gc_content;
		walk_signal[i_nuc+1] = cur_cumul_walk_sig;
		nucleotide_signal[i_nuc+1] = nuc_2_num(toupper(sequence[i_nuc]));
	} // i_nuc loop.

	fprintf(stderr, "GC percentage: %lf\n", gc_count / l_sequenced);
}

// Given the 0-based seq_i indexing, returns the genome index with the assumed index convention.
bool get_genome_i_per_seq_i(char strand, 
	int genome_start, 
	int genome_end, 
	int& genome_i, int seq_i)
{
	int l_seq = genome_end - genome_start + 1;
	genome_i = seq_i + genome_start;

	// Revert the index if the mirna is on negative strand.
	if(strand == '-')
	{
		//seq_i = l_seq - seq_i - 1;
		genome_i = l_seq + genome_start - seq_i - 1;
	}

	return(true);
}

// This returns 0 based index for the sequence.
bool get_seq_i_per_genome_i(char strand, 
	int genome_start, 
	int genome_end, 
	int genome_i, int& seq_i)
{
	// genome_i must be within the element of interest.
	if(genome_i < genome_start ||
		genome_i > genome_end)
	{
		return(false);
	}

	// genome_i must be within the element of interest: If it is not, fix it.
	if(genome_i < genome_start)
	{
		genome_i = genome_start;
	}
	
	if(genome_i > genome_end)
	{
		genome_i = genome_end;
	}

	int l_seq = genome_end - genome_start + 1;
	seq_i = genome_i - genome_start;

	// Revert the index if the mirna is on negative strand.
	if(strand == '-')
	{
		seq_i = l_seq - seq_i - 1;
	}

	return(true);
}

//void binarize_buffered_fasta_file(char* fasta_fp, char* bin_dir)
void binarize_fasta_file(char* fasta_fp, char* bin_dir)
{
	printf("Binarizing %s.\n", fasta_fp);
	FILE* f_fasta = open_f(fasta_fp, "r");

	char* cur_line = NULL;
	char* cur_entry_buffer = new char[250 * 1000 * 1000];
	int cur_entry_i = 0;
	char cur_entry_id[1000];

	while(1)
	{
		// Process the current buffer.
		cur_line = getline(f_fasta);

		if(cur_line == NULL)
		{
			// File ended, dump the last entry if there are values in it.
			if(cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				fclose(f_bin);
			}
			break;
		}
		else if(cur_line[0] == '>')
		{
			if(cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin", bin_dir, cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				fclose(f_bin);
			}

			// Update the id, reset the counter.
			strcpy(cur_entry_id, &(cur_line[1]));
			cur_entry_i = 0;
		}
		else
		{
			// Concatenate the current sequence line.
			int l_cur_line = t_string::string_length(cur_line);
			for(int i = 0; i < l_cur_line; i++)
			{
				cur_entry_buffer[cur_entry_i] = cur_line[i];
				cur_entry_i++;
			} // i loop.
		}

		delete [] cur_line;
	} // file reading loop.

	fclose(f_fasta);

	delete [] cur_entry_buffer;
}

//char* numerize_sequence_signal(char* seq_signal, int l_seq)
//{
//	char* numerized_seq_buffer = new char[l_seq+1];
//
//	for(int i_nuc = 0; i_nuc < l_seq; i_nuc++)
//	{
//		numerized_seq_buffer[i_nuc] = nuc_2_num(seq_signal[i_nuc]);
//	} // i_nuc loop.
//
//	return(numerized_seq_buffer);
//}

char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq)
{
	FILE* f_bin_seq = open_f(bin_seq_fp, "rb");
	int l_bin_seq = 0;

	// Read the length.
	fread(&l_bin_seq, sizeof(int), 1, f_bin_seq);

	l_seq = l_bin_seq;

	// Read the sequence.
	char* bin_seq = new char[l_bin_seq+2];
	fread(bin_seq, sizeof(char), l_bin_seq, f_bin_seq);
	fclose(f_bin_seq);

	return(bin_seq);
}

vector<t_annot_region*>* load_BED_with_ancestral_derived_seqs(char* bed_w_seq_fp)
{
	vector<t_annot_region*>* regs_w_seqs = new vector<t_annot_region*>();

	FILE* f_bed_w_seq = open_f(bed_w_seq_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_bed_w_seq);
		if(cur_line == NULL)
		{
			break;
		}

		char chrom[100];
		int start;
		int end;
		char strand;
		char* name = new char[t_string::string_length(cur_line) + 2];
		char** cur_anc_der_seqs = new char*[2];
		char* cur_ancestral_seq = new char[t_string::string_length(cur_line) + 2];
		char* cur_derived_seq = new char[t_string::string_length(cur_line) + 2];
		cur_anc_der_seqs[0] = cur_ancestral_seq;
		cur_anc_der_seqs[1] = cur_derived_seq;

		if(sscanf(cur_line, "%s %d %d %s %*s %c %s %s", chrom, &start, &end, name, &strand, 
			cur_ancestral_seq, cur_derived_seq) != 7)
		{
			fprintf(stderr, "Could not parse the sequence bed file line:\n%s\n", cur_line);
			exit(0);
		}

		t_annot_region* cur_region = get_empty_region();
		cur_region->chrom = t_string::copy_me_str(chrom);
		cur_region->start = start;
		cur_region->end = end;
		cur_region->name = t_string::copy_me_str(name);
		cur_region->strand = strand;

		// The sequence is stored in the data entry.
		cur_region->data = (void*)cur_anc_der_seqs;

		regs_w_seqs->push_back(cur_region);

		delete [] name;
	} // file reading loop.

	fclose(f_bed_w_seq);

	return(regs_w_seqs);
}

// Load the bed file with sequence information that is dumped by sequence extraction.
vector<t_annot_region*>* load_BED_with_sequences(char* bed_w_seq_fp)
{
	vector<t_annot_region*>* regs_w_seqs = new vector<t_annot_region*>();

	FILE* f_bed_w_seq = open_f(bed_w_seq_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_bed_w_seq);
		if(cur_line == NULL)
		{
			break;
		}

		char chrom[100];
		int start;
		int end;
		char strand;
		char name[100];
		char* cur_seq = new char[t_string::string_length(cur_line) + 2];		

		if(sscanf(cur_line, "%s %d %d %s %*s %c %s", chrom, 
			&start, &end, name, &strand, cur_seq) != 6)
		{
			fprintf(stderr, "Could not parse the sequence bed file line:\n%s\n", cur_line);
			exit(0);
		}

		t_annot_region* cur_region = get_empty_region();
		cur_region->chrom = t_string::copy_me_str(chrom);
		cur_region->start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_region->end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_region->strand = strand;

		// The sequence is stored in the data entry.
		cur_region->data = (void*)cur_seq;

		regs_w_seqs->push_back(cur_region);
	} // file reading loop.

	fclose(f_bed_w_seq);

	return(regs_w_seqs);
}

void batch_extract_interval_sequences(char* chr_data_dir, vector<t_annot_region*>* regions)
{
	t_sorted_annot_region_lists* retstr_regions = restructure_annot_regions(regions);

	for(int i_chr = 0; i_chr < (int)retstr_regions->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", retstr_regions->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_regions = retstr_regions->regions_per_chrom[i_chr];
		for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
		{
			// Process the current transcript: Merge the exons
			vector<t_annot_region*>* merged_exons = merge_annot_regions(cur_chr_regions->at(i_reg)->intervals, 0, true);

			// Sort the merged regions.
			sort(merged_exons->begin(), merged_exons->end(), sort_regions);

			t_string* cur_transcript_seq = new t_string();

			// Extract the sequence for all the exons.
			for(int i_ex = 0; i_ex < (int)merged_exons->size(); i_ex++)
			{
				if(merged_exons->at(i_ex)->strand == '-')
				{
					merged_exons->at(i_ex)->strand = '+';
				}

				char* cur_exon_seq = extract_sequence(chr_data_dir, retstr_regions->chr_ids->at(i_chr), merged_exons->at(i_ex)->start, merged_exons->at(i_ex)->end, merged_exons->at(i_ex)->strand);
				cur_transcript_seq->concat_string(cur_exon_seq);
				delete [] cur_exon_seq;
			} // i_ex loop.

			// If the transcript is on negative strand, revert the sequence.
			char* transcript_seq = t_string::copy_me_str(cur_transcript_seq->str());
			if(cur_chr_regions->at(i_reg)->strand == '-')
			{
				//cur_transcript_seq->revert();
				delete [] transcript_seq;
				transcript_seq = get_reverse_complement(cur_transcript_seq->str(), 0, t_string::string_length(cur_transcript_seq->str())-1);
			}

			// Copy the data.
			cur_chr_regions->at(i_reg)->data = transcript_seq;

			delete(cur_transcript_seq);
			// Memory for merged exons must be freeed.
		} // i_reg loop.
	} // i_chr loop.
}

/*
Extract the sequences for multiple regions at one time.
Start and end are 1 based.
*/
//vector<t_annot_region*>* batch_extract_region_sequences(char* chr_data_dir, char* regions_bed_fp)
void batch_extract_region_sequences(char* chr_data_dir, vector<t_annot_region*>* regions)
{
	//vector<t_annot_region*>* regions = load_BED(regions_bed_fp);

	vector<char*>* chr_ids = get_chr_ids(regions);

	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_regions = get_regions_per_chromosome(regions, chr_ids->at(i_chr));

		fprintf(stderr, "Extracting sequences in %s\n", chr_ids->at(i_chr));

		// Sort the regions.
		sort(cur_chr_regions->begin(), cur_chr_regions->end(), sort_regions);

		char chr_seq_fp[1000];
		sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, chr_ids->at(i_chr));
		if(check_file(chr_seq_fp))
		{
			FILE* f_chr = open_f(chr_seq_fp, "rb");

			// Go over all the regions.
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				fprintf(stderr, "%d (%ld)             \r", i_reg, cur_chr_regions->size());

				int l_reg = cur_chr_regions->at(i_reg)->end - cur_chr_regions->at(i_reg)->start + 1;
				int f_start = cur_chr_regions->at(i_reg)->start - 1 + sizeof(int);
				int f_end = f_start + l_reg - 1;

				int length = f_end - f_start + 1;

				char* seq_buffer = new char[length + 2];
				memset(seq_buffer, 0, sizeof(char) * (length + 1));

				// Read the data.
				fseek(f_chr, f_start, SEEK_SET);

				// Fill the buffer.
				fread(seq_buffer, length, 1, f_chr);

				// If the requested strand is on negative strand, do reverse complementing.
				if(cur_chr_regions->at(i_reg)->strand == '-')
				{
					// Replace the sequence buffer with its reverse complement.
					char* reverse_complement_seq = get_reverse_complement(seq_buffer, 0, t_string::string_length(seq_buffer) - 1);
					delete [] seq_buffer;
					seq_buffer = reverse_complement_seq;
				}

				//fprintf(stderr, "%s\t%d\t%d\t.\t.\t%c\t%s\n", cur_chr_regions->at(i_reg)->chrom, 
				//	cur_chr_regions->at(i_reg)->start, cur_chr_regions->at(i_reg)->end, 
				//	cur_chr_regions->at(i_reg)->strand, seq_buffer);
			
				// Save the sequence to data entry of region.
				cur_chr_regions->at(i_reg)->data = (void*)seq_buffer;
			} // i_reg.

			fprintf(stderr, "\n");

			fclose(f_chr);
		}
		else
		{
			fprintf(stderr, "Could not find %s\n", chr_seq_fp);
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				cur_chr_regions->at(i_reg)->data = NULL;
			} // i_reg loop.
		}
	} // i_chr loop.
}

char* load_whole_chromosome_sequence_per_bin_file(char* chr_fp, int& n_nucs)
{
	FILE* f_chr = open_f(chr_fp, "rb");
	
	// Get the file size.
	fseek (f_chr, 0 , SEEK_END);
	n_nucs = (int)(ftell(f_chr));
	rewind(f_chr);

	char* nucs = new char[n_nucs+2]; 

	fread(nucs, sizeof(char), n_nucs, f_chr);

	fclose(f_chr);

	return(nucs);
}

/*
Start and end are 1 based.
*/
char* extract_sequence(char* chr_data_dir, char* raw_chr_name, int start, int end, char strand, bool upper_case)
{
	char chr_seq_fp[1000];
	sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, raw_chr_name);
	FILE* f_chr = fopen(chr_seq_fp, "rb");

	if(f_chr == NULL)
	{
		// Try to fix the chromosome file path: Add 'chr' to beginning.
		sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, raw_chr_name);
		f_chr = fopen(chr_seq_fp, "rb");

		if(f_chr == NULL)
		{
			fprintf(stderr, "Could not identify the conforming chromosome name from id %s, returning NULL.\n", raw_chr_name);
			return(NULL);
		}
	}

	// Read the binary sequence length.
	int l_bin_seq = 0;
	fread(&l_bin_seq, sizeof(int), 1, f_chr);

	int f_start = start - 1 + sizeof(int);
	//int f_end = end - 1;
	int f_end = f_start + (end-start);

	char* seq_buffer = NULL;

	if(f_start < f_end)
	{
		int length = (end - start + 1);

		seq_buffer = new char[length + 2];
		memset(seq_buffer, 0, sizeof(char) * (length + 1));

		// Go to the starting position.
		fseek(f_chr, f_start, SEEK_SET);

		// Fill the buffer.
		fread(seq_buffer, sizeof(char), length, f_chr);

		seq_buffer[length] = 0;
	} // start, end coparison.
	else
	{
		// Get the file size.
		fseek (f_chr, 0 , SEEK_END);
		int length = (int)(ftell(f_chr));
		
		seq_buffer = new char[length+2]; 
		memset(seq_buffer, 0, sizeof(char) * (length + 1));

		rewind(f_chr);

		fread(seq_buffer, sizeof(char), length, f_chr);

		seq_buffer[length] = 0;
	} // whole chromosome reading check.

	// If the requested strand is on negative strand, do reverse complementing.
	if(strand == '-')
	{
		// Replace the sequence buffer with its reverse complement.
		char* reverse_complement_seq = get_reverse_complement(seq_buffer, 0, t_string::string_length(seq_buffer) - 1);
		delete [] seq_buffer;
		seq_buffer = reverse_complement_seq;
	}

	fclose(f_chr);

	if(upper_case)
	{
		t_string::to_upper(seq_buffer);
	}

	return(seq_buffer);
}

/*
For reverse strand, the transcribed RNA is exactly same as forward strand with the exception that T(t)'s
are converted to U(u)'s.
*/
char* get_reverse_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* reverse_RNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(reverse_RNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	int i_rna_nuc = 0;
	for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	{
		if(cur_aln_line[i_nuc] == 'T')
		{
			reverse_RNA[i_rna_nuc] = 'U';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			reverse_RNA[i_rna_nuc] = 'u';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'A' || cur_aln_line[i_nuc] == 'a' ||
			cur_aln_line[i_nuc] == 'C' || cur_aln_line[i_nuc] == 'c' ||
			cur_aln_line[i_nuc] == 'G' || cur_aln_line[i_nuc] == 'g' ||
			cur_aln_line[i_nuc] == '-')
		{
			reverse_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
		else
		{
			printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			reverse_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
	} // i_nuc loop.

	return(reverse_RNA);
}

char* get_reverse_complement(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* rev_comp_DNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(rev_comp_DNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	//printf("Forward strand RNA for strand:\n");
	//for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	//{
	//	printf("%c", cur_aln_line[i_nuc]);
	//}
	//printf("\n");

	int i_dna_nuc = 0;

	// The forward strand RNA is transcribed from 3' end.
	for(int i_nuc = i_nuc_max; i_nuc >= i_nuc_min; i_nuc--)
	{
		if(cur_aln_line[i_nuc] == 'A')
		{
			rev_comp_DNA[i_dna_nuc] = 'T';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'T')
		{
			rev_comp_DNA[i_dna_nuc] = 'A';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'G')
		{
			rev_comp_DNA[i_dna_nuc] = 'C';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'C')
		{
			rev_comp_DNA[i_dna_nuc] = 'G';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'a')
		{
			rev_comp_DNA[i_dna_nuc] = 't';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			rev_comp_DNA[i_dna_nuc] = 'a';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'g')
		{
			rev_comp_DNA[i_dna_nuc] = 'c';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'c')
		{
			rev_comp_DNA[i_dna_nuc] = 'g';
			i_dna_nuc++;
		}
		//else if(cur_aln_line[i_nuc] == '-')
		//{
		//	rev_comp_DNA[i_dna_nuc] = '-';
		//	i_dna_nuc++;
		//}
		else
		{
			//printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			rev_comp_DNA[i_dna_nuc] = cur_aln_line[i_nuc];
			i_dna_nuc++;
		}
	} // i_nuc loop.

	//printf("%s", rev_comp_DNA);
	//getc(stdin);
	
	return(rev_comp_DNA);
}

char* get_forward_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* forward_RNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(forward_RNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	//printf("Forward strand RNA for strand:\n");
	//for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	//{
	//	printf("%c", cur_aln_line[i_nuc]);
	//}
	//printf("\n");

	int i_rna_nuc = 0;

	// The forward strand RNA is transcribed from 3' end.
	for(int i_nuc = i_nuc_max; i_nuc >= i_nuc_min; i_nuc--)
	{
		if(cur_aln_line[i_nuc] == 'A')
		{
			forward_RNA[i_rna_nuc] = 'U';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'T')
		{
			forward_RNA[i_rna_nuc] = 'A';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'G')
		{
			forward_RNA[i_rna_nuc] = 'C';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'C')
		{
			forward_RNA[i_rna_nuc] = 'G';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'a')
		{
			forward_RNA[i_rna_nuc] = 'u';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			forward_RNA[i_rna_nuc] = 'a';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'g')
		{
			forward_RNA[i_rna_nuc] = 'c';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'c')
		{
			forward_RNA[i_rna_nuc] = 'g';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == '-')
		{
			forward_RNA[i_rna_nuc] = '-';
			i_rna_nuc++;
		}
		else
		{
			printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			forward_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
	} // i_nuc loop.

	//printf("%s", forward_RNA);
	//getc(stdin);
	
	return(forward_RNA);
}






