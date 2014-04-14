#ifndef __NUCLEOTIDE__
#define __NUCLEOTIDE__

// Basic functionality about the nucleotide naming, pairing, numbering, etc.
bool check_rna_pairing(char nuc1, char nuc2);
char num_2_nuc(int num);
int nuc_2_num(char nuc);
bool is_valid_nuc_char(char nuc);
char get_transcribed_rna_nuc_per_dna_nuc(char dna_nuc);
char get_transcribing_dna_nuc_per_rna_nuc(char rna_nuc);

char get_complementary_dna_nuc_per_dna_nuc(char dna_nuc);
char get_complementary_rna_nuc_per_rna_nuc(char rna_nuc);

// Convert nucleotide symbols into indices: XACGUI -> 012345 
// refer to IUPAC nucleotide symbols for more information:
// http://www.mun.ca/biochem/courses/3107/symbols.html
void map_nuc_IUPAC_code(char raw_nuc, 
								char& trans_nuc, 
								int& num, 
								bool& force_unpaired);

char get_dna_pair(char nuc);

#endif // __NUCLEOTIDE__