#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "ms_nucleotide.h"
#include <string.h>

char get_dna_pair(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return('T');
	}
	else if(toupper(nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(nuc) == 'T')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

bool check_rna_pairing(char nuc1, char nuc2)
{
	if(toupper(nuc1) == 'A')
	{
		if(toupper(nuc2) == 'U' ||
			toupper(nuc2) == 'T')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'C')
	{
		if(toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'G')
	{
		if(toupper(nuc2) == 'U' || 
			toupper(nuc2) == 'T' || 
			toupper(nuc2) == 'C')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'T')
	{
		if(toupper(nuc2) == 'A' ||
			toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	if(toupper(nuc1) == 'U')
	{
		if(toupper(nuc2) == 'A' ||
			toupper(nuc2) == 'G')
		{
			return(true);
		}
	}

	return(false);
}

bool is_valid_nuc_char(char nuc)
{
	if(toupper(nuc) == 'A' ||
		toupper(nuc) == 'C' ||
		toupper(nuc) == 'G' ||
		toupper(nuc) == 'U' ||
		toupper(nuc) == 'T')
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

int nuc_2_num(char nuc)
{
	if(toupper(nuc) == 'A')
	{
		return(0);
	}
	else if(toupper(nuc) == 'C')
	{
		return(1);
	}
	else if(toupper(nuc) == 'G')
	{
		return(2);
	}
	else if(toupper(nuc) == 'U' ||
		toupper(nuc) == 'T')
	{
		return(3);
	}
	else
	{
		return(4);
		//fprintf(stderr, "Cannot convert %c to number.\n", nuc);
		//exit(0);
	}
}

char num_2_nuc(int num)
{
	char nucs[] = "ACGU";

	char nuc = nucs[num];

	return(nuc);
}

// Get the rna nucleotide corresponding to a given dna nucleotide.
char get_transcribed_rna_nuc_per_dna_nuc(char dna_nuc)
{
	if(toupper(dna_nuc) == 'A' || 
		toupper(dna_nuc) == 'C' || 
		toupper(dna_nuc) == 'G')
	{
		return(toupper(dna_nuc));
	}
	else if(toupper(dna_nuc) == 'T')
	{
		return('U');
	}
	else
	{
		return(toupper(dna_nuc));
	}
}

char get_transcribing_dna_nuc_per_rna_nuc(char rna_nuc)
{
	if(toupper(rna_nuc) == 'A' || 
		toupper(rna_nuc) == 'C' || 
		toupper(rna_nuc) == 'G')
	{
		return(toupper(rna_nuc));
	}
	else if(toupper(rna_nuc) == 'U')
	{
		return('T');
	}
	else
	{
		return(toupper(rna_nuc));
	}
}

char get_complementary_dna_nuc_per_dna_nuc(char dna_nuc)
{
	if(toupper(dna_nuc) == 'A')
	{
		return('T');
	}
	else if(toupper(dna_nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(dna_nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(dna_nuc) == 'T' ||
		toupper(dna_nuc) == 'U')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

char get_complementary_rna_nuc_per_rna_nuc(char rna_nuc)
{
	if(toupper(rna_nuc) == 'A')
	{
		return('U');
	}
	else if(toupper(rna_nuc) == 'C')
	{
		return('G');
	}
	else if(toupper(rna_nuc) == 'G')
	{
		return('C');
	}
	else if(toupper(rna_nuc) == 'U')
	{
		return('A');
	}
	else
	{
		return('N');
	}
}

/*

A 
A adenine 

C 
C cytosine 

G 
G guanine 

T 
T thymine 

U 
U uracil 
   

R 
A or G purine 

Y 
C or T (U) pyrimidine 
   

M 
A or C amino 

K 
G or T (U) keto 

S 
C or G strong (3 H bonds) 

W 
A or T (U) weak (2 H bonds) 
   

B 
C or G or T (U) not A 

D 
A or G or T (U) not C 

H 
A or C or T (U) not G 

V 
A or C or G not T (U) 
   

N 
A or C or G or T (U) any nucleotide 

*/
void map_nuc_IUPAC_code(char raw_nuc, 
						char &trans_nuc, 
						int &num, 
						bool& force_unpaired)
{
	if(raw_nuc == 'a' || raw_nuc == 'c' || raw_nuc == 'g' || raw_nuc == 'u' || raw_nuc == 't')
	{
		force_unpaired = true;
	}
	else
	{
		force_unpaired = false;
	}

	if (toupper(raw_nuc) == 'A') 
	{
		trans_nuc = raw_nuc;
		num=1;
	}
	else if(toupper(raw_nuc) == 'B')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'C')
	{
		trans_nuc = raw_nuc;
		num = 2;
	}
	else if(toupper(raw_nuc) == 'D')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'G')
	{
		trans_nuc = raw_nuc;
		num = 3;
	}
	else if(toupper(raw_nuc) == 'H')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'I')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'K')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'M')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'N')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'R')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'S')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'T')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'U')
	{
		trans_nuc = raw_nuc;
		num = 4;
	}
	else if(toupper(raw_nuc) == 'V')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else if(toupper(raw_nuc) == 'W')
	{
		trans_nuc = 'N';
		num = 0;
	}
        else if(toupper(raw_nuc) == 'X')
        {
                trans_nuc = 'N';
                num = 0;
        }
	else if(toupper(raw_nuc) == 'Y')
	{
		trans_nuc = 'N';
		num = 0;
	}
	else
	{
		trans_nuc = 'N';
		num = 0;
	}

	if(num == 0)
	{
		//printf("Found %c\n", raw_nuc);
		//getc(stdin);
	}
}
