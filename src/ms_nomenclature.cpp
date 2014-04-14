#include <stdio.h>
#include <stdlib.h>
#include "ms_nomenclature.h"
#include "ms_ansi_string.h"

void normalize_chr_id(char* raw_chr_id)
{
	if(t_string::starts_with(raw_chr_id, "chr") || 
		t_string::starts_with(raw_chr_id, "Chr") ||
		t_string::starts_with(raw_chr_id, "CHR"))
	{
		// Get rid of the first 3 characters.
		char temp[1000];
		t_string::copy(temp, &raw_chr_id[3]);
		t_string::copy(raw_chr_id, temp);
	}
	else
	{
		// No need to change.
		return;
	}
}

