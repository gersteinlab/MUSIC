#include <stdio.h>
#include <stdlib.h>
#include "ms_nomenclature.h"
#include "ms_ansi_string.h"
#include <string.h>

void normalize_chr_id(char* raw_chr_id)
{
	t_string::replace_avoid_list(raw_chr_id, "+-./\*^$#@!&()~+=?\"`][}{", '_');

	if(t_string::starts_with(raw_chr_id, "chr") || 
		t_string::starts_with(raw_chr_id, "Chr") ||
		t_string::starts_with(raw_chr_id, "CHR"))
	{
		// Get rid of the first 3 characters.
		char temp[1000];
		strcpy(temp, &raw_chr_id[3]);
		strcpy(raw_chr_id, temp);
	}
	else
	{
		// No need to change.
		return;
	}
}

