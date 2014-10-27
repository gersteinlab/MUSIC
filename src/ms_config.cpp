#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ms_ansi_string.h"
#include "ms_utils.h"
#include "ms_config.h"
#include <string.h>

t_config::t_config(const char* config_fp)
{
	FILE* f_conf = open_f(config_fp, "r");
	if(f_conf == NULL)
	{
		printf("Could not open configuration file %s\n", config_fp);
		exit(0);
	}

	this->ids = new vector<char*>();
	this->vals = new vector<vector<char*>*>();
	//char cur_id[1000];
	//char cur_val[2000];

	char cur_line[5000];
	while(fgets(cur_line, 5000, f_conf) != NULL)
	{
		// Get rid of the new line.
		int l_line = strlen(cur_line);
		if(cur_line[l_line-1] == '\n')
		{
			cur_line[l_line-1] = 0;
		}

		if(cur_line[0] != '#')
		{
			t_string* line_str = new t_string(cur_line);
			t_string_tokens* line_tokens = line_str->tokenize_by_chars(" \t");

			// Add all the values in this line as a new entry.
			if((int)line_tokens->size() < 2)
			{
				//printf("Empty entry: %s\n", cur_line);
			}
			else
			{
				char* new_id = new char[strlen(line_tokens->at(0)->str()) + 2];
				strcpy(new_id, line_tokens->at(0)->str());
				vector<char*>* new_val_list = new vector<char*>();

				// Add all the values to the value list.
				for(int i_val = 1; i_val < (int)line_tokens->size(); i_val++)
				{
					char* new_val = new char[strlen(line_tokens->at(i_val)->str()) + 2];
					strcpy(new_val, line_tokens->at(i_val)->str());
					new_val_list->push_back(new_val);
				} // i_val loop.

				// Add the new entries.
				this->ids->push_back(new_id);
				this->vals->push_back(new_val_list);
			}
/*
			if(sscanf(cur_line, "%s %s", cur_id, cur_val) == 2)
			{
				char* new_id = new char[strlen(cur_id) + 2];
				char* new_val = new char[strlen(cur_val) + 2];
				strcpy(new_id, cur_id);
				strcpy(new_val, cur_val);
				this->ids->push_back(new_id);
				this->vals->push_back(new_val);
			} // Skip the comments in the configuration file.
*/
		} // Skip comments.
	} // File reading loop.

	fclose(f_conf);
}

t_config::~t_config()
{
}

bool t_config::get_double_val(const char* id, double& d_val)
{
	for(int i_id = 0; i_id < (int)this->ids->size(); i_id++)
	{
		if(strcmp(this->ids->at(i_id), id) == 0)
		{
			d_val = atof(this->vals->at(i_id)->at(0));
			return(true);
		}
	}
	return(false);
}

bool t_config::get_str_val(const char* id, char* val_buff)
{
        for(int i_id = 0; i_id < (int)this->ids->size(); i_id++)
        {
                if(strcmp(this->ids->at(i_id), id) == 0)
                {
			strcpy(val_buff, this->vals->at(i_id)->at(0));
                        return(true);
                }
        }
	return(false);
}

bool t_config::get_int_val(const char* id, int& int_val)
{
        for(int i_id = 0; i_id < (int)this->ids->size(); i_id++)
        {
                if(strcmp(this->ids->at(i_id), id) == 0)
                {
                        int_val = atoi(this->vals->at(i_id)->at(0));
                        return(true);
                }
        }
	return(false);
}

vector<vector<char*>*>* t_config::get_all_entries_per_id(const char* id)
{
	vector<vector<char*>*>* all_entries = new vector<vector<char*>*>();
    for(int i_id = 0; i_id < (int)this->ids->size(); i_id++)
    {
            if(strcmp(this->ids->at(i_id), id) == 0)
            {
				all_entries->push_back(this->vals->at(i_id));
            }
    }
    return(all_entries);
}

vector<char*>* t_config::get_val_list(const char* id)
{
    for(int i_id = 0; i_id < (int)this->ids->size(); i_id++)
    {
            if(strcmp(this->ids->at(i_id), id) == 0)
            {
                    return(this->vals->at(i_id));
            }
    }
    return(NULL);
}



