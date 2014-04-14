#include <stdio.h>
#include <stdlib.h>
#include "ms_ansi_cli.h"
#include "ms_ansi_string.h"

/* 
INP: Standalone Inputs without any options, indicated by the preceding prefix.
OPT: Options, or id's, for values.
VAL: Value for a certain option.
*/ 
t_ansi_cli::t_ansi_cli(int argc, 
					   char* argv[], 
					   const char* _option_prefix)
{
	this->options = new vector<char*>();
	this->values = new vector<char*>();
	this->inputs = new vector<char*>();

	if(_option_prefix != NULL)
	{
		this->option_prefix = (char*)malloc(sizeof(char) * (t_string::string_length(_option_prefix) + 2));
		t_string::copy(this->option_prefix, _option_prefix);
	}
	else
	{
		printf("Cannot use NULL option prefix for parsing command line @ %s(%d).\n", __FILE__, __LINE__);
		exit(0);
	}

	// Copy the exe file name.
	this->exe_cmd = (char*)malloc(sizeof(char) * (t_string::string_length(argv[0]) + 2));
	t_string::copy(this->exe_cmd, argv[0]);

	// Go over all the arguments.
	//for(int i_arg = 1; i_arg < argc; i_arg++)
	enum{START, INP, VAL, OPT};
	int i_arg = 1;
	char state = START;
	while(i_arg < argc)
	{
		// New argument to be parsed.
		char* cur_arg = argv[i_arg];

		if(state == START)
		{
			// Is this an option?
			if(this->is_option(cur_arg))
			{
				char* new_option = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_option, cur_arg);
				this->options->push_back(new_option); // Must also push a new value for this option.

				state = OPT;
				i_arg++;
			}
			else
			{
				// Has to be an input.
				char* new_input = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_input, cur_arg);
				this->inputs->push_back(new_input); // Must also push a new value for this option.

				state = INP;
				i_arg++;
			}
		}
		else if(state == INP)
		{
			// Option or a new input.
			// Is this an option?
			if(this->is_option(cur_arg))
			{
				char* new_option = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_option, cur_arg);
				this->options->push_back(new_option); // Must also push a new value for this option.

				state = OPT;
				i_arg++;
			}
			else
			{
				// Has to be an input.
				char* new_input = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_input, cur_arg);
				this->inputs->push_back(new_input); // Must also push a new value for this option.

				state = INP;
				i_arg++;
			}
		}
		else if(state == VAL)
		{
			// Option or a new input.
			// Is this an option?
			if(this->is_option(cur_arg))
			{
				char* new_option = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_option, cur_arg);
				this->options->push_back(new_option); // Must also push a new value for this option.

				state = OPT;
				i_arg++;
			}
			else
			{
				// Has to be an input.
				char* new_input = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_input, cur_arg);
				this->inputs->push_back(new_input); // Must also push a new value for this option.

				state = INP;
				i_arg++;
			}
		}
		else if(state == OPT)
		{
			// Option or a new value. Cannot go to an input!
			// Is this an option?
			if(this->is_option(cur_arg))
			{
				char* new_option = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_option, cur_arg);
				this->options->push_back(new_option); // Must also push a new value for this option.

				this->values->push_back(NULL); // OPT->OPT transition: Implies that current option has no values, push NULL for the value of that option.

				state = OPT;
				i_arg++;
			}
			else
			{
				// Has to be a value.
				char* new_value = (char*)malloc(sizeof(char) * (t_string::string_length(cur_arg) + 2));
				t_string::copy(new_value, cur_arg);
				this->values->push_back(new_value); // Must also push a new value for this option.

				state = VAL;
				i_arg++;
			}
		}
		else
		{
			printf("Could not understand the state while parsing command line.\n");
			exit(0);
		}
	} // i_arg loop.

	if(state == OPT)
	{
		this->values->push_back(NULL);
	}

	if(this->options->size() != this->values->size())
	{
		printf("The size of options is not same as size of values %ld, %ld!\n", this->options->size(), this->values->size());
		exit(0);
	}

	//printf("Options - Values: (%d)\n", this->options->size());
	//for(int i = 0; i < this->options->size(); i++)
	//{
	//	printf("%s: %s\n", this->options->at(i), this->values->at(i));
	//}

	//printf("Inputs: (%d)\n", this->inputs->size());
	//for(int i = 0; i < this->inputs->size(); i++)
	//{
	//	printf("%s\n", this->inputs->at(i));
	//}
}

bool t_ansi_cli::is_option(const char* arg)
{
	t_string* test_str = new t_string(arg);
	if(test_str->starts_with(this->option_prefix))
	{
		delete(test_str);
		return(true);
	}
	else
	{
		delete(test_str);
		return(false);
	}
}

t_ansi_cli::~t_ansi_cli()
{
	for(int i_opt = 0; i_opt < (int)this->options->size(); i_opt++)
	{
		free(this->options->at(i_opt));
	}
	delete(this->options);

	for(int i_val = 0; i_val < (int)this->values->size(); i_val++)
	{
		free(this->values->at(i_val));
	}
	delete(this->values);

	for(int i_ip = 0; i_ip < (int)this->inputs->size(); i_ip++)
	{
		free(this->inputs->at(i_ip));
	}
	delete(this->inputs);

	if(this->option_prefix != NULL)
	{
		free(this->option_prefix);
	}

	if(this->exe_cmd != NULL)
	{
		free(this->exe_cmd);
	}
}

/*
Argument "success" indicates if the option is found or not.
*/
char* t_ansi_cli::get_value_by_option(const char* option2search, bool& success)
{
	success = false;
	for(int i_opt = 0; i_opt < (int)this->options->size(); i_opt++)
	{
		if(t_string::compare_strings(this->options->at(i_opt), option2search))
		{
			success = true;
			return(this->values->at(i_opt)); // Can be NULL!
		}
	}

	success = false;
	return(NULL);
}

bool t_ansi_cli::is_flag_set(const char* flag)
{
	bool success = false;

	char* ret = NULL;
	ret = this->get_value_by_option(flag, success);

	// If this option exists with no values, it is a flag.
	if(success && ret == NULL)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

vector<char*>* t_ansi_cli::get_input_list()
{
	return(this->inputs);
}


