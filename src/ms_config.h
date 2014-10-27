#ifndef __CONFIG__
#define __CONFIG__

#include <vector>
using namespace std;

class t_config
{
public:
	t_config(const char* config_fp);
	~t_config();

	// ID's and values for the entries in the configuration file.
	vector<char*>* ids;
	vector<vector<char*>*>* vals; // A value list for each identifier.

	// Can oly access the entries using th id's.
	bool get_double_val(const char* id, double& d_val);
	bool get_str_val(const char* id, char* val_buff);
	bool get_int_val(const char* id, int& int_val);
	vector<char*>* get_val_list(const char* id);

	vector<vector<char*>*>* get_all_entries_per_id(const char* id);
};

#endif // __CONFIG__
