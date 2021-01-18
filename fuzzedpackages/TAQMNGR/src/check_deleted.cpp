#include <taq.h>

int check_deleted(string file, vector<string> deleted)
{
	int n_deleted = deleted.size();
	
	if(n_deleted == 0) return 1;
	else
	{
		for(int j = 0; j < n_deleted; j++)
		{
			if(!file.compare(deleted[j])) return 0;
		}
	}
	return 1;
}
	
	
