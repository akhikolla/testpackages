#include <taq_aggregate.h>

int ListaFileToBeAggregated(string dir_root, string dir_cleaned, string symbol, string seconds, vector<string> &solofile, vector<string> &dirfile, int useAggregated)
{
	string dir_aggregated = dir_root + separator + "Aggregate" + separator + symbol + separator + seconds;
	vector<string> aggregated, to_be_aggregated, aggregatedDir, to_be_aggregatedDir;
	int nAgr, n2bAgr;
	int i, j;
	int ind;
	int count = 0;
		
	if(VerificaDir(dir_aggregated.c_str()) || useAggregated == 0)
	{
		return ListaFile(dir_cleaned, solofile, dirfile);
	}
		
	if(ListaFile(dir_aggregated, aggregated, aggregatedDir) == -1) return -1;
	if(ListaFile(dir_cleaned, to_be_aggregated, to_be_aggregatedDir) == -1) return -1;
		
	nAgr = aggregated.size();
	n2bAgr = to_be_aggregated.size();
	
	for(i = 0; i < n2bAgr; i++)
	{
		ind = 1;
		for(j = 0; j < nAgr; j++)
		{
			if(!((to_be_aggregated[i]).compare(aggregated[j])))
			{
				ind = 0;
				break;
			}
		}
		if(ind)
		{
			solofile.push_back(to_be_aggregated[i]);
			dirfile.push_back(to_be_aggregatedDir[i]);
			count++;
		}
	}
	return count;
}
