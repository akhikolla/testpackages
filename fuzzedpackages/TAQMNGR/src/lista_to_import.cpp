#include <taq_import.h>

int ListaToImport(string Directory, int StrYear, int StrMonth, int StrDay, int EndYear, int EndMonth, int EndDay, int &nToImport, vector<string> &ToImport, int &nMissing, vector<string> &Missing)
{
	vector<string> SoloFiles, DirFiles, NeedFiles;
	int nPresent, nNeed, Dum;
	int i, j;
	
	nMissing = 0;
	nToImport = 0;
	
	nPresent = ListaFile(Directory, SoloFiles, DirFiles);
	if(nPresent == -1) return 1;
	nNeed = NeededList(StrYear, StrMonth, StrDay, EndYear, EndMonth, EndDay, NeedFiles);
	
	for(i = 0; i < nNeed; i++)
	{
		Dum = 0;
		for(j = 0; j < nPresent; j++)
		{
			if(!((NeedFiles[i]).compare(SoloFiles[j])))
			{
				Dum = 1;
				break;
			}
		}
		if(Dum)
		{
			ToImport.push_back(DirFiles[j]);
			nToImport++;
		}
		else
		{
			Missing.push_back(NeedFiles[i]);
			nMissing++;
		}
	}
	return 0;
}
