#include <taq.h>

int ListaConfFile(string dir, string ListofCleaned, vector<string> &solofileConf, vector<string> &dirfileConf)
{
	DIR *dp = NULL;
	struct dirent *dirp = NULL;
	string word;
	vector<string> solofile, dirfile, solofilecln;
	int count = 0, countConf = 0;
	int i, j;

	dp = opendir(dir.c_str());
	if (dp == NULL)
	{
		cout << "FATAL ERROR: unable to list the file in " << dir << endl;
		return -1;
	}


	while ((dirp = readdir(dp)) != NULL)
	{
		word = string(dirp->d_name);
		if(word.compare(".") != 0 && word.compare("..") != 0)
		{
			solofile.push_back(word);
			dirfile.push_back(dir + separator + word);
			count++;
		}
	}
	closedir(dp);

	ifstream incln(ListofCleaned.c_str());
	while(incln.good())
	{
		getline(incln, word, '\n');
		if(!word.empty()) solofilecln.push_back(word);
	}

	int ind, dimcln = solofilecln.size();
	for(i = 0; i < count; i++)
	{
		ind = 1;
		for(j = 0; j < dimcln; j++)
		{
			if(!((solofile[i]).compare(solofilecln[j])))
			{
				ind = 0;
				break;
			}
		}
		if(ind)
		{
			solofileConf.push_back(solofile[i]);
			dirfileConf.push_back(dirfile[i]);
			countConf++;
		}
	}
	return countConf;
}
