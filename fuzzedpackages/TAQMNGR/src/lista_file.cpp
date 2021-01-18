#include <taq.h>

int ListaFile(string dir, vector<string> &solofile, vector<string> &dirfile)
{
	DIR *dp = NULL;
	struct dirent *dirp = NULL;
	string word;
	int count = 0;

	dp = opendir(dir.c_str());
	if (dp == NULL)
	{
		cout << "FATAL ERROR: unable to list the file in " << dir << endl;
		return -1;
	}


	while ((dirp = readdir(dp)) != NULL)
	{
		word = string(dirp->d_name);
		if(word.compare(".") != 0 && word.compare("..") != 0 && word.compare(".daily.log"))
		{
			solofile.push_back(word);
			dirfile.push_back(dir + separator + word);
			count++;
		}
	}
	closedir(dp);
	return count;
}
