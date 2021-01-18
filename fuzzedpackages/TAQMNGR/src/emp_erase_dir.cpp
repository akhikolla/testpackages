#include <taq.h>


int EmpEraseDir(string DirPath)
{
	DIR *dp = NULL;
	struct dirent *dirp = NULL;
	string item, erase;

	dp = opendir(DirPath.c_str());
	if (dp == NULL)
	{
		cout << "FATAL ERROR: unable to open the directory " << DirPath << endl;
		return -1;
	}

	while ((dirp = readdir(dp)) != NULL)
	{
		item = string(dirp->d_name);
		if(item.compare(".") != 0 && item.compare("..") != 0)
		{
			erase = DirPath + separator + item;
			if(remove(erase.c_str()) != 0)
			{
				cout << "FATAL ERROR: unable to delete file" << erase << endl;
				return -1;
			}
		}
	}
	if(rmdir(DirPath.c_str()) != 0)
	{
		cout << "FATAL ERROR: unable to remove the directory" << DirPath << endl;
		return -1;
	}
	return 0;
}
