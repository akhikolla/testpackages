#include <taq.h>

int create_log(string dirOut, vector<string> sym, vector<string> date, vector<double> total, vector<double> notcorr_delayed, vector<double> BrownGallo, int cleaned_list_flag, vector<string> &allready_deleted)
{
	int day_cleaned = total.size();
	string file_name;
	string useless = dirOut + separator + ".useless.log";
	
	ofstream out;
	out.open(useless.c_str());
	hide_file(useless.c_str());
	
	for(int j = 1; j < day_cleaned; j++)
	{
		if((sym[j]).compare(sym[j-1]))
		{
			out.close();
			file_name = dirOut + separator + "Trade" + separator + (sym[j]) + separator + ".daily.log";

			if((!VerificaDir(file_name.c_str())) && cleaned_list_flag == 0 && check_deleted(file_name, allready_deleted))
			{
				if(remove(file_name.c_str()) != 0)
				{
					cout << "FATAL ERROR: unable to remove the file " << file_name << "." << endl;
					return 1;
				}
				allready_deleted.push_back(file_name);
			}
			
			out.open(file_name.c_str(), ios::app);
			hide_file(file_name.c_str());
		}
		
		out << date[j] << " " << total[j] << " " << notcorr_delayed[j] << " " << BrownGallo[j] << endl;
	}
	out.close();
	if(remove(useless.c_str()))
	{
		cout << "FATAL ERROR: unable to remove the file " << useless << "." << endl;
		return 1;
	}
	return 0;
}
