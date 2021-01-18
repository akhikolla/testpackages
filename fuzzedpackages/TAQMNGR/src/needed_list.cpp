#include <taq_import.h>

int NeededList(int StaY, int StaM, int StaD, int EndY, int EndM, int EndD, vector<string> &NeedFile)
{
	int conta = 0, Year = StaY, Month = StaM, Day = StaD;
	stringstream ss;
	
	do
	{
		if(Month < 10) ss << Year << "0" << Month;
		else ss << Year << Month;
		if(Day < 10) ss << "0" << Day << ".txt.gz";
		else ss << Day << ".txt.gz";
		NeedFile.push_back(ss.str());
		conta++;
		if(Year < EndY && Month < 12 && Day == 31)
		{
			Day = 1;
			Month++;
		}
		else if(Year < EndY && Month == 12 && Day == 31)
		{
			Day = 1;
			Month = 1;
			Year++;
		}
		else if(Year == EndY && Month < EndM && Day == 31)
		{
			Day = 1;
			Month++;
		}
		else if(Year == EndY && Month == EndM && Day == EndD) Year++;
		else Day++;
		ss.clear();
		ss.str(string());
	}while(Year <= EndY);
	
	return conta;
}
