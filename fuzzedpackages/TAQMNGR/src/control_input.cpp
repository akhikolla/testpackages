#include <taq_import.h>

int ControlInput(string DirIn, string symbol, string seconds, int StYear, int StMonth, int StDay, int EYear, int EMonth, int EDay, int MissCont, vector<string> &ToImp)
{
	vector<string> Mesi;
	int NowYear;
	
	Mesi.push_back("January");
	Mesi.push_back("February");
	Mesi.push_back("March");
	Mesi.push_back("April");
	Mesi.push_back("May");
	Mesi.push_back("June");
	Mesi.push_back("July");
	Mesi.push_back("August");
	Mesi.push_back("September");
	Mesi.push_back("October");
	Mesi.push_back("November");
	Mesi.push_back("December");
	
	time_t now = time(0);
	tm *ltm = localtime(&now);
	NowYear = 1900 + ltm->tm_year;
	
	string DirAgg = DirIn + separator + "Aggregate";
	if(VerificaDir(DirAgg.c_str()))
	{
		cout << "FATAL ERROR: folder " << DirAgg << " not found." << endl;
		cout << "REMEMBER TO AGGREGATE DATA BEFORE IMPORTING." << endl;
		return 0;
	}
	string DirSym = DirAgg + separator + symbol;
	if(VerificaDir(DirSym.c_str()))
	{
		cout << "FATAL ERROR: folder " << DirSym << " not found." << endl;
		cout << "REMEMBER TO AGGREGATE DATA BEFORE IMPORTING." << endl;
		return 0;
	}
	string DirSec = DirSym + separator + seconds;
	if(VerificaDir(DirSec.c_str()))
	{
		cout << "FATAL ERROR: folder " << DirSec << " not found." << endl;
		cout << "REMEMBER TO AGGREGATE DATA BEFORE IMPORTING." << endl;
		return 0;
	}

	int err = 0;
	if(StYear > NowYear)
	{
		cout << "FATAL ERROR: for starting year insert an integer value less then or equal to " << NowYear << "." << endl;
		err++;
	}
	if(StMonth < 1 || StMonth > 12)
	{
		cout << "FATAL ERROR: for starting month insert an integer value between 1 and 12." << endl;
		err++;
	}
	if(EYear > NowYear)
	{
		cout << "FATAL ERROR: for ending year insert an integer value less then or equal to " << NowYear << "." << endl;
		err++;
	}
	if(EMonth < 1 || EMonth > 12)
	{
		cout << "FATAL ERROR: for ending month insert an integer value between 1 and 12." << endl;
		err++;
	}
	if(EYear < StYear)
	{
		cout << "FATAL ERROR: the ending date precedes the starting date." << endl;
		err++;
	}
	else if(StYear == EYear && EMonth < StMonth)
	{
		cout << "FATAL ERROR: the ending date precedes the starting date." << endl;
		err++;
	}
	if(err) return 0;
	
	int nImp, nMiss;
	vector<string> Miss;
	
	if(ListaToImport(DirSec, StYear, StMonth, StDay, EYear, EMonth, EDay, nImp, ToImp, nMiss, Miss)) return 0;
	if(nImp == 0)
	{
		cout << "FATAL ERROR: no importable files found." << endl;
		return 0;
	}
	return nImp;
}
