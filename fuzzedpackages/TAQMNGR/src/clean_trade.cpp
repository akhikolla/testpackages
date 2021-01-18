#include <taq.h>


int CleanTrade(const char *NomeFile, string DirOut, string DirTemp, int win, double delta, double gra, int flag, vector<string> &deleted_log)
{
	vector<string> header;
	header.push_back("SYMBOL");
	header.push_back("DATE");
	header.push_back("TIME");
	header.push_back("PRICE");
	header.push_back("G127");
	header.push_back("CORR");
	header.push_back("COND");
	header.push_back("EX");
	header.push_back("SIZE");
	
	vector<string> MtrRows;
	vector<double> MtrPrc, SrtPrc, unit;
	vector<int> IndOrd;
	int HlfWin = win/2, NObs = 2*HlfWin + 1, EndStg = HlfWin + 1, sz, szMat;
	double nrighe, price, nzside = floor(NObs*(delta/2)), UpBound = NObs - nzside, TrMean = 0, StDev = 0;
	double perc, avanza = 0;
	int subcount = 0;
	char ConfSym[32], ConfDate[32], *NomeOut;
	char **lineWords = new char*[9];
	for(int i = 0; i < 9; i++) lineWords[i] = new char[32];
	char *oraOut = new char[32];
	int CntFlds = 0, outcome, esito, Stage = 1, RowNum = 0;
	MyGzipDec *inFile;
	inFile = new MyGzipDec(NomeFile);
	NomeOut = new char[DirTemp.length() + 16];
	ogzstream out;
	string NomeTemp, useless, RowOut;
	for(int i = 0; i < nzside; i++) unit.push_back(0);
	for(int i = nzside; i < UpBound; i++) unit.push_back(1);
	for(int i = UpBound; i < NObs; i++) unit.push_back(0);
	
	double count_total = 0, count_notcorr_delayed = 0, count_BrownGallo = 0;
	vector<string> log_symbol, log_date;
	vector<double> day_total, day_notcorr_delayed, day_BrownGallo;
	
	NomeTemp = DirTemp + separator + "temp.taq";
	strcpy(ConfSym, "unifi");
	strcpy(ConfDate, "00000000");
	strcpy(NomeOut, DirTemp.c_str());
	strcat(NomeOut, separator);
	strcat(NomeOut, "inutile.taq");
	
	log_symbol.push_back("useless");
	log_date.push_back("useless");
	
	useless = NomeOut;
	nrighe = conteggio(NomeFile);
	out.open(NomeTemp.c_str());
	
	do
	{
		avanza++;
		subcount++;
		if(subcount == 1000)
		{
			perc = (avanza/nrighe)*100;
			if(perc <= 100)
			{
				printf("\r%4.2f completed.",perc);
				subcount = 0;
			}
			else
			{
				printf("\rTerminating the cleaning procedure...");
				subcount = 0;
			}
		}
		outcome = inFile->GetLineWords(lineWords, CntFlds, 9);
		
		if(CntFlds == 1 || !strcmp(lineWords[0], "SYMBOL")) continue;
		if(CntFlds != 8 && CntFlds != 9)
		{
			cout << "File anomaly detected at line " << avanza << ":" << endl;
			for(int j = 0; j < CntFlds; j++) cout << "\t" << lineWords[j];
			cout << endl;
			return 1;
		}
		
		count_total++;
		
		if(CntFlds == 8)
		{
			strcpy(lineWords[8], lineWords[7]);
			strcpy(lineWords[7], lineWords[6]);
			strcpy(lineWords[6], "@");
		}
		
		if(!strcmp(lineWords[5], "0") && strcmp(lineWords[6], "Z"))
		{
			price = atof(lineWords[3]);
		}
		else
		{
			count_notcorr_delayed++;
			continue;
		}
		
		if(strcmp(lineWords[0], ConfSym) || strcmp(lineWords[1], ConfDate))
		{
			Stage = 1;
			RowNum = 0;
			sz = SrtPrc.size();
			
			if(sz)
			{
				bubsort(SrtPrc, IndOrd, sz);
				trstdev(SrtPrc, unit, sz, TrMean, StDev);
				szMat = MtrRows.size();
				for(int t = 0; t < szMat; t++)
				{
					if((abs(MtrPrc[t] - TrMean)) < (3*StDev + gra))
					{
						out << MtrRows[t];
					}
					else
					{
						count_BrownGallo++;
					}
				}
			}
			
			MtrRows.clear();
			MtrPrc.clear();
			SrtPrc.clear();
			IndOrd.clear();
			strcpy(ConfSym, lineWords[0]);
			strcpy(ConfDate, lineWords[1]);
			out.close();
			
			if(!VerificaDir(NomeOut))
			{
				if(remove(NomeOut))
				{
					cout << "FATAL ERROR: unable to overwrite the file " << NomeOut << "." << endl;
					return 1;
				}
			}
			
			esito = rename(NomeTemp.c_str(), NomeOut);
			if(esito != 0)
			{
				cout << "FATAL ERROR: unable to create the file " << NomeOut << "." << endl;
				return 1;
			}
			
			day_total.push_back(count_total);
			day_notcorr_delayed.push_back(count_notcorr_delayed);
			day_BrownGallo.push_back(count_BrownGallo);
			
			count_total = 0;
			count_notcorr_delayed = 0;
			count_BrownGallo = 0;
			
			delete[] NomeOut;
			NomeOut = new char[DirOut.length() + 64];
			strcpy(NomeOut, DirOut.c_str());
			strcat(NomeOut, separator);
			strcat(NomeOut, "Trade");
			if(VerificaDir(NomeOut))
			{
				esito = create_directory(NomeOut);
				if(esito)
				{
					cout << "FATAL ERROR: unable to create the folder " << NomeOut << "." << endl;
					return 1;
				}
			}
			strcat(NomeOut, separator);
			strcat(NomeOut, ConfSym);
			if(VerificaDir(NomeOut))
			{
				esito = create_directory(NomeOut);
				if(esito)
				{
					cout << "FATAL ERROR: unable to create the folder " << NomeOut << "." << endl;
					return 1;
				}
			}
			log_symbol.push_back(ConfSym);
			
			strcat(NomeOut, separator);
			strcat(NomeOut, ConfDate);
			strcat(NomeOut, ".txt.gz");
			
			log_date.push_back(ConfDate);
			
			out.open(NomeTemp.c_str());
			for(int h = 0; h < 9; h++) out <<  "\t" << header[h];
			out << endl;
		}
		
		if(Stage)
		{
			orario(lineWords[2], oraOut);
			strcpy(lineWords[2], oraOut);
			RowOut = "\n";
			for(int k = 0; k < 9; k++)
			{
				RowOut += "\t";
				RowOut += lineWords[k];
			}
			RowNum++;
			MtrRows.push_back(RowOut);
			MtrPrc.push_back(price);
			SrtPrc.push_back(price);
			IndOrd.push_back(RowNum);
			if(RowNum == NObs)
			{
				bubsort(SrtPrc, IndOrd, NObs);
				trstdev(SrtPrc, unit, NObs, TrMean, StDev);
				
				for(int t = 0; t < EndStg; t++)
				{
					if((abs(MtrPrc[0] - TrMean)) < (3*StDev + gra))
					{
						out << MtrRows[0];
					}
					else
					{
						count_BrownGallo++;
					}
					MtrRows.erase(MtrRows.begin());
					MtrPrc.erase(MtrPrc.begin());
				}
				Stage = 0;
			}
		}
		else
		{
			orario(lineWords[2], oraOut);
			strcpy(lineWords[2], oraOut);
			RowOut = "\n";
			for(int k = 0; k < 9; k++)
			{
				RowOut += "\t";
				RowOut += lineWords[k];
			}
			MtrRows.push_back(RowOut);
			MtrPrc.push_back(price);			
			resort(SrtPrc, IndOrd, price, NObs);
			trstdev(SrtPrc, unit, NObs, TrMean, StDev);
			if((abs(MtrPrc[0] - TrMean)) < (3*StDev + gra))
			{
				out << MtrRows[0];
			}
			else
			{
				count_BrownGallo++;
			}
			MtrRows.erase(MtrRows.begin());
			MtrPrc.erase(MtrPrc.begin());
		}
	}while(outcome);
	sz = SrtPrc.size();
	if(sz)
	{
		bubsort(SrtPrc, IndOrd, sz);
		trstdev(SrtPrc, unit, sz, TrMean, StDev);
		szMat = MtrRows.size();
		for(int t = 0; t < szMat; t++)
		{
			if((abs(MtrPrc[t] - TrMean)) < (3*StDev + gra))
			{
				out << MtrRows[t];
			}
			else
			{
				count_BrownGallo++;
			}
		}
	}
	printf("\rTerminating the cleaning procedure...");
	cout << "\nCleaning procedure completed." << endl;
	MtrRows.clear();
	MtrPrc.clear();
	SrtPrc.clear();
	IndOrd.clear();
	out.close();
	
	if(!VerificaDir(NomeOut))
	{
		if(remove(NomeOut))
		{
			cout << "FATAL ERROR: unable to overwrite the file " << NomeOut << "." << endl;
			return 1;
		}
	}
	
	esito = rename(NomeTemp.c_str(), NomeOut);
	if(esito != 0)
	{
		cout << "\nFATAL ERROR: unable to create the file " << NomeOut << "." << endl;
		return 1;
	}
	
	day_total.push_back(count_total);
	day_notcorr_delayed.push_back(count_notcorr_delayed);
	day_BrownGallo.push_back(count_BrownGallo);
	
	if(create_log(DirOut, log_symbol, log_date, day_total, day_notcorr_delayed, day_BrownGallo, flag, deleted_log))
	{
		cout << "WARNING: log file creation failed." << endl;
	}
	
	delete inFile;
	delete[] NomeOut;
	for(int d = 0; d < 9; d++) delete[] lineWords[d];
	delete[] lineWords;
	delete[] oraOut;
	if(remove(useless.c_str()))
	{
		cout << "FATAL ERROR: unable to remove the file " << useless << "." << endl;
		return 1;
	}
	return 0;
}
