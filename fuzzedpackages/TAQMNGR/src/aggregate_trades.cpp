#include <taq_aggregate.h>

RcppExport SEXP Aggregate(SEXP DirIn, SEXP Symbol, SEXP IntervalSeconds, SEXP useAggregated)
{
	string DirCleaned, DirectoryTemp, symbol, DirSymT, DirSymQ, secondi;
	vector<string> DownTrade, DownTradeD;
	int esito, t, nFile;
	int useAgg = Rcpp::as<int>(useAggregated);
	
	Rcpp::NumericVector ret(1);
	ret[0] = 0;
	
	DirCleaned = Rcpp::as<string>(DirIn);
	if(VerificaDir(DirCleaned.c_str()))
	{
		cout << "FATAL ERROR: folder " << DirCleaned << " not found." << endl;
		return ret;
	}
	
	symbol = Rcpp::as<string>(Symbol);
	DirSymT = DirCleaned + separator + "Trade" + separator + symbol;
	if(VerificaDir(DirSymT.c_str()))
	{
		cout << "FATAL ERROR: folder " << DirSymT << " not found." << endl;
		cout << "Data must be cleaned before being aggregated." << endl;
		return ret;
	}
	if(VerificaDir((DirCleaned + separator + "Aggregate").c_str()))
	{
		esito = create_directory((DirCleaned + separator + "Aggregate").c_str());
		if(esito)
		{
			cout << "FATAL ERROR: unable to create " << (DirCleaned + separator + "Aggregate") << "." << endl;
			return ret;
		}
	}
	DirectoryTemp = DirCleaned + separator + "Aggregate" + separator + ".TempTaqGarbage";
	if(VerificaDir(DirectoryTemp.c_str()))
	{
		esito = create_directory(DirectoryTemp.c_str());
		if(esito)
		{
			cout << "FATAL ERROR: unable to create " << DirectoryTemp << "." << endl;
			return ret;
		}
		hide_file(DirectoryTemp.c_str());
	}
	else
	{
		if(EmpEraseDir(DirectoryTemp) != 0)
		{
			cout << "FATAL ERROR: unable to delete " << DirectoryTemp << "." << endl;
			return ret;
		}
		esito = create_directory(DirectoryTemp.c_str());
		if(esito)
		{
			cout << "FATAL ERROR: unable to create " << DirectoryTemp << "." << endl;
			return ret;
		}
		hide_file(DirectoryTemp.c_str());
	}
	
	secondi = Rcpp::as<string>(IntervalSeconds);
	
	cout << endl << "Aggregating " << symbol << " data:" << endl;
	
	if(ListaFileToBeAggregated(DirCleaned, DirSymT, symbol, secondi, DownTrade, DownTradeD, useAgg) == -1)
	{
		EmpEraseDir(DirectoryTemp);
		return ret;
	}
	nFile = DownTrade.size();
	
	if(nFile == 0)
	{
		cout << "The folder doesn't contain files to be aggregated." << endl;
		EmpEraseDir(DirectoryTemp);
		return ret;
	}
	double nF = (double) nFile, avanza = 0, perc;
	
	for(t = 0; t < nFile; t++)
	{
		avanza++;
		perc = (avanza/nF)*100;
		printf("\r%4.2f completed.",perc);
		if(AggregateT(DirCleaned, symbol, DownTradeD[t], DownTrade[t], secondi)) return 0;
	}
	cout << "\rProcedure completed." << endl;
	EmpEraseDir(DirectoryTemp);
	return ret;
}
