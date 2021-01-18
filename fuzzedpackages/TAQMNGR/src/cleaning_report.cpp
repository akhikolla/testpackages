#include <taq.h>
#include <sstream>

RcppExport SEXP CleaningReport(SEXP DirIn, SEXP Symbol)
{
	vector<int> day_date;
	vector<double> day_total, day_notcorrected_delayed, day_BrownGallo;
	int date;
	double total, notcorrected_delayed, BrownGallo;
	
	Rcpp::NumericVector ret(1);
	ret[0] = 0;
	
	string log_file, symbol, directory, directory_symbol, row;
	stringstream ss;
	
	directory = Rcpp::as<string>(DirIn);
	
	if(VerificaDir(directory.c_str()))
	{
		cout << "FATAL ERROR: folder " << directory << " not found." << endl;
		return ret;
	}
	
	symbol = Rcpp::as<string>(Symbol);
	directory_symbol = directory + separator + "Trade" + separator + symbol;
	
	if(VerificaDir(directory_symbol.c_str()))
	{
		cout << "FATAL ERROR: folder " << directory_symbol << " not found." << endl;
		return ret;
	}
	
	log_file = directory + separator + "Trade" + separator + symbol + separator + ".daily.log";
	if(VerificaDir(log_file.c_str()))
	{
		cout << "FATAL ERROR: no cleaning report available for the folder " << directory_symbol << ".\n";
		return ret;
	}
	
	ifstream input(log_file.c_str());
	while(input.good())
	{
		getline(input, row, '\n');
		if(row.empty()) continue;
		ss << row;
		ss >> date >> total >> notcorrected_delayed >> BrownGallo;
		
		ss.clear();
		ss.str(string());
		
		day_date.push_back(date);
		day_total.push_back(total);
		day_notcorrected_delayed.push_back(notcorrected_delayed);
		day_BrownGallo.push_back(BrownGallo);
	}
	
	int dim = day_date.size();
	double tot_total = 0, tot_notcorrected_delayed = 0, tot_BrownGallo = 0;
	bubsort_cleaning_report(day_date, day_total, day_notcorrected_delayed, day_BrownGallo, dim);
	
	cout << "#################################" << endl;
	cout << "#     DAILY CLEANING REPORT     #" << endl;
	cout << "#################################" << endl << endl;
	
	cout << "Directory: " << directory << endl;
	cout << "Symbol: " << symbol << endl << endl;
	
	cout << "DATE \t #TRADES \t NOTCORR_DELAY \t BROWN_GALLO" << endl;
	for(int j = 0; j < dim; j++)
	{
		cout << day_date[j] << "\t" << day_total[j] << "\t" << day_notcorrected_delayed[j] << "\t" << day_BrownGallo[j] << endl;
		tot_total += day_total[j];
		tot_notcorrected_delayed += day_notcorrected_delayed[j];
		tot_BrownGallo += day_BrownGallo[j];
	}
	
	cout << endl << "#################################" << endl;
	cout << "#     TOTAL CLEANING REPORT     #" << endl;
	cout << "#################################" << endl << endl;
	
	cout << "Directory: " << directory << endl;
	cout << "Symbol: " << symbol << endl << endl;
	
  // (2016-12-09) Replaced \% with %% to fix compilation warnings on some platforms
	cout << "#TRADES \t NOTCORR_DELAY \t BROWN_GALLO \t NOTCORR_DELAY(%%) \t BROWN_GALLO(%%)" << endl;
	cout << tot_total << "\t" << tot_notcorrected_delayed << "\t" << tot_BrownGallo << "\t" << tot_notcorrected_delayed / tot_total * 100 << "\t" << tot_BrownGallo / tot_total * 100 << endl << endl;
	return ret;
}
