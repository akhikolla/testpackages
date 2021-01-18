#include <taq_import.h>

RcppExport SEXP Import(SEXP directory, SEXP symbol, SEXP Campi, SEXP StartYear, SEXP StartMonth, SEXP StartDay, SEXP EndYear, SEXP EndMonth, SEXP EndDay, SEXP InterSeconds)
{	
	string DirIn, Symbol;
	int StartingYear, StartingMonth, StartingDay, EndingYear, EndingMonth, EndingDay, IntervalSeconds, MissC = 1;
	vector<string> FieldString;
	
	DirIn = Rcpp::as<string>(directory);
	Symbol = Rcpp::as<string>(symbol);
	StartingYear = Rcpp::as<int>(StartYear);
	StartingMonth = Rcpp::as<int>(StartMonth);
	StartingDay = Rcpp::as<int>(StartDay);
	EndingYear = Rcpp::as<int>(EndYear);
	EndingMonth = Rcpp::as<int>(EndMonth);
	EndingDay = Rcpp::as<int>(EndDay);
	IntervalSeconds = Rcpp::as<int>(InterSeconds);
	
	Rcpp::NumericVector ret(1);
	ret[0] = 0;
	
	Rcpp::CharacterVector FieldChar(Campi);
	int nCampi = FieldChar.size();
	for(int h = 0; h < nCampi; h++) FieldString.push_back(string(FieldChar[h]));
	
	string Sec;
	stringstream outstr;
	outstr << IntervalSeconds;
	Sec = outstr.str();
	
	vector<string> Fields;
	Fields.push_back("FIRST");
	Fields.push_back("MIN");
	Fields.push_back("MAX");
	Fields.push_back("LAST");
	Fields.push_back("SIZE");
	Fields.push_back("#TRADES");
	Fields.push_back("VWAP");

	vector<string> ToImpList;
	int nToImp = ControlInput(DirIn, Symbol, Sec, StartingYear, StartingMonth, StartingDay, EndingYear, EndingMonth, EndingDay, MissC, ToImpList);
	
	if(nToImp == 0) return ret;
	vector<int> FieldsIn, FieldsReady;
	
	for(int j = 0; j < nCampi; j++)
	{
		for(int h = 0; h < 7; h++)
		{
			if(!FieldString[j].compare(Fields[h]))
			{
				FieldsIn.push_back(h+2);
				break;
			}
		}
	}
	BubsortDelDouble(FieldsIn, FieldsReady);
	int contaRighe = 0;
	igzstream in;
	string row;
	
	for(int f = 0; f < nToImp; f++)
	{
		in.open((ToImpList[f]).c_str());
		while(in.good())
		{
			getline(in, row, '\n');
			if(!row.empty()) contaRighe++;
		}
		in.close();
		in.clear();
	}
	
	int IndCol, r = 0, c = 0, nRighe = contaRighe - nToImp, nFieldsReady = FieldsReady.size(), nColonne = nFieldsReady + 2;
	Rcpp::NumericMatrix dataOut(nRighe, nColonne);
	string riga, SingleRow[9];
	stringstream ss;
	igzstream dataIn;

	for(int f = 0; f < nToImp; f++)
	{
		dataIn.open((ToImpList[f]).c_str());
		do
		{
			getline(dataIn, riga,'\n');
		}while(riga.empty());
		while(dataIn.good())
		{
			getline(dataIn, riga, '\n');
			if(riga.empty()) continue;
			
			ss << riga;
			for(int k = 0; k < 9; k++) ss >> SingleRow[k];
			ss.clear();
			ss.str(string());
			
			dataOut(r, c) = atof((SingleRow[0]).c_str());
			c++;
			dataOut(r, c) = atof((SingleRow[1]).c_str());
			for(int g = 0; g < nFieldsReady; g++)
			{
				IndCol = FieldsReady[g];
				c++;
				if(!(SingleRow[IndCol]).compare("NaN")) dataOut(r, c) = NA_REAL;
				else dataOut(r, c) = atof((SingleRow[IndCol]).c_str());
			}
			r++;
			c = 0;
		}
		dataIn.close();
		dataIn.clear();
	}
	Rcpp::CharacterVector FieldsOut(nColonne);
	FieldsOut[0] = "DATE";
	FieldsOut[1] = "TIME";
	
	for(int g = 0; g < nFieldsReady; g++)
	{
		IndCol = FieldsReady[g];
		FieldsOut[g+2] = (Fields[IndCol - 2]).c_str();
	}
	Rcpp::CharacterVector rows(nRighe);
	Rcpp::List dimNames = Rcpp::List::create(rows, FieldsOut);
	dataOut.attr("dimnames") = dimNames;
	return dataOut;
}
