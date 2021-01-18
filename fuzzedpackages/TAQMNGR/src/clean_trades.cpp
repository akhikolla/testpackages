#include <taq.h>

RcppExport SEXP CleanTickByTick(SEXP DirIn, SEXP DirOut, SEXP Window, SEXP DeltaTrimmed, SEXP Granularity, SEXP UseCleaned)
{
	string DirectoryIn = Rcpp::as<string>(DirIn);
	string DirectoryOut = Rcpp::as<string>(DirOut);
	int Win = Rcpp::as<int>(Window);
	double DelTrim = Rcpp::as<double>(DeltaTrimmed);
	double Granul = Rcpp::as<double>(Granularity);
	int UseDefault = Rcpp::as<int>(UseCleaned);
	
	Rcpp::NumericVector ret(1);
	ret[0] = 0;

	int nFile, tipo, esito;
	string scelta, sceltac, CleanedList = DirectoryOut + separator + ".cleaned.taq";
	string DirectoryTemp, EraseTemp, ComFile;
	vector<string> Files, SoloFiles;
	

  if(VerificaDir(DirectoryOut.c_str()))
  {
	  esito = create_directory(DirectoryOut.c_str());
	  if(esito)
	  {
		  Rcpp::Rcout << "FATAL ERROR: unable to create " << DirectoryOut << "." << endl;
		  return ret;
	  }
  }
  
  DirectoryTemp = DirectoryOut + separator + ".TempTaqGarbage";
  if(VerificaDir(DirectoryTemp.c_str()))
  {
	  esito = create_directory(DirectoryTemp.c_str());
	  if(esito)
	  {
		  Rcpp::Rcout << "FATAL ERROR: unable to create " << DirectoryTemp << "." << endl;
		  return ret;
	  }
	  hide_file(DirectoryTemp.c_str());
  }
  else
  {
	  if(EmpEraseDir(DirectoryTemp) != 0) return ret;
	  esito = create_directory(DirectoryTemp.c_str());
	  if(esito)
	  {
		  Rcpp::Rcout << "FATAL ERROR: unable to create " << DirectoryTemp << "." << endl;
		  return ret;
	  }
	  hide_file(DirectoryTemp.c_str());
  }
  
  if(UseDefault)
  {
	  if(!VerificaDir(CleanedList.c_str())) nFile = ListaConfFile(DirectoryIn, CleanedList, SoloFiles, Files);
	  else
	  {
		  Rcpp::Rcout << "WARNING: no cleaned file list found." << endl;
		  nFile = ListaFile(DirectoryIn, SoloFiles, Files);
	  }
  }
  else
  {
	  if(!VerificaDir(CleanedList.c_str()))
	  {
		  if(remove(CleanedList.c_str()) != 0)
		  {
			  Rcpp::Rcout << "FATAL ERROR: unable to remove '.cleaned.taq' file." << endl;
			  return ret;
		  }
	  }
	  nFile = ListaFile(DirectoryIn, SoloFiles, Files);
  }
  

  if(nFile == -1) return ret;
  else if(nFile == 0)
  {
	  Rcpp::Rcout << "The folder doesn't contain files to be cleaned." << endl << endl;
	  return ret;
  }

  Rcpp::Rcout << "The folder contains " << nFile << " files to be cleaned.\n";
  vector<string> deleted_files_log;
  ofstream ClOut(CleanedList.c_str(), ios::app);
  hide_file(CleanedList.c_str());

  for(int i=0; i < nFile; i++)
  {
	  ComFile = Files[i];
	  Rcpp::Rcout << endl << "Processing file " << i+1 << " of " << nFile << ".\n";
	  tipo = IsQuote(ComFile.c_str());
	  if(tipo == 0)
	  {
		  Rcpp::Rcout << ComFile << " is a quote report.\n";
		  Rcpp::Rcout << "No available cleaning procedure for quotes (yet).\n";
		  ClOut << SoloFiles[i] << endl;
	  }
	  else if(tipo == 1)
	  {
		  Rcpp::Rcout << Files[i] << " is a trade report.\n";
		  Rcpp::Rcout << "Cleaning trade...\n";
		  
		  esito = CleanTrade(ComFile.c_str(), DirectoryOut, DirectoryTemp, Win, DelTrim, Granul, UseDefault, deleted_files_log);
		  if(esito) return ret;
		  ClOut << SoloFiles[i] << endl;
	  }
	  else if(tipo == 5)
	  {
		  Rcpp::Rcout << Files[i] << " is a trade report (millisecond version).\n";
		  Rcpp::Rcout << "Cleaning trade (millisecond version)...\n";
		  
		  esito = CleanTradeMS(ComFile.c_str(), DirectoryOut, DirectoryTemp, Win, DelTrim, Granul, UseDefault, deleted_files_log);
		  if(esito) return ret;
		  ClOut << SoloFiles[i] << endl;
	  }
	  else
	  {
		  Rcpp::Rcout << "FATAL ERROR: unable to identify the file " << Files[i] << endl;
		  Rcpp::Rcout << "Please make sure you downloaded all the fields." << endl << endl;
		  return ret;
	  }
  }
  ClOut.close();
  EmpEraseDir(DirectoryTemp);
  return ret;
}
