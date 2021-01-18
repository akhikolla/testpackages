
#ifndef __TAQHEADER_AGGREGATE_H
#define __TAQHEADER_AGGREGATE_H

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include "gzstream.h"
#include <Rcpp.h>

#ifdef _WIN32
#define separator "\\"
#else
#define separator "/"
#endif

#define cout Rcpp::Rcout
#define printf Rprintf

using namespace std;

class MyGzipDec
{
	public:
	gzFile ComFile;
	
	MyGzipDec(const char *FileName);
	~MyGzipDec();
	
	int GetLineWords(char **BufLineWords, int &Nfld, int dim);
};

int AggregateT(string DirClean, string sym, string TradeD, string TradeF, string Secs);
int ListaFileToBeAggregated(string dir_root, string dir_cleaned, string symbol, string seconds, vector<string> &solofile, vector<string> &dirfile, int useAggregated);
int ListaFile(string dir, vector<string> &solofile, vector<string> &dirfile);
int VerificaDir(const char *ndir);
int EmpEraseDir(string DirPath);
int TimeStamp(int interval, vector<int> &tempo);
int create_directory(const char *name);
void hide_file(const char *name);

#endif
