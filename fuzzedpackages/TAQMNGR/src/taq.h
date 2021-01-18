
#ifndef __TAQHEADER_H
#define __TAQHEADER_H

#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <gzstream.h>
#include <zlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <Rcpp.h>

#ifdef _WIN32
#undef Realloc
#undef Free
#include <windows.h>
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

class MyGzipDecMillisecond
{
	public:
	gzFile ComFile;
	int FieldsWidth[13];
	int FWsum;
	
	MyGzipDecMillisecond(const char *FileName);
	~MyGzipDecMillisecond();
	
	int GetLineWords(char **BufLineWords, int &Nfld);
	int isHeader(char *row);
	void InitFieldsWidth(char *headerRow);
	char *trimwhitespace(char *str);
};

double conteggio(const char *NomeFile);
int CleanTrade(const char *NomeFile, string DirOut, string DirTemp, int win, double delta, double gra, int flag, vector<string> &deleted_log);
int CleanTradeMS(const char *NomeFile, string DirOut, string DirTemp, int win, double delta, double gra, int flag, vector<string> &deleted_log);
int VerificaDir(const char *ndir);
void bubsort(vector<double> &nums, vector<int> &ordr, int size);
void bubsort_cleaning_report(vector<int> &date, vector<double> &total, vector<double> &notcorrected_delayed, vector<double> BrownGallo, int N);
void resort(vector<double> &nums, vector<int> &ordr, double prezzo, int size);
void trstdev(vector<double> vec, vector<double> ones, int size, double &trMean, double &StdDev);
void orario(char *Tempo, char *StrOut);
void itoa(int a, char *ChOut);
int ListaConfFile(string dir, string ListofCleaned, vector<string> &solofileConf, vector<string> &dirfileConf);
int ListaFile(string dir, vector<string> &solofile, vector<string> &dirfile);
int EmpEraseDir(string DirPath);
int IsQuote(const char *NomeFile);
int create_log(string dirOut, vector<string> sym, vector<string> date, vector<double> total, vector<double> notcorr_delayed, vector<double> BrownGallo, int cleaned_list_flag, vector<string> &allready_deleted);
int check_deleted(string file, vector<string> deleted);
int create_directory(const char *name);
void hide_file(const char *name);

#endif
