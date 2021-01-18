#ifndef __TAQHEADER_IMPORT_H
#define __TAQHEADER_IMPORT_H

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
#include <gzstream.h>
#include <time.h>
#include <Rcpp.h>

#ifdef _WIN32
#define separator "\\"
#else
#define separator "/"
#endif

#define cout Rcpp::Rcout

using namespace std;

int ListaToImport(string Directory, int StrYear, int StrMonth, int StrDay, int EndYear, int EndMonth, int EndDay, int &nToImport, vector<string> &ToImport, int &nMissing, vector<string> &Missing);
int ControlInput(string DirIn, string symbol, string seconds, int StYear, int StMonth, int StDay, int EYear, int EMonth, int EDay, int MissCont, vector<string> &ToImp);
int ListaFile(string dir, vector<string> &solofile, vector<string> &dirfile);
int NeededList(int StaY, int StaM, int StaD, int EndY, int EndM, int EndD, vector<string> &NeedFile);
int VerificaDir(const char *ndir);
void BubsortDelDouble(vector<int> &nums, vector<int> &numsOut);
#endif
