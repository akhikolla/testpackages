//***********************************************************************

//
//
//                              COMMON.CPP
//
//
// Author(s) : Serge GUILLAUME, Brigitte CHARNOMORDIC 
// FISPRO Version : 3.2 - Copyright INRA - Cemagref - 2002 - IDDN.FR.001.030024.000.R.P.2005.003.31235
// Licence http://www.cecill.info/licences/Licence_CeCILL_V2-en.html/
// Last modification date:  December 1, 2007
// File : Out of class functions used by FISPRO, part of library fispro

//**********************************************************************



//*********************************************************************
//
//
//                Common functions for FIS classes
//
//
//*********************************************************************


#include "common.h"
#include <math.h>  
#include <stdexcept>
#include <stdlib.h>

char ErrorMsg[300];
char DataFilePath[DATAFILEMAXPATH] = "";
char ** VarNameG = NULL;
int NbVarG = 0;


int MaxLineSize( ifstream &f )
  //**************************
{
  int max_size = 0, size = 1;
  f.seekg( 0, ios::end );
  long endal_position = f.tellg();
  for( long i=0 ; i<=endal_position ; i++ )
    {
      f.seekg( i, ios::beg );
      char car = (char)f.peek();
      if( car == '\n' )
	{
	  max_size = ( max_size > size ) ? max_size : size;
	  size = 1;
	}
      else size++;
    }
  f.seekg( 0, ios::beg );
  f.clear();
  
  return max_size;
}

int MaxLineSize( const char *FileName )
  //*****************************
{
  int max_size = 0, size = 1;
  char car;
  FILE *f = fopen( FileName, "rt" );
  while( (car = (int)fgetc(f)) != EOF )
    {
    if( car == '\n' )
      {
      max_size = ( max_size > size ) ? max_size : size;
      size = 1;
      }
    else size++;
    }
  fclose(f);
  return max_size;
}

//******************************************************************
//
//           Random generation (Park and Miller version)
//
//******************************************************************
#ifdef _OPENMP
#include <omp.h>
static long int seed = omp_get_thread_num();
#pragma omp threadprivate(seed)
#else
static long int seed = 1;
#endif

int setseed(long int newseed)
  //*************************
{
  //if ((newseed < 1) || (newseed > (long int) FIS_RAND_MAX)) return 1;
  if (newseed == 0)
    newseed=time(NULL);
  if ((newseed < 0) || (newseed > (long int) FIS_RAND_MAX)) return 1;
  
  seed = newseed; 
 
  return 0;
}

double randpm(void)
  //***************
{
  long int a=16807,  q=127773, r=2836;
  int lo, hi, test;
  
  hi = seed / q;
  lo = seed % q;
  test = a * lo - r * hi;
  seed = (test > 0) ? test : test + FIS_RAND_MAX;
  return ((double) seed / FIS_RAND_MAX);
}

double FisRand(void)      
  //****************
{
  return randpm();
}

double GaussDice( double sig, int n )
  //*********************************
{
  // returns pseudo random normal element 
  // normal distribution mean 0 and range n*sig
  // default n=12 - default range of result: between -n*sig and n*sig
  int i;
  double s;
  if(!n) return 0.0;
  s = 0.0;
  for ( i = 0; i < n ; i++ ) s += FisRand();//s varies between 0 and n
  s = s - (double) n/2;// s varies between -n/2 and n/2
  return ( s * sig ) ;//to limit range of variations
}

//******************************************************************
//
//           Missing value Handling
//
//******************************************************************
double FisMknan(void)     
  //*****************      
{
  return sqrt(-1.);   
}

int FisIsnan(double x)    
  //******************
{
  return (x != x);        
}

//******************************************************************
//
//           Sort functions
//
//******************************************************************

int CmpDbl(const void * a, const void * b)
  //**************************************
{
  if( *(double *)a > *(double *)b ) return -1;
  if( *(double *)a < *(double *)b ) return 1;
  return 0;
}

int CmpDblAsc(const void * a, const void * b)
  //*****************************************
{
  if( (FisIsnan(*(double *)a)) && (FisIsnan(*(double *)b)) ) return 0;
  if(FisIsnan(*(double *)a)) return 1;
  if(FisIsnan(*(double *)b)) return -1;

  if( *(double *)a > *(double *)b ) return 1;
  if( *(double *)a < *(double *)b ) return -1;
  return 0;
}

void StatArray(double * T, int Tsize, int nb, double & median, double & mean, 
	       double & delta, double & min, double & max, int estim)
  //*************************************************************************
{
  double tmp; 
  int i;

  int dof;  // degrees of freedom 

  
  if(Tsize - 2*nb < 1)  // not enough items
    {  median = FisMknan(); mean = median; delta = median;   return; }


  qsort(T, Tsize, sizeof(double), CmpDbl); 
  
  min = T[Tsize - 1];
  max = T[0];
  median = T[Tsize/2];
  
  tmp = 0.0;
  for(i = nb; i < Tsize - nb; i++)  tmp += T[i]; 
  mean = tmp / (Tsize-2*nb);

  tmp = 0.0;
  dof = Tsize- 2*nb;
  
  if (estim) dof --;
  
  for(i = nb; i < Tsize - nb; i++)  tmp += (T[i] - mean) * (T[i] - mean);
  if(dof) delta = sqrt(tmp / dof);  
  else delta = 0.0;
}

void StatArrayQuart(double * T, int Tsize, int nb, double & firstq, double & median, double & thirdq, double & min, double & max)
  //*************************************************************************
{
  if(Tsize - 2*nb < 1)  // not enough items
    {  median = FisMknan(); firstq= median; thirdq = median;   return; }


  qsort(T, Tsize, sizeof(double), CmpDbl); 
  //returns sorted T

  min = T[Tsize - 1];
  max = T[0];
  median = T[(int)Tsize/2];
  thirdq = T[(int)(Tsize*0.75)];
  firstq = T[(int)Tsize/4];
}

void InitUniq(double *T, int n, double * & ValPos, int & NPos)
  //**********************************************************
{
#define IMPOSSIBLE (-INFINI-0.0005) 
  
  int j, k;
  double *val;
 
  NPos = 0;
  if(n <= 0) return; 
  val = new double [n];

  for(j = 0; j < n; j ++) val[j] = IMPOSSIBLE;

  for(j = 0; j < n; j ++)
    {
      for(k = 0; k < NPos; k++)
	{
	  if(fabs(T[j]- val[k])< EPSILON)  break;
	}
     
      if(k == NPos)
	{
	  val[NPos] = T[j];
	  NPos++;
	}
    }	
  ValPos = new double [NPos];
  for(k = 0; k < NPos; k++) 	ValPos[k] = val[k];
 
  delete [] val;
}


int SortUniq(double * T, int Tsize, double ** U, int ** p, int &nu, double thres)
  //*****************************************************************************
{
  double * V; 
  int *pV, *ind;
  int i, j;
  int stop, ret;
  double *sum;

  V = new double [Tsize];
  pV = new int [Tsize];
  ind = new int [Tsize]; 
  sum = new double [Tsize];

  ret = 0;
 
  for(i = 0; i < Tsize; i++) V[i] = T[i];
 
  qsort(V, Tsize, sizeof(double), CmpDblAsc); 

  j = 0;
  pV[0] = 1;      // occurrence number
  ind[0] = 0;     // unique value index
  stop = false;
  sum[0] = V[0];

  for(i = 1; i < Tsize && stop == false; i++) 
    {
      if (V[i] <= V[ind[j]] + thres)  
	// One occurrence more 
	{
	  pV[j]++;  
	  sum[j]+=V[i];
	}
      else // New unique value
	{
	  sum[j]/=pV[j];
	  j++; 
	  pV[j] = 1;
	  ind[j] = i;
	  sum[j]=V[i];
	  if(FisIsnan(V[i]))  // All the next are NaN
	    {
	      pV[j] = Tsize - i;
	      ret = -1;
	      stop = true;
	    }
	}
    }
 
  sum[j]/=pV[j];
 
  nu = j + 1; // Number of unique  values
  (*U) = new double [nu];
  (*p) = new int [nu];

  for(i = 0; i < nu; i++)
    {
      (*U)[i] = sum[i]; 
      (*p)[i] = pV[i];
    }
 
  delete [] pV;
  delete [] V;
  delete [] ind;
  delete [] sum;
  return ret;
} // End of SortUniq

//******************************************************************
//
//           Search functions
//
//******************************************************************

int SearchStr(const char * source, char * chaine, char Delim)
  //***************************************************
{
  chaine[0] = 0;
  const char * begin  = strchr(source, Delim);
  if(begin == NULL) return 1;
  int IndexDebut = begin - source;
  const char * end = strchr(source + IndexDebut + 1, Delim); 
  if(end  == NULL) return 1;
  int IndexFin = end - source;
  strncat(chaine, source + IndexDebut + 1, IndexFin - IndexDebut - 1);

  return 0;
}

int SearchNb(const char * source, double * val, int size, 
	     char separator, int DelimBeg, int DelimEnd)
  //****************************************************
{
  char *tmp;
  int endi, i, sourcelen;
  int IndexDebut=0, IndexFin=0;
  const char *end;
  const char *begin;
  double tmp_val;
  char c[5];

  tmp = new char[strlen(source)+1];

  if(DelimBeg == DELIM_START) IndexDebut = 0;
  else
    {
      begin  = strchr(source, DelimBeg);
      if(begin == NULL) return -1;
      IndexDebut = begin - source + 1;
    }

  end = strchr(source + IndexDebut + 1, DelimEnd);
  endi = end - source;
  sourcelen=strlen(source);
  i = 0;

  while (1)
    {
      if (IndexDebut>=sourcelen) break;
      end = strchr(source + IndexDebut + 1, separator);
      if(end != NULL)
	IndexFin = end - source;

      // the next separator does not belong to source
      if(end != NULL && (end - source > endi))
	{
	  delete [] tmp;
	  return i;
	}
      else if (end  == NULL) // no further separator
	{
	  IndexFin = endi;
	  while((source[IndexDebut] == ' ' || source[IndexDebut] == '\t'|| source[IndexDebut] == '\r') 
		&& IndexDebut < IndexFin) IndexDebut++;
	  if((IndexFin - IndexDebut) < 1) // empty string
	    {
	      delete [] tmp;
	      return i;
	    }
	}
      tmp[0] = 0;
      strncat(tmp, source + IndexDebut, IndexFin - IndexDebut);

      if(strstr(tmp, CODE_NAN)) val[i++] = FisMknan();          
      else
	{
	  if(sscanf( tmp, "%lf %4s", &tmp_val, c ) != 1)
	    {
	      sprintf(ErrorMsg, "~NotaNumber~:  %.50s", tmp);
	      throw std::runtime_error( ErrorMsg );
	    }
	  else val[i++] = tmp_val;    
	}
      if(i == size) endi = 1;
      IndexDebut = IndexFin + 1;
    }
  delete [] tmp;
  return i;
}// End of SearchNb(double *)

int SearchVarNames(char * source, int size, char sep)
  //*************************************************
{
  unsigned int j,len;
  int begin, end, IndexBegin;

  len = strlen(source);
  j = 0;
  IndexBegin = 0;
  NbVarG = 0;

  VarNameG = new char *[size]; 
  begin = end = false;
  for(j = 0; j < len+1; j++)
    {
      if((begin && !end && !isalnum(source[j]) && (source[j] != '_') ) ||
	 (begin && j == len && (isalnum(source[j]) || (source[j] == '_')) )) // end of line
	{ 
	  end = true; 
	  VarNameG[NbVarG] = new char [j - IndexBegin + 1];
	  VarNameG[NbVarG][0] = '\0';     
	  strncat(VarNameG[NbVarG], source + IndexBegin, j - IndexBegin);
	  NbVarG++;
	}
      else if(!begin && (isalnum(source[j]) || (source[j] == '_'))) 
	{ begin = true; IndexBegin = j; }      
      if(source[j] == sep)
	{ begin = false; end = false; }
      if(NbVarG == size) break;
    } 
  return NbVarG;
}

int FileNameIndex(const char * source)
  //****************************
{
#define UNIX_CASE    1
#define WINDOWS_CASE 2

  int ostype = UNIX_CASE;
  char d;

  const char * begin  = strchr(source, '/'); // unix case

  if(begin == NULL) 
    {
    begin = strchr(source, '\\');      // windows case
    if(begin != NULL) ostype = WINDOWS_CASE;
    else return 0;
    }  
  
  if(ostype == UNIX_CASE) d = '/';
  else d = '\\';

  const char *last = begin;
  while(last != NULL)
    {
      begin = last;
      last = strchr(begin + 1, d);  
    }

  return begin - source + 1;
}


//******************************************************************
//
//
//                      Sample file reading
//
//
//******************************************************************

double ** ReadSampleFile(const char *N,  int &NCol, int &NRow)
  //****************************************************
{
  int i, hdr, len;
  double ** T;
  char s;
  len = 0;

  s = ReadSeparator(N, hdr);
  SampleFileSize(N,  NCol, NRow, len, s, hdr);
  T = new double * [NRow];
  try
    {
      for(i = 0; i < NRow; i ++)    T[i] = NULL;
      for(i = 0; i < NRow; i ++)    T[i] = new double [NCol];
      ReadItems(N,  NCol, NRow, T, len, s, hdr);
    }
  catch(std::exception &e)
    {
      for(i = 0; i < NRow; i ++)  delete []  T[i];
      delete []  T;
      throw;
    }
  return T;
}

char ReadSeparator(const char *FileName, int & hdr)
  //*****************************************
{
  ifstream f(FileName);
  char *buf, ret;
  unsigned i, Msize;

  // correct bug with minus and plus values - August 30 2007
  if (f.fail())
    {
      sprintf( ErrorMsg, "~CannotOpenDataFile~: %.100s~", FileName);
      throw std::runtime_error( ErrorMsg );
    }

  Msize = MaxLineSize( FileName );
      
  buf = new char[Msize]; 
  hdr = false;
  i = 0;
  ret = 0;
  f.getline(buf, Msize);
  while ( isspace(buf[i]) ) i++; 
  if (!isdigit(buf[i]) && buf[i]!='-' && buf[i]!='+')   // header line 
    {    
      hdr = true;
      f.getline(buf, Msize);
      i = 0;
    }
  while ( isdigit(buf[i]) || isspace(buf[i]) || (buf[i]== '.')|| (buf[i]== '-')|| (buf[i]== '+') )   i++;
  if(i < strlen(buf)) ret = *(buf+i);
  else  ret = SEPARE;
  // When no valid field separator has been found, we use SEPARE (,) to read one-column files
  delete [] buf;
  return ret;
}

void SampleFileSize(const char *FileName, int &NbCol, int &NbRow, int & LineLen, char sep, int hdr)
  //*****************************************************************************************
{
  ifstream f(FileName);
  int NbN;
  char * buf;

  if (f.fail())
    {
      sprintf( ErrorMsg, "~CannotOpenDataFile~: %.100s~", FileName);
      throw std::runtime_error( ErrorMsg );
    }

  LineLen = MaxLineSize( FileName );
      
  buf = new char[LineLen]; 
  NbCol = 0;
  NbRow = 0;

  if(hdr)  f.getline(buf, LineLen);

  while (!f.eof())
    {
      f.getline(buf, LineLen);
      NbN = CntNbs(buf, sep);
      NbCol = (NbCol < NbN ? NbN : NbCol);
      if( strlen(buf) && (buf[0]!=0x0D) )  NbRow++;  // skip empty lines
    }
  delete [] buf;
}


int ReadOneItem(ifstream & flot, int bufsize, char sep, double * Ind, int N)
  //************************************************************************
{
  char * buf;
  int ret;

  ret = -1;
  buf = new char[bufsize]; 
  try
    {
      flot.getline(buf, bufsize);
      if((strlen(buf) != 0) && (buf[0]!=0x0D))
	ret = SearchNb(buf, Ind, N, sep);
  
      delete [] buf;
      return ret;
    }

  catch( std::exception &e )
    {
      sprintf( ErrorMsg, "~ErrorInDataFile~\n~ErrorInReadOneItem~:%.50s\n%.100s", buf, e.what() );
      delete [] buf;
      throw std::runtime_error( ErrorMsg );
    }
}



void ReadItems(const char *FileName, int NbVals, int NbInd, double ** vv, int bufsize, char sep, int hdr)
  //***********************************************************************************************
{
  int j;
  char * buf=NULL;
  ifstream f(FileName);

  j = 0;

  if (f.fail())
    {
      sprintf( ErrorMsg, "~CannotOpenDataFile~: %.100s~", FileName);
      throw std::runtime_error( ErrorMsg );
    }

  buf = new char[bufsize]; 
  try
    {
      if(hdr)
	{
	  if(VarNameG != NULL) 
	    {
	      for(j = 0; j < NbVarG; j++)
		if(VarNameG[j] != NULL)  delete [] VarNameG[j];
	      delete [] VarNameG;
	      VarNameG = NULL;
	    }
	  f.getline(buf, bufsize);
	  if(SearchVarNames(buf, NbVals, sep) != NbVals)
	    {
	      sprintf( ErrorMsg,"~ErrorInDataFile~: %.100s\n~UnexpectedNumberOfColumnsInLineOne ~", FileName);
	      //	      delete [] buf;
	      throw std::runtime_error(ErrorMsg);
	    }
  	}

      for(j = 0; j < NbInd; j++)
	{
	  f.getline(buf, bufsize);
	  if((strlen(buf) != 0) && (buf[0] != 0x0D))
	    if(SearchNb(buf, vv[j], NbVals, sep)  != NbVals )
	      {
		sprintf( ErrorMsg, "~ErrorInDataFile~: %.100s\n~UnexpectedNumberOfColumnsInLine~ %d~", FileName, j + 1);
		//delete [] buf;
		throw std::runtime_error( ErrorMsg );
	      }
	}
      delete [] buf;
    }
  catch( std::exception &e )
    {
      delete [] buf;
      sprintf( ErrorMsg, "~ErrorInDataFile~\n~ErrorInLine~: %d\n%.100s", j+1, e.what() );
      throw std::runtime_error( ErrorMsg );
    }
}// End of ReadItems



int CntNbs( char *buf, char separator, char begin, char end )
  //*********************************************************
{
  int nb_separator = 1; // the # of numbers is equal to the # of separators in the buffer + 1
  int ibegin=0, iend=strlen(buf); //  start and end indices for search
  // search for  start character in buffer
  if( begin != 0 )
    {
      for( ibegin=0 ; ibegin<(int)strlen(buf) ; ibegin++ )
	if( buf[ibegin] == begin )   break;
      //  exception
    }
  // search for  end character in buffer
  if( end != 0 )
    {
      for( iend=ibegin ; iend<(int)strlen(buf) ; iend++ )
	if( buf[iend] == end ) 	  break;
      //  exception
    }      
  // search for separators in buffer
  for( int i=ibegin ; i<iend ; i++ )
    if( buf[i] == separator )  nb_separator++;

  return nb_separator;
}

void GetColumn(double ** T, int n, int col, double *C)
  //**************************************************
{
  int i;

  for(i = 0; i < n; i++)  C[i] = T[i][col];
}

//******************************************************************
//
//           k-means clustering 
//
//******************************************************************

int Kmeans(double *T, int ni, double *C, int nc, int norm)
  //******************************************************
{
  double change;    // value to be compared to threshold
  double * tmp;     // centers being built
  int * nb;         // count in each group
  int k;            // number of iterations
  int groupe;       // group number to assign to item
  int i, j;

  if(norm)
    {
      double min, max, diff;
      min = max = T[0];

      for(i = 1; i < ni; i++)
	{
	  if(T[i] < min) min = T[i];
	  if(T[i] > max) max = T[i];
	}

      diff = max - min;
      for(i = 0; i < ni; i++) T[i] = (T[i] - min) / diff;
    }

  tmp = new double [nc];
  nb = new int [nc];

  change = CVG_THRES + 1.;
  k = 0;

  while(change > CVG_THRES)
    {
      k++;
      change = 0.;
      for(j = 0; j < nc; j++)
	{
	  tmp[j] = 0.;
	  nb[j] = 0;
	}

      for(i = 0; i < ni; i++)
	{
	  groupe = AssignClas(T[i], C, nc);
	  tmp[groupe] += T[i];
	  nb[groupe]++;
	}

      for(j = 0; j < nc; j++)
	{
	  if(nb[j])
	    { 
	      tmp[j] /= nb[j];
	      change += (tmp[j] - C[j])*(tmp[j] - C[j]);
	      C[j] = tmp[j];
	    }
	}
    }

  delete [] tmp;
  delete [] nb;
  return k;
}

int AssignClas(double v, double *cg, int ng)
  //****************************************
{
  int i, gr;
  double val, min;

  min = 1e20;
  gr = -1;

  for(i = 0; i < ng; i++)
    {
      val = (v-cg[i])*(v-cg[i]);
      if(val < min) 
	{
	  gr = i;
	  min = val;
	}
    }
  return gr;
}

void InitCentres(double *&C, int n, double min, double max)
  //*******************************************************
{
  double inc;
  int i;

  C = new double[n];
  inc = (max - min)/(n - 1);   
  for(i = 0; i < n ;i++)
      C[i] = min + i*inc;
}

int KmeansNE(double *T, int ni, double * C, int & nc)
  //*************************************************
{
  int i, j, empty, * eff;
  eff = new int [nc];
  empty = 0;

  for(i = 0; i < nc ; i++) eff[i] = 0;

  for(i = 0; i < ni; i++) eff[AssignClas(T[i],C,nc)]++;

  for(i = 0; i < nc - empty; i++) 
    if (!eff[i]) 
      {
	empty ++;
	for(j = i; j < nc - empty - 1; j++)
	  {
	    C[j] = C[j + 1];
	    eff[j] = eff[j + 1];
	    C[nc - empty] = INFINI;
	    eff[nc - empty] = 0;
	  }
      }
  delete [] eff;
  nc -= empty; // number of non empty groups
  return empty;
}


// Multidimensional version

int Kmeans(double **T, int ni, double **C, int nc, int dim, int norm)
  //*****************************************************************
{
  double change;    //  value to compare to threshold
  double ** tmp;    // centers being built
  int * nb;         // count in each group
  int k;            // number of iterations
  int groupe;       // group number to assign to item
  int i, j, d;

  if(norm)
    {
      double min, max, diff;
      for(d = 0; d < dim; d++)
	{
	  min = max = T[0][d];
	  
	  for(i = 1; i < ni; i++)
	    {
	      if(T[i][d] < min) min = T[i][d];
	      if(T[i][d] > max) max = T[i][d];
	    }

	  diff = max - min;
	  for(i = 0; i < ni; i++) T[i][d] = (T[i][d] - min) / diff;
	}
    }

  tmp = new double * [nc];
  for(j = 0; j < nc; j++) tmp[j] = new double [dim];



  nb = new int [nc];

  change = CVG_THRES + 1.;
  k = 0;

  while(change > CVG_THRES)
    {
      k++;
      change = 0.;
      for(j = 0; j < nc; j++)
	{
	  for(d = 0; d < dim; d++) tmp[j][d] = 0.;
	  nb[j] = 0;
	}

      for(i = 0; i < ni; i++)
	{
	  groupe = AssignClas(T[i], dim, C, nc);
	  for(d = 0; d < dim; d++) tmp[groupe][d] += T[i][d];
	  nb[groupe]++;
	}

      for(j = 0; j < nc; j++)
	{
	  if(nb[j])
	    { 
	      for(d = 0; d < dim; d++)
		{
		  tmp[j][d] /= nb[j];
		  change += (tmp[j][d] - C[j][d])*(tmp[j][d] - C[j][d]);
		  C[j][d] = tmp[j][d];
		}
	    }
	}
    }

  for(j = 0; j < nc; j++) delete [] tmp[j];
  delete [] tmp;
  delete [] nb;
  return k;
}

int AssignClas(double *v, int dim, double **cg, int ng)
  //***************************************************
{
  int i, d, gr;
  double val, min;

  min = 1e20;
  gr = -1;

  for(i = 0; i < ng; i++)
    {
      val = 0.;
      for(d = 0; d < dim; d++)
	val += (v[d]-cg[i][d])*(v[d]-cg[i][d]);

      if(val < min) 
	{
	  gr = i;
	  min = val;
	}
    }
  return gr;
}

int KmeansNE(double **T, int ni, double ** C, int & nc, int dim)
  //************************************************************
{
  int i, j, d, empty, * eff;
  eff = new int [nc];
  empty = 0;

  for(i = 0; i < nc ; i++) eff[i] = 0;

  for(i = 0; i < ni; i++) eff[AssignClas(T[i], dim, C, nc)]++;

  for(i = 0; i < nc - empty; i++) 
    if (!eff[i]) 
      {
	empty ++;
	for(j = i; j < nc - empty - 1; j++)
	  {
	    for(d = 0; d < dim; d++)
	      {
		C[j][d] = C[j + 1][d];
		C[nc - empty][d] = INFINI;
	      }
	    eff[j] = eff[j + 1];
	    eff[nc - empty] = 0;
	  }
      }
  delete [] eff;
  nc -= empty; // number of non empty groups
  return empty;
}

// Global data normalisation
void Normalize ( double ** SampleData, int col, int nbrow, double min, double max)
  //******************************************************************************
{
  double ratio = max - min;
  for ( int cpt = 0 ; cpt < nbrow ; cpt++ ) 
    SampleData[cpt][col]=(SampleData[cpt][col] - min)/ratio;
  return;
}

void UnNormalize ( double ** SampleData, int col, int nbrow, double min, double max)
  //********************************************************************************
{
  double ratio = max - min;
  for ( int cpt = 0 ; cpt < nbrow ; cpt++ ) 
    SampleData[cpt][col]=SampleData[cpt][col] * ratio +  min;
  return;
}

void ReadTemplate(char * file, double &KW, double &SW)
  //********************************************************************************
{
  double ** templ=NULL;
  int ncol,nrow;

  templ=ReadSampleFile(file,ncol,nrow);
  
  if (ncol !=2)
    {
      sprintf( ErrorMsg, "~#columns~must~be~equal~to~two");
      throw std::runtime_error( ErrorMsg );
    }
  if (nrow<1)
    {
      sprintf( ErrorMsg, "no~rows~in~template~file");
      throw std::runtime_error( ErrorMsg );
    }
  KW=templ[0][0];
  SW=templ[0][1];

  if (templ!=NULL)
    {
      for (int j=0;j<nrow;j++)
	delete [] templ[j];
      delete [] templ;
    }
}

void WriteTemplate(char * file, double KW, double SW)
  //********************************************************************************
{
  FILE * f=NULL;
  
  try
    {
      f = fopen(file, "wt");
      fprintf(f,FORMAT_DOUBLE,KW);
      fprintf(f,"%c",SEPARE);
      fprintf(f,FORMAT_DOUBLE,SW);
      fprintf(f,"\n");
      
      if (f!=NULL) fclose(f);
    }
  catch(std::exception &e)
    {
      if (f!=NULL) fclose(f);
      sprintf( ErrorMsg, "problem~in~writing~template~file: %.100s~",file);
      throw std::runtime_error( ErrorMsg );
    }
 
}
//***************************  Common.cpp  **********************

