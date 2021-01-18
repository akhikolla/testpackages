//*--------------------------------------------------------------------------*//
//  File:               routines.h
//  Created by:         Pavlo Mozharovskyi
//  First released:     01.11.2018
//
//  Header to the file routines.cpp that contains routines for computing Tukey
//  curve depth and related functions.
//
//  Subsequent changes are listed below:
//  01.11.2018 (Pavlo Mozharovskyi): First version released.
//*--------------------------------------------------------------------------*//


#include <iostream>
#include <memory.h>
#include "hypermatrix.h"

// How much more can a function contain voxels than grid discretes
const int timesGrid = 10;
const double eps = 1e-8;

// Sorting record
struct RecSort{
  int index;
  double weight;
  double value;
};

// Sorting structures and routines
struct RecMatrix{
  double val;
  int iRow;
  int iCol;
};

int Compare(RecSort &rs1, RecSort &rs2);

void Swap(RecSort* rs1, RecSort* rs2);

int Compare(RecMatrix &rm1, RecMatrix &rm2);

void Swap(RecMatrix* rm1, RecMatrix* rm2);

/* -------------------------------------------------------------------------- */
/* quickSort from http://www.proggen.org/doku.php?id=algo:quicksort           */
/* (modified, templated)                                                      */
/* 'left' and 'right' are the elements' numbers to eb accessed                */
/* -------------------------------------------------------------------------- */
template<typename T>
void quick_sort(T *values, int left, int right, int(*cmp)(T& x, T& y),
                void(*swap)(T* x, T* y)) {
  int i = left, j = right; // Z?hlindizes
  // Pivot-Element (Array-Mitte) bestimmen
  T pivot = values[(left + right) >> 1];
  // Solange Paare suchen und tauschen, bis sich die Z?hlindizes "getroffen"
  // haben
  do {
    // Paar finden, dass getauscht werden muss
    while (cmp(values[i], pivot)) { ++i; }
    while (cmp(pivot, values[j])) { --j; }
    // Wenn sich die Z?hlindizes noch nicht "getroffen" haben, dann
    // tauschen und weiterz?hlen. Sollten die Z?hlindizes gleich sein, dann
    // z?hle einfach weiter und unterbrich die Schleife
    if (i < j) {
      swap(&values[i], &values[j]); ++i; --j;
    }
    else { if (i == j) { ++i; --j; break; } }
  } while (i <= j);
  // Wenn die Teillisten mehr als ein Element enthalten, dann wende quickSort
  // auf sie an
  if (left < j) { quick_sort(values, left, j, cmp, swap); }
  if (i < right) { quick_sort(values, i, right, cmp, swap); }
}
template<typename T>
void quick_sort(T *values, int left, int right) {
  quick_sort(values, left, right, Compare, Swap);
}

// Type of curve
enum CurveType{
  EMPTYCURVE = 0,
  LINSEGMENTS = 1, // Curve defined by connected line segments
  VOXELS = 2, // Curve defined by voxel numbers
  VOXCOORDS = 3, // Curve defined by voxel coordinates
  IMAGEDENSITY = 4 // Density in the space
};

// Curve structure where fields are used depending on the curve type
struct Curve{
  CurveType type;
  int n;
  int d;
  double* args;
  double** vals;
  int** voxels;
  double* rawVals;
  int* rawVoxels;
  bool closed;
  Curve(){
    // std::cout << "Constructor Curve() called." << std::endl;
    this->type = EMPTYCURVE;
    this->d = 0;
    this->n = 0;
    this->args = 0;
    this->vals = 0;
    this->voxels = 0;
    this->rawVals = 0;
    this->rawVoxels = 0;
    this->closed = false;
  }
  Curve(int d){
    // std::cout << "Constructor Curve(int d) called." << std::endl;
    this->type = EMPTYCURVE;
    this->d = d;
    this->n = 0;
    this->args = 0;
    this->vals = 0;
    this->voxels = 0;
    this->rawVals = 0;
    this->rawVoxels = 0;
    this->closed = false;
  }
  Curve(int d, int n, double* args, double* rawVals){
    this->type = LINSEGMENTS;
    this->d = d;
    this->n = n;
    this->closed = false;
    this->args = new double[this->n];
    memcpy(this->args, args, n * sizeof(double));
    this->rawVals = new double[this->n * this->d];
    memcpy(this->rawVals, rawVals, n * d * sizeof(double));
    this->vals = new double*[this->d];
    for (int i = 0; i < this->d; i++){
      this->vals[i] = this->rawVals + i * this->n;
    }
  }
  Curve(int d, int n, int* rawVoxels){
    this->type = VOXELS;
    this->d = d;
    this->n = n;
    this->closed = false;
    this->rawVoxels = new int[this->n * this->d];
    memcpy(this->rawVoxels, rawVoxels, n * d * sizeof(int));
    this->voxels = new int*[this->n];
    for (int i = 0; i < this->n; i++){
      this->voxels[i] = this->rawVoxels + i * this->d;
    }
  }
  Curve(int d, int n, double* rawVals){
    this->type = VOXCOORDS;
    this->d = d;
    this->n = n;
    this->closed = false;
    this->rawVals = new double[this->n * this->d];
    memcpy(this->rawVals, rawVals, n * d * sizeof(double));
    this->vals = new double*[this->n];
    for (int i = 0; i < this->n; i++){
      this->vals[i] = this->rawVals + i * this->d;
    }
  }
  double distHausdorff(const Curve &curve);
  double distCurve(const Curve &curve, bool oneWay);
  ~Curve(){
    //Rcout << "Curve ";
    if (this->vals){
      //Rcout << "destruction";
      delete[] this->rawVals;
      delete[] this->vals;
    }
    //Rcout << std::endl;
  }
};

// A function that provides doule-indexing.
template<typename T> T* asMatrix(T arr, int n, int d){
  T* mat = new T[n];
  for (int i = 0; i < n; i++)
    mat[i] = arr + i*d;
  return mat;
}

// Tranpose for IntegerMatrix / NumericMatrix, see array.c in R
template <typename T> T transpose(const T & m){
  int k = m.rows();
  int n = m.cols();
  T z(n, k);
  int sz1 = n * k - 1;
  typename T::const_iterator mit ;
  typename T::iterator zit;
  for (mit = m.begin(), zit = z.begin(); mit != m.end(); mit++, zit += n){
    if (zit >= z.end()){
      zit -= sz1;
    }
    *zit = *mit;
  }
  return(z);
}

class ImageDensity : public typeHypermatrix<double>, public Curve{
public:
  int d;
  ImageDensity(int d, int* ns) : typeHypermatrix<double>(d, ns), Curve(d) {
    // std::cout << "Constructor ImageDensity(int d, int* ns) called." << std::endl;
    this->type = IMAGEDENSITY;
    this->d = d;
    this->n = this->size;
    // std::cout << this->n << " " << this->rawVals << " " << this->vals << std::endl;
    this->rawVals = new double[this->n * this->d];
    this->vals = new double*[this->n];
    for (int i = 0; i < this->n; i++){
      this->vals[i] = this->rawVals + i * this->d;
    }
    // The array of voxels' numbers
    int* counters = new int[this->d];
    for(int i = 0; i < this->d - 1; i++){
      counters[i] = 0;
    }
    counters[this->d - 1] = -1;
    // Go through all cells
    for (int i = 0; i < this->n; i++){
      // Choose current cell
      counters[this->d - 1]++;
      for (int j = this->d - 1; j > 0; j--){
        if (counters[j] == this->ns[j]){
          counters[j] = 0;
          counters[j - 1]++;
        }else{
          break;
        }
      }
      // for (int j = 0; j < d; j++){
      //   std::cout << counters[j] << " ";
      // }
      // std::cout << std::endl;
      // For the current cell do:
      for (int j = 0; j < this->d; j++){
        // Define voxels' centers
        vals[i][j] = ((double)counters[j] + 0.5) / (double)this->ns[j];
        // std::cout << vals[i][j] << " ";
      }
    }
    // std::cout << std::endl;
    delete[] counters;
  }
  ImageDensity(int d, int* ns, double* voxelDns) : ImageDensity(d, ns){
    // std::cout << "Constructor ImageDensity(int d, int* ns, double* voxelDns) called." << std::endl;
    memcpy(this->body, voxelDns, this->n * sizeof(double));
  }
  Curve* toCurve();
  Curve* toOrderedCurve();
};

class EmpDist : public Curve{
public:
  double* weights;
  EmpDist(int n, int d, bool zeroInit);
  EmpDist(EmpDist &empDist);
  EmpDist(Curve* curves, int nCurves, double precision);
  EmpDist(ImageDensity* imageDensities, int nImageDensities, double precision);
  void updateWeights(bool dropZeros);
  ~EmpDist(){
    //Rcout << "EmpDist ";
    if (weights){
      //Rcout << "destruction";
      delete[] weights;
    }
    //Rcout << std::endl;
  }
private:
  EmpDist();
};

double calcOneDepth(Curve curve, double** curvePrj, int nDirs,
                    Curve* curves, double*** curvePrjs, int n);

double calcOneDepth(ImageDensity &object, double** objectPrj, int nDirs,
                    ImageDensity &data, double** dataPrj, bool sbSmpl);

double curvePortion(double* curvePrj, double pointPrj, int nVoxels);

double imagePortion(double* imagePrj, double* imageDns, double pointPrj,
                    int nVoxels);

int generateDirections(int seed, int k, int d, double** dirs);

int projectCurveVoxels(const Curve &curve, int nDirs, double** dirs,
                       double** prjs);

double getMinMaxDist(int n1, int n2, double* dists);

double calcOneDepth(EmpDist &curEmpDist, EmpDist &genEmpDist, double** dirs,
                    int nDirs, int d);

double calcOneDepth(EmpDist & refEmpDist, EmpDist &curEmpDist,
                    EmpDist &genEmpDist, double** dirs,
                    int nDirs, int d);

double empDistPortion(double* empDistPrj, double* empDistWeights,
                      double pointPrj, int nVoxels);

double approxOneDepth(EmpDist &refEmpDist, EmpDist &curEmpDist,
                      EmpDist &genEmpDist, double** dirs,
                      int nDirs, int d, double curMinMass, double genMinMass);

double calcOneDepth(EmpDist &refEmpDist, EmpDist &curEmpDist,
                    EmpDist &genEmpDist, int d,
                    double minMassObj, double minMassDat);

double calcExPointDepth(double* point, EmpDist &curEmpDist,
                        EmpDist &genEmpDist,
                        double curMinMass, double genMinMass);

double calcExPointDepthRec(double* point, EmpDist &curEmpDist,
                           EmpDist &genEmpDist,
                           double curMinMass, double genMinMass);

double calcExPointDepth2D(double* point, EmpDist &curEmpDist,
                          EmpDist &genEmpDist,
                          double curMinMass, double genMinMass,
                          double curPosMass, double genPosMass,
                          double curNegMass, double genNegMass);

int updateDepth(double* nhMassCur, double nzMassCur, double curPosMass,
                double curNegMass, double* nhMassGen, double nzMassGen,
                double genPosMass, double genNegMass,
                double curMinMass, double genMinMass, double* depth);

double norm2(double* x, int d);
