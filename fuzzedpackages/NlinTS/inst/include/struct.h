/*
 * @author Youssef Hmamouche
 *
 * @brief some useful structures and functions
 *
 * @date 01-04-2016
 *
 */

#ifndef STRUCT_H
#define STRUCT_H


#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "exception.h"


namespace Struct
{

class CVChar : public std::vector<char>
{
public:
    CVChar (unsigned Size) : vector<char> (Size) {}
    CVChar () : vector<char> () {}
};

class CMChar : public std::vector<CVChar>
{
public:
    CMChar (unsigned Size) : vector<CVChar> (Size) {}
    CMChar () : vector<CVChar> () {}
};
class CVDouble : public std::vector<double>
{
public:
    CVDouble            (unsigned Size = 0) : vector<double> (Size) {}
    //CVDouble            () : vector<double> () {}
    double Mean         () const;    //throw (Exception);
    bool   Contains     (unsigned & x);
    double CMean        () const; //    throw (Exception);
    double StdDev       () const; //    throw (Exception);
    double Min          () const; //   throw (Exception);
    double Max          () const; //    throw (Exception);
    void Standardise    (); //          throw (Exception);
    void Normalise  (); //              throw (Exception);
    void Add            (double & m);
    bool NBR_NAN        () const;
};
class CMatDouble : public std::vector<CVDouble>
{
public:
    CMatDouble (unsigned Size /* = 0 */) : vector<CVDouble>(Size) {}
    CMatDouble () : vector<CVDouble>() {}
    void Init_Mat( const vector< vector<double> > & M);
    vector< vector<double> > to_Mat ();
    void Standardise  ();
    CMatDouble Normalise  ();
    void Denormalising (const CMatDouble & minMax);
    void Interpol () ; //             throw (Exception);
};


void permute( CVDouble &X , CVDouble &Y);

CMatDouble Trans (const CMatDouble & M); // Transpose Matrix

bool Trig( CMatDouble & A , CMatDouble & B); // Triangularisation

double Det (const CMatDouble & M); // Determinant

void Resolve(const CMatDouble  & A , const CVDouble  & B, CVDouble & X); // Resolution of the system A * X = B

bool Inverse ( const CMatDouble & B, CMatDouble & M); // throw (Exception); //  Inversed Matrix

double Quartil_1 (const CVDouble &T);

double Quartil_3 (const CVDouble &T);

void boxPlotOutliersDetection (CMatDouble &M,
                               unsigned int fstd);
void algebraicOutlier (CMatDouble & M, unsigned fstd);

} // namespace Struct

#endif // STRUCT_H
