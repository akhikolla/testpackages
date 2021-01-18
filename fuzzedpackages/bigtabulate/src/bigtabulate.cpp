#include <math.h>

#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <sstream>

#include <Rcpp.h>
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/isna.hpp"


#include "bigmemory/util.h"

template<typename T>
string ttos(T i)
{
  stringstream s;
  s << i;
  return s.str();
}

template<>
string ttos<char>(char i)
{
  stringstream s;
  s << static_cast<short>(i);
  return s.str();
}

SEXP StringVec2RChar( const vector<string> &strVec )
{
  if (strVec.empty())
    return R_NilValue;
  SEXP ret = Rf_protect(Rf_allocVector(STRSXP, strVec.size()));
  vector<string>::size_type i;
  for (i=0; i < strVec.size(); ++i)
  {
    SET_STRING_ELT(ret, i, Rf_mkChar(strVec[i].c_str()));
  }
  Rf_unprotect(1);
  return ret;
}

std::vector<std::string> RDouble2StringVec( SEXP numerics)
{
  vector<string> ret( Rf_length(numerics) );
  vector<string>::size_type i;
  for (i=0; i < ret.size(); ++i)
  {
    ret[i] = ttos(REAL(numerics)[i]);
  }
  return ret;
}

std::vector<std::string> RInteger2StringVec( SEXP numerics)
{
  vector<string> ret( Rf_length(numerics) );
  vector<string>::size_type i;
  for (i=0; i < ret.size(); ++i)
  {
    ret[i] = ttos(INTEGER(numerics)[i]);
  }
  return ret;
}


template<typename T>
class Mapper
{
  public:
    typedef std::size_t size_type;

  public:
    virtual int to_index( const T value ) const=0;

    size_type size() const {return _size;}

  protected:
    size_type _size;
};

template<typename T>
class IndexMapper : public Mapper<T>
{
  public:
    // These should be given in sorted order.
    // The NA should appear at the end.
    // useNA indicates that there is an NA
    // and it should be used.
    IndexMapper( T *pFirst, T *pLast, bool useNA )
    {
      _begin = pFirst;
      _end = pLast;
      _useNA = useNA;
      Mapper<T>::_size = distance(_begin, _end);
    } 

    virtual int to_index( const T value ) const
    {
      if (isna(value))
      {
        return _useNA ? std::distance(_begin, _end)+1 : -1;
      }  
      return distance(_begin, std::lower_bound( _begin, 
        _end - static_cast<std::size_t>(_useNA), value));
    }

  protected:
    T* _begin;
    T* _end;
    bool _useNA;
};

template<typename T>
class BreakMapper : public Mapper<T>
{
  public:
    BreakMapper( double min, double max, double numBreaks, bool useNA ) 
      : _min(min), _useNA(useNA)
    {
      _width = max - _min;
      _breakWidth = _width / (numBreaks-1);
      _totalBreaks= numBreaks-1;
      _naIndex = static_cast<index_type>(_totalBreaks+1);
      Mapper<T>::_size = static_cast<std::size_t>(_totalBreaks) +
        static_cast<std::size_t>(_useNA);

    }

    virtual int to_index( const T value ) const
    {
      if (isna(value))
      {
        return _useNA ? _naIndex : -1;
      }
      int bin = static_cast<int>(
        (static_cast<double>(value)-_min) / _breakWidth);
      if (bin < 0 || bin > _totalBreaks)
        return -1;
      return bin;
    }

  protected:
    double _min;
    double _width;
    double _breakWidth;
    double _totalBreaks; // The total number of valid breaks (not including NA).
    bool _useNA;
    index_type _naIndex;
    std::vector<int> _bins;
};

template<typename T>
struct NAMaker;

template<>
struct NAMaker<char>
{char operator()() const {return NA_CHAR;}};

template<>
struct NAMaker<short>
{short operator()() const {return NA_SHORT;}};

template<>
struct NAMaker<int>
{int operator()() const {return NA_INTEGER;}};

template<>
struct NAMaker<double>
{double operator()() const {return NA_REAL;}};

template<typename ValueType, typename InputIter>
std::vector<ValueType> get_unique( const InputIter itStart, 
  const InputIter itEnd, const int includeNA )
{
  InputIter it;
  typedef std::vector<ValueType> Values;
  bool naAdded=false;
  NAMaker<ValueType> make_na;
  Values v;
  if (itStart == itEnd)
    return v;
  for (it = itStart; it != itEnd; ++it)
  {
    if (isna(*it))
    { 
      if (includeNA > 0 && !naAdded)
      {
        v.push_back(make_na());
        naAdded=true;
      }
    }
    else
    {
      typename Values::iterator cit = std::lower_bound( v.begin(), 
        v.end()- static_cast<std::size_t>(naAdded), *it );
      // If we can't find it, we need to add it.
      if (cit == v.end() || *cit != *it) 
      {
        v.insert(cit, *it);
      }
    }
  }
  if (includeNA==2 && !naAdded)
  {
    v.push_back(make_na());
  }
  return v;
}

template<typename RType, typename MatrixAccessorType>
SEXP UniqueGroups( MatrixAccessorType m, SEXP columns, 
  SEXP breakSexp, SEXP useNA )
{
  double *pBreaks = REAL(breakSexp);
  NewVec<RType> RNew;
  VecPtr<RType> RData;
  NAMaker<typename MatrixAccessorType::value_type> make_na;
  typedef std::vector<typename MatrixAccessorType::value_type> Values;
  index_type i,j;
  MatrixAccessor<double> breaks( pBreaks, 3 );
  SEXP ret = Rf_protect(Rf_allocVector(VECSXP, Rf_length(columns)));
  int protectCount = 1;
  Values v;
  const index_type breaksIndex = 2;

  index_type column;
  for (i=0; i < Rf_length(columns); ++i)
  {
    SEXP sv;
    column = static_cast<index_type>(REAL(columns)[i])-1;
    if ( !isna(breaks[i][0]) )
    {
      v.resize(static_cast<std::size_t>(breaks[i][breaksIndex]));
      for (j=0; j < breaks[i][breaksIndex]; ++j)
      {
        v[j] = static_cast<typename MatrixAccessorType::value_type>(j);
      }
      if (INTEGER(useNA)[0] == 1)
      {
        for (j=0; j < m.nrow(); ++j)
        {
          if (isna(m[column][j]))
          {
            v.push_back(make_na());
            break;
          }
        }
      }
      else if (INTEGER(useNA)[0] == 2)
      { 
        v.push_back(make_na());
      }
    }
    else
    {
      v = get_unique<typename MatrixAccessorType::value_type>( 
        (m[column]), (m[column] + m.nrow()), INTEGER(useNA)[0] );
    }
    sv = RNew(v.size());
    std::copy( v.begin(), v.end(), RData(sv) );
    SET_VECTOR_ELT( ret, i, sv );
  }
  Rf_unprotect(protectCount);
  return(ret);
}

template<typename T>
double stable_mean( T* pv, const std::vector<double> &rows, 
  LDOUBLE mean1 )
{
  std::size_t i;
  if ( R_FINITE(static_cast<double>(mean1) ) )
  {
    LDOUBLE t= 0.;
    T v;
    for (i=0; i < rows.size(); ++i)
    {
      v = pv[static_cast<index_type>(rows[i])-1];
      if (!isna(v)) t += static_cast<LDOUBLE>(v) - mean1;
    }
    return mean1 + t / static_cast<LDOUBLE>(rows.size());
    
  }
  else
  {
    return mean1;
  }
}

template<typename T>
double var( T* pv, const std::vector<double> &rows, double mean )
{
  if (rows.size() == 0)
    return NA_REAL;

  LDOUBLE s=0.0;
  std::size_t i, naCount=0;
  T v;
  for (i=0; i < rows.size(); ++i)
  {
    v = pv[static_cast<index_type>(rows[i])-1];
    if (!isna(v)) 
    {
      s+=(static_cast<double>(v)-mean)*(static_cast<double>(v)-mean);
    }
    else
    {
      ++naCount;
    }
  }
  return static_cast<double>(s/(static_cast<LDOUBLE>(rows.size()-naCount)-1.0));
}

template<typename RType, typename MatrixAccessorType, typename MappersType>
std::string MakeIndexLevelName( MatrixAccessorType &m, 
  index_type i, std::vector<bool> &isBreakMapper, 
  MappersType &mappers, SEXP columns )
{
  RType val = static_cast<RType>(
    (m[static_cast<index_type>(REAL(columns)[0]-1)][i]));
  string ret( isBreakMapper[0] ? ttos(mappers[i]->to_index(val)) : ttos(val) );
  int j;
  for (j=1; j < Rf_length(columns); ++j)
  {
    val = static_cast<RType>(
      (m[static_cast<index_type>(REAL(columns)[j]-1)][i]));
    ret += ":" + 
      (isBreakMapper[j] ? ttos(mappers[j]->to_index(val)) : ttos(val));
  }
  return ret;
}

template<typename T>
struct zero_size : public std::unary_function<T, bool>
{
  bool operator()( const T &vec ) const {return vec.size() == 0;}
};

typedef std::vector<std::string> strings;

strings interact( const strings &s1, SEXP s2 )
{
  std::vector<std::string> ret( s1.size() * Rf_length(s2), string("") );
  size_t i, j, k;
  k=0;
  for (i=0; i < s1.size(); ++i)
  {
    if (Rf_isInteger(s2))
    {
      for (j=0; j < static_cast<size_t>(Rf_length(s2)); ++j)
      {
        ret[k++] = ttos(INTEGER(s2)[j]) + ":" + s1[i];
      }
    }
    else
    {
      for (j=0; j < static_cast<size_t>(Rf_length(s2)); ++j)
      {
        ret[k++] = ttos(REAL(s2)[j]) + ":" + s1[i];
      }
    }
  }
  return ret;
}

template<typename RType, typename MatrixAccessorType>
SEXP TAPPLY( MatrixAccessorType m, SEXP columns, SEXP breakSexp,
  SEXP returnTable, SEXP useNA, 
  SEXP returnSummary, SEXP processColumns, SEXP summaryNARM,
  SEXP splitcol, SEXP splitlist )
{
  std::vector<std::string> retNames;
  retNames.push_back(std::string("levels"));
  SEXP uniqueGroups=Rf_protect(UniqueGroups<RType>(m, columns, breakSexp, useNA));
  int i, j, k;
  strings groupNames;
  if (Rf_isInteger(VECTOR_ELT(uniqueGroups, Rf_length(uniqueGroups)-1)))
  {
    groupNames = 
      RInteger2StringVec(VECTOR_ELT(uniqueGroups, Rf_length(uniqueGroups)-1));
  }
  else
  {
    groupNames = 
      RDouble2StringVec(VECTOR_ELT(uniqueGroups, Rf_length(uniqueGroups)-1));
  }
  for (i=Rf_length(uniqueGroups)-2; i >= 0; --i)
  {
    groupNames = interact( groupNames, VECTOR_ELT(uniqueGroups, i) );
  }
  std::map<std::string, int> lmi;
  i=0;
  lmi["levels"] = i++;
  if ( splitcol != R_NilValue )
  {
    lmi["split"] = i++;
    retNames.push_back(std::string("split"));
  }
  if ( LOGICAL(returnTable)[0] )
  {
    lmi["table"] = i++;
    retNames.push_back(std::string("table"));
  }
  if ( LOGICAL(returnSummary)[0] )
  {
    lmi["summary"] = i++;
    retNames.push_back(std::string("summary"));
  }
  SEXP ret = Rf_protect(Rf_allocVector(VECSXP, i));
  Rf_setAttrib( ret, R_NamesSymbol, StringVec2RChar( retNames ) );
  SET_VECTOR_ELT( ret, lmi[string("levels")], uniqueGroups );
  MatrixAccessor<double> breaks( REAL(breakSexp), 3 );
  typedef boost::shared_ptr<Mapper<RType> > MapperPtr;
  typedef std::vector<MapperPtr> Mappers;
  VecPtr<RType> RData;
  NewVec<RType> RNew;
  Mappers mappers;
  int protectCount=2;
  int totalListSize =0;
  std::vector<int> accMult;
  // Create the data structures that map values to indices for each of the
  // columns.
  std::vector<bool> isBreakMapper(Rf_length(uniqueGroups), false);
  accMult.resize(Rf_length(uniqueGroups));
  int lastVecLen=0;
  for (i=0; i < Rf_length(uniqueGroups); ++i)
  {
    SEXP vec = VECTOR_ELT(uniqueGroups, i);
    int vecLen = Rf_length(vec);
    if (!isna(breaks[i][0]))
    {
      mappers.push_back( MapperPtr(
        new BreakMapper<RType>( breaks[i][0], breaks[i][1], 
          breaks[i][2], INTEGER(useNA)[0] > 0 ) ) );
      isBreakMapper[i]=true;
    }
    else
    {
      mappers.push_back( MapperPtr( 
        new IndexMapper<RType>( RData(vec), 
          RData(vec) + vecLen, INTEGER(useNA)[0] > 0 ) ) );
    }
    totalListSize = (totalListSize == 0 ? vecLen : totalListSize*vecLen);
    if (i==0)
    {
      accMult[i] = 1;
      lastVecLen = vecLen;
    }
    else
    {
      accMult[i] = accMult[i-1] * lastVecLen;
      lastVecLen = vecLen;
    }
  }
  typedef std::vector<double> Indices;
  typedef std::vector<Indices> TableIndices;
  TableIndices tis;

  typedef std::vector<RType> Values;
  typedef std::vector<Values> TableIndexValues;
  TableIndexValues tiv;

  Indices tvs(totalListSize, 0);
  
  typedef std::vector<double> TableSummary;
  typedef std::vector<TableSummary> TableSummaries;
  std::vector<TableSummaries> ts;
  
  typedef std::vector<index_type> ProcessColumns;
  ProcessColumns procCols;

  if ( splitcol != R_NilValue || LOGICAL(returnSummary)[0] )
  {
    if ( isna(REAL(splitcol)[0]) || LOGICAL(returnSummary)[0] )
    {
      tis.resize(totalListSize);
    }
    else
      tiv.resize(totalListSize);
  }
  if ( LOGICAL(returnTable)[0] || LOGICAL(returnSummary)[0] )
  {
    tvs.resize(totalListSize);
    std::fill( tvs.begin(), tvs.end(), 0 );
  }
  if ( LOGICAL(returnSummary)[0] )
  {
    procCols.resize(Rf_length(processColumns));
    for (k=0; k < static_cast<index_type>(procCols.size()); ++k)
    {
      procCols[k] = static_cast<index_type>(REAL(processColumns)[k])-1;
    }
    ts.resize( procCols.size() );
    // min, max, sum, sum^2
    std::fill( ts.begin(), ts.end(), 
      TableSummaries(totalListSize, TableSummary(6, 0.)) );
  }
  // Get the indices for each of the column-value combinations.
  for (i=0; i < m.nrow(); ++i)
  {
    int tableIndex=0;
    int mapperVal;
    for (j=1; j < Rf_length(columns); ++j)
    {
      mapperVal = mappers[j]->to_index( static_cast<RType>(
          (m[static_cast<index_type>(REAL(columns)[j]-1)][i])) );
      if (mapperVal == -1)
      {
        tableIndex = -1;
        break;
      }
      tableIndex += accMult[j] * mapperVal;
    }
    mapperVal = mappers[0]->to_index( static_cast<RType>(
        (m[static_cast<index_type>(REAL(columns)[0]-1)][i])) );
    if (tableIndex == -1 || mapperVal == -1)
    {
      continue;
    }
    tableIndex += mapperVal;
    if ( splitcol != R_NilValue || LOGICAL(returnSummary)[0] )
    {
      if ( isna(REAL(splitcol)[0]) || LOGICAL(returnSummary)[0] )
      {
        tis[tableIndex].push_back(i+1);
      }
      else
      {
        tiv[tableIndex].push_back( 
          m[static_cast<index_type>(REAL(splitcol)[0])-1][i] );
      }
    }
    if ( LOGICAL(returnTable)[0] || LOGICAL(returnSummary)[0] )
    {
      ++tvs[tableIndex];
    }
    if ( LOGICAL(returnSummary)[0] )
    {
      for (k=0; k < static_cast<index_type>(ts.size()); ++k)
      {
        // index 5 is the na count
        // index 4 means first min value is found
        // index 3 means first max value is found
        // index 2 is the running sum
        TableSummaries &ss = ts[k];
        double matVal = static_cast<double>(
          m[static_cast<index_type>(procCols[k])][i]);
        if (isna(static_cast<typename MatrixAccessorType::value_type>(matVal)))
        {
          ++ss[tableIndex][5];
          continue;
        }
        if (ss[tableIndex][3] == 0) 
        {
          ss[tableIndex][0] = matVal;
          ss[tableIndex][3] = 1;
        }
        else
        {
          if (ss[tableIndex][0] > matVal)
            ss[tableIndex][0] = matVal;
        }
        if (ss[tableIndex][4] == 0) 
        {
          ss[tableIndex][1] = matVal;
          ss[tableIndex][4] = 1;
        }
        else
        {
          if (ss[tableIndex][1] < matVal)
            ss[tableIndex][1] = matVal;
        }
        ss[tableIndex][2] += matVal;
      }
    }
  }
  if ( splitcol != R_NilValue )
  { 
    SEXP mapRet;
    if ( INTEGER(splitlist)[0] == 1)
    {
      SEXP vec;
      // Copy to a list of vectors that R can read.
      mapRet = Rf_protect(Rf_allocVector(VECSXP, groupNames.size() ));
      ++protectCount;
      if ( isna(REAL(splitcol)[0]) )
      {
        for (i=0; i < static_cast<index_type>(tis.size()); ++i)
        {
          Indices &ind = tis[i];
          vec = Rf_allocVector(REALSXP,ind.size());
          std::copy( ind.begin(), ind.end(), REAL(vec) );
          SET_VECTOR_ELT( mapRet, i, vec );
        }
      }
      else
      {
        for (i=0; i < static_cast<index_type>(tiv.size()); ++i)
        {
          Values &ind = tiv[i];
          vec = RNew(tiv[i].size());
          std::copy( ind.begin(), ind.end(), RData(vec) );
          SET_VECTOR_ELT( mapRet, i, vec );
        }
      }
      Rf_setAttrib(mapRet, R_NamesSymbol, StringVec2RChar(groupNames));
      groupNames.clear();
      groupNames.reserve(0);
    }
    else if ( INTEGER(splitlist)[0] == 2)
    {
      std::size_t zeroCount;
      if (tis.size() > 0)
      {
        zeroCount=std::count_if( tis.begin(), tis.end(),
          zero_size<typename TableIndices::value_type>() );
      }
      else
      {
        zeroCount=std::count_if( tiv.begin(), tiv.end(),
          zero_size<typename TableIndexValues::value_type>() );
      }

      SEXP vec;
      // Copy to a list of vectors that R can read.
      mapRet = Rf_protect(Rf_allocVector(VECSXP,groupNames.size()-zeroCount ));
      ++protectCount;
      SEXP mapNames = Rf_protect(Rf_allocVector(STRSXP, Rf_length(mapRet)));
      ++protectCount;
      j=0;
      if ( isna(REAL(splitcol)[0]) )
      {
        for (i=0; i < static_cast<index_type>(tis.size()); ++i)
        {
          if (tis[i].size() > 0)
          {
            Indices &ind = tis[i];
            vec = Rf_allocVector(REALSXP,tis[i].size());
            std::copy( ind.begin(), ind.end(), REAL(vec) );
            SET_VECTOR_ELT( mapRet, j, vec );
            SET_STRING_ELT( mapNames, j, Rf_mkChar(groupNames[i].c_str()) );
            ++j;
          }
        }
      }
      else
      {
        for (i=0; i < static_cast<index_type>(tiv.size()); ++i)
        {
          if (tiv[i].size() > 0)
          {
            Values &ind = tiv[i];
            vec = RNew(tiv[i].size());
            std::copy( ind.begin(), ind.end(), RData(vec) );
            SET_VECTOR_ELT( mapRet, j, vec );
            SET_STRING_ELT( mapNames, j, Rf_mkChar(groupNames[i].c_str()) );
            ++j;
          }
        }
      }
      Rf_setAttrib(mapRet, R_NamesSymbol, mapNames);
      groupNames.clear();
      groupNames.reserve(0);
    }
    else // INTEGER_VALUE(splitlist) == 0 
    {
      SEXP mn = Rf_protect(Rf_allocVector(STRSXP, m.nrow()));
      ++protectCount;
      mapRet = Rf_allocVector(REALSXP,m.nrow());
      double *pmr = REAL(mapRet);
      for (i=0; i < static_cast<index_type>(tis.size()); ++i)
      {
        Indices &inds = tis[i];
        for (j=0; j < static_cast<index_type>(inds.size()); ++j)
        {
          SET_STRING_ELT(mn, static_cast<int>(inds[j]-1), 
            Rf_mkChar(groupNames[i].c_str()));
          pmr[static_cast<index_type>(inds[j])-1] = i+1;
        }
      }
      groupNames.clear();
      groupNames.reserve(0);
      Rf_setAttrib(mapRet, R_NamesSymbol, mn);
    }
    SET_VECTOR_ELT(ret, lmi[string("split")], mapRet);
  }
  if ( LOGICAL(returnTable)[0] )
  {
    SEXP tableRet = Rf_allocVector(INTSXP,tvs.size());
    std::copy( tvs.begin(), tvs.end(), INTEGER(tableRet) );
    SET_VECTOR_ELT(ret, lmi[string("table")], tableRet);
  } 
  if ( LOGICAL(returnSummary)[0] )
  {
    // If we change the data structures holding the summaries, can
    // we get better performance.
    std::vector<std::string> colnames;
    colnames.push_back(string("min"));
    colnames.push_back(string("max"));
    colnames.push_back(string("mean"));
    colnames.push_back(string("sd"));
    colnames.push_back(string("NAs"));
    SEXP summaryRet = Rf_protect(Rf_allocVector(VECSXP, ts[0].size()));
    ++protectCount;
    for (i=0; i < Rf_length(summaryRet); ++i)
    {
      SEXP dimnames = Rf_allocVector(VECSXP, 2);
      SET_VECTOR_ELT(dimnames, 0, R_NilValue );
      SET_VECTOR_ELT(dimnames, 1, StringVec2RChar(colnames) );
      
      SEXP retMat = Rf_allocMatrix(REALSXP, ts.size(), 5);
      Rf_setAttrib(retMat, R_DimNamesSymbol, dimnames);
      MatrixAccessor<double> rm( REAL(retMat), ts.size() );
      LDOUBLE temp;
      for (j=0; j < static_cast<index_type>(ts.size()); ++j)
      {
        if (tvs[i] > 0)
        {
          if (LOGICAL(summaryNARM)[0] || ts[j][i][5] == 0)
          {
            rm[0][j] = ts[j][i][0];
            rm[1][j] = ts[j][i][1];
            temp = stable_mean( m[procCols[j]], tis[i], 
              static_cast<LDOUBLE>(ts[j][i][2]) / 
                static_cast<LDOUBLE>(tvs[i]-ts[j][i][5]));
            rm[2][j] = static_cast<double>(temp);
            rm[3][j] = sqrt( var( m[procCols[j]], tis[i], temp ) );
            rm[4][j] = ts[j][i][5];
          }
          else
          {
            rm[0][j] = NA_REAL;
            rm[1][j] = NA_REAL;
            rm[2][j] = NA_REAL;
            rm[3][j] = NA_REAL;
            rm[4][j] = ts[j][i][5];
          }
        }
        else
        {
          rm[0][j] = NA_REAL;
          rm[1][j] = NA_REAL;
          rm[2][j] = NA_REAL;
          rm[3][j] = NA_REAL;
          rm[4][j] = 0;
        }
      }
      SET_VECTOR_ELT(summaryRet, i, retMat);
    }
    SET_VECTOR_ELT(ret, lmi[string("summary")], summaryRet);
  }
  Rf_unprotect( protectCount );
  return ret;
}

extern "C"
{

SEXP RNumericTAPPLY( SEXP numericMatrix , SEXP columns, SEXP breaks,
  SEXP returnTable, SEXP useNA, SEXP returnSummary, 
  SEXP processColumns, SEXP summaryNARM, SEXP splitcol, SEXP splitlist )
{
  return TAPPLY<double>( MatrixAccessor<double>( REAL(numericMatrix),
    static_cast<index_type>(Rf_nrows(numericMatrix)) ), 
    columns, breaks, returnTable, useNA, 
    returnSummary, processColumns, summaryNARM, splitcol, splitlist );
}

SEXP RIntTAPPLY( SEXP numericMatrix , SEXP columns, SEXP breaks,
  SEXP returnTable, SEXP useNA, SEXP returnSummary, 
  SEXP processColumns, SEXP summaryNARM, SEXP splitcol, SEXP splitlist )
{
  return TAPPLY<int>( MatrixAccessor<int>( INTEGER(numericMatrix),
    static_cast<index_type>(Rf_nrows(numericMatrix)) ), 
    columns, breaks, returnTable, useNA, 
    returnSummary, processColumns, summaryNARM, splitcol, splitlist );
}

// Return both the table indices in a list and the unique levels
//useNa =0 "no", 1 "ifany", 2 "always".
//ccols breaks, boolean return.map, boolean table, integer useNA,
//boolean summary?, numeric summary columns, boolean summary.na.rm
SEXP BigMatrixTAPPLY( SEXP bigMatAddr, SEXP columns, SEXP breaks, 
  SEXP returnTable, SEXP useNA, SEXP returnSummary, 
  SEXP processColumns, SEXP summaryNARM, SEXP splitcol, SEXP splitlist )
{
  BigMatrix *pMat = reinterpret_cast<BigMatrix*>(
    R_ExternalPtrAddr(bigMatAddr));
  if (pMat->separated_columns())
  {
    switch (pMat->matrix_type())
    {
      case 1:
        return TAPPLY<int>( SepMatrixAccessor<char>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 2:
        return TAPPLY<int>( SepMatrixAccessor<short>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 4:
        return TAPPLY<int>( SepMatrixAccessor<int>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 8:
        return TAPPLY<double>( SepMatrixAccessor<double>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
    }
  }
  else
  {
    switch (pMat->matrix_type())
    {
      case 1:
        return TAPPLY<int>( MatrixAccessor<char>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 2:
        return TAPPLY<int>( MatrixAccessor<short>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 4:
        return TAPPLY<int>( MatrixAccessor<int>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
      case 8:
        return TAPPLY<double>( MatrixAccessor<double>(*pMat), 
          columns, breaks, returnTable, useNA,
          returnSummary, processColumns, summaryNARM, splitcol, splitlist );
    }
  }
  return R_NilValue;
}

}
