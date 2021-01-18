/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  DManager
 * Purpose:  Declaration of the class ImportExportToCsv.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 */

/** @file STK_ExportToCsv.h
 *  @brief In this file we define the class ExportToCsv.
 **/

#ifndef STK_EXPORTCSV_H
#define STK_EXPORTCSV_H

#include "STK_DataFrame.h"
#include "STK_ReadWriteCsv.h"

namespace STK
{

namespace Csv
{
  /** @ingroup DManager
   *  Defines the default prefix to used when naming the columns
   *  of the ReadWriteCsv.
   **/
  static const String DEFAULT_COLUMN_PREFIX  =  _T("Var");
}


/** @ingroup DManager
 *  @brief Export data to a Csv.
 *
 * An ExportToCsv object creates a @c ReadWriteCsv from a container of data
 * like a DataFrame, a Vector, a point or a ArrayXX. The data are stored in a
 * String format in the @c ReadWriteCsv struct.
 **/
class ExportToCsv
{
  public:
    /** Default constructor. Create an instance of ExportToCvs. */
    inline ExportToCsv(): p_data_(new ReadWriteCsv()) {}
    /** Constructor : create an instance of ExportToCvs with a DataFrame.
     *  @param df the DataFrame to export
     **/
    inline ExportToCsv( DataFrame const& df): p_data_(new ReadWriteCsv())
    {
      p_data_->setReserve(df.sizeRows());
      p_data_->resize(df.sizeRows(), df.sizeCols());
      // for each field Try a String conversion
      for(int iVar = df.beginCols(), irw = p_data_->begin(); iVar<df.endCols(); iVar++, irw++)
      { if (df.elt(iVar)) df.elt(iVar)->exportAsString(p_data_->var(irw));}
    }
    /** Constructor : create an instance of ExportToCvs with a ReadWriteCsv.
     *  @param rw the TReadWriteCsv to export
     **/
    template<class Type>
    ExportToCsv( TReadWriteCsv<Type> const& rw): p_data_(new ReadWriteCsv())
    {
      p_data_->setReserve(rw.sizeRows());
      p_data_->resize(rw.sizeRows(), rw.sizeCols());
      p_data_->setWithNames(rw.withNames());
      // for each field Try a String conversion
      for(int iRw = rw.beginCols(), iData = p_data_->begin(); iRw<rw.endCols(); iRw++, iData++)
      { rw.var(iRw).exportAsString(p_data_->var(iData));}
    }

    /** Instantiates an instance of ExportToCvs with an Array1D, a list1D, etc....
     *  @param A the 1D container to export
     *  @param byCol export the container as a column vector or a row vector ?
     *  @param prefix the prefix of the name to set to the variable
     **/
    template < class Container>
    ExportToCsv( ITContainer1D<Container> const& A
               , bool byCol = true
               , String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
               : p_data_(new ReadWriteCsv())
    {
#ifdef STK_DMANAGER_DEBUG
      stk_cout << "Entering ExportToCsv( ITContainer1D<Container> const& A, byCol, prefix)\n";
      stk_cout << "A.range()= " << A.range()  << "\n";
#endif
      // add an empty string variable (an empty column)
      p_data_->setWithNames(true);
      if (byCol)
      { // add to ReadWriteCsv a new variable
        //p_data_->
        p_data_->push_back(Variable<String>(A.range(), prefix));

        for(int i = A.begin(); i<A.end(); i++)
        { p_data_->back()[i] = typeToString(A.at(i));}
      }
      else // by row
      {
        for(int i = A.begin(); i<A.end(); i++)
        { p_data_->push_back( Variable<String>(1, typeToString(A.at(i)), prefix+typeToString(i)) );}
      }
    }

    /** Instantiates an instance of ExportToCvs with a vector.
     *  @param A the ITContainer to export
     *  @param byCol export the container as a column vector or a raw vector ?
     *  @param prefix the prefix of the name to set to the variable
     **/
    template < class Container>
    ExportToCsv( ITContainer<Container, Arrays::vector_> const& A
               , bool byCol = true
               , String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
               : p_data_(new ReadWriteCsv())
    {
        // add an empty string variable (an empty column)
        p_data_->setWithNames(true);
        if (byCol)
        {
          p_data_->push_back(Variable<String>(A.range(), prefix));

          for(int i = A.begin(); i<A.end(); i++)
          { p_data_->back()[i] = typeToString(A.at(i));}
        }
        else
        {
          for(int i = A.begin(); i<A.end(); i++)
          {
            p_data_->push_back(Variable<String>(1, prefix+typeToString(i)));
            p_data_->back().front() = typeToString(A.at(i));
          }
        }
    }
    /** Instantiates an instance of ExportToCvs with a point.
     *  @param A the ITContainer1D to export
     *  @param byCol export the container as a column vector or a raw vector ?
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class Container>
    ExportToCsv( ITContainer<Container, Arrays::point_> const& A
               , bool byCol = false
               , String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
               : p_data_(new ReadWriteCsv())
    {
        // add an empty string variable (an empty column)
        p_data_->setWithNames(true);
        if (byCol)
        {
          p_data_->push_back(Variable<String>(A.range(), prefix));

          for(int i = A.begin(); i<A.end(); i++)
          { p_data_->back()[i] = typeToString(A.at(i));}
        }
        else
        {
          for(int i = A.begin(); i<A.end(); i++)
          {
            p_data_->push_back(Variable<String>(1, prefix+typeToString(i)));
            p_data_->back().front() = typeToString(A.at(i));
          }
        }
    }
    /** Instantiates an instance of ExportToCvs with a general array
     *  @param A the IArray2d to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template <class Container >
    ExportToCsv( ITContainer<Container> const& A
               , String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
               : p_data_(new ReadWriteCsv())
    {
      p_data_->setWithNames(true);
      for(int iVar = A.beginCols(); iVar<A.endCols(); iVar++)
      {
        // add an empty string variable (an empty column)
        p_data_->push_back(Variable<String>(A.rows(), prefix));
        for (int iRow=A.beginRows(); iRow<A.endRows(); iRow++)
        { p_data_->back()[iRow] = typeToString(A.at(iRow,iVar));}
      }
    }

    /** destructor.
     *  The protected field p_data_ will be liberated.
     **/
    inline ~ExportToCsv(){ if (p_data_) delete p_data_;}
    /** Append a vector.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class Container>
    void append( ITContainer<Container, Arrays::vector_> const& A, String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
    {
      // add an empty string variable
      p_data_->push_back(Variable<String>(A.range(), prefix));
      // add strings to the String variable
      for(int i = A.begin(); i< A.end(); i++)
      { p_data_->back()[i] = typeToString(A.elt(i));}
    }
    /** Append a point (as a column vector).
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class Container>
    void append( ITContainer<Container, Arrays::point_> const& A, String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
    {
      // add an empty string variable
      p_data_->push_back(Variable<String>(A.range(), prefix));
      // add strings to the String variable
      for(int i = A.begin(); i< A.end(); i++)
      { p_data_->back()[i] = typeToString(A.elt(i));}
    }
    /** Append a 2D container.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class Container>
    void append( ITContainer<Container, Arrays::array2D_> const& A, String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
    {
      // for each field try a String conversion
      for(int iVar = A.beginCols(), iNum=1; iVar<A.endCols(); iVar++, iNum++)
      {
        // add an empty string variable (an empty column)
        p_data_->push_back(Variable<String>(A.rows(), prefix+typeToString(iNum)));
        for (int iRow=A.beginRows(); iRow<A.endRows(); iRow++)
        { p_data_->back()[iRow] = typeToString(A.elt(iRow,iVar));}
      }
    }
    /** Append a 2D container.
     *  @param A the container to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class Container>
    void append( ITContainer<Container, Arrays::square_> const& A, String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
    {
      // for each field try a String conversion
      for(int iVar = A.beginCols(), iNum=1; iVar<A.endCols(); iVar++, iNum++)
      {
        // add an empty string variable (an empty column)
        p_data_->push_back(Variable<String>(A.rows(), prefix+typeToString(iNum)));
        for (int iRow=A.beginRows(); iRow<A.endRows(); iRow++)
        { p_data_->back()[iRow] = typeToString(A.elt(iRow,iVar));}
      }
    }
    /** Append a value.
     *  @param A the value to export
     *  @param prefix the prefix ot the name to set to the variable
     **/
    template < class TYPE>
    void appendData( TYPE const& A, String const& prefix=Csv::DEFAULT_COLUMN_PREFIX)
    {
      // add an empty string variable
      p_data_->push_back(Variable<String>( prefix + typeToString(p_data_->lastIdx())) );
      // add strings to the String variable
      p_data_->back().push_back(typeToString(A));
    }
    /** Set a name to each column of the ReadWriteCsv using the form
     *  prefix + number.
     *  @param prefix the prefix to use for the names of the columns
     **/
    inline void setColumnsNames(String const& prefix = Csv::DEFAULT_COLUMN_PREFIX)
    {
      for(int i = p_data_->begin(); i<p_data_->end(); i++)
      { p_data_->setName(i, prefix + typeToString(i)) ;}
    }

    /** Accesor. Return a ptr on the the ReadWriteCsv. */
    inline ReadWriteCsv* const p_readWriteCsv() const { return p_data_;}
    /** release the ReadWriteCsv. It will be freed by the user. */
    inline void release() { p_data_ =0;}

  protected:
    /** ptr on the ReadWriteCsv containing the data. */
    ReadWriteCsv* p_data_;
};

} // namespace STK

#endif /*STK_EXPORTCSV_H*/
