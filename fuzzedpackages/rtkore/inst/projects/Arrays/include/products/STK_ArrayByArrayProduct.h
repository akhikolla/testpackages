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
 * Project:  stkpp::Arrays
 * created on: 30 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayByArrayProduct.h
 *  @brief In this file we implement the General Array by Array product.
 **/


#ifndef STK_ARRAYBYARRAYPRODUCT_H
#define STK_ARRAYBYARRAYPRODUCT_H

namespace STK
{

namespace hidden
{
//template<typename Lhs, typename Rhs, typename Result, bool Orient_> struct bp;
template<typename Lhs, typename Rhs, typename Result> struct BlockByPanel;
template<typename Lhs, typename Rhs, typename Result> struct PanelByBlock;
template<typename Lhs, typename Rhs, typename Result> struct BlockPanelProduct;
template<typename Lhs, typename Rhs, typename Result> struct PanelBlockProduct;

/** @ingroup hidden
 *  Methods to use for C=AB with A divided in blocks and B divided in panels.
 *  The structure BlockByPanel use data cache and contains only static method
 *  and typedef and should normally not be used directly.
 *
 *  @sa BlockPanelProduct
 **/
template<typename Lhs, typename Rhs, typename Result>
struct BlockByPanel
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  typedef hidden::CopySubArrayImpl<Lhs, Type> CopyLhsImpl;
  typedef hidden::CopySubArrayImpl<Rhs, Type> CopyRhsImpl;

  /** Main method for matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
#ifdef STK_ARRAYS_DEBUG
     stk_cout << _T("Entering BlockByPanel::run()\n");
#endif
    // compute dimensions
    int nbInnerLoop = lhs.sizeCols()/blockSize_; // = rhs.sizeRows()/blockSize_;
    int nbBlocks    = lhs.sizeRows()/blockSize_;
    int nbPanels    = rhs.sizeCols()/panelSize_;
    // remaining sizes in the matrices
    int pSize = rhs.sizeCols() - panelSize_*nbPanels;
    int bSize = lhs.sizeRows() - blockSize_*nbBlocks;
    int tSize = lhs.sizeCols() - blockSize_*nbInnerLoop;
              // = rhs.sizeRows() -  rhs.sizeRows()/blockSize_
    int iLastRow = lhs.beginRows() + nbBlocks * blockSize_;
    int jLastCol = rhs.beginCols() + nbPanels * panelSize_;
    int kLastPos = lhs.beginCols() + blockSize_ * nbInnerLoop;
    // start
    if (nbInnerLoop)
    {
      // create panels and blocks
      Panel<Type>* tabPanel = new Panel<Type>[nbPanels+1];
      Block<Type>* tabBlock = new Block<Type>[nbBlocks+1];
      // start blocks by panel
      for (int k = 0,kPos = lhs.beginCols(); k<nbInnerLoop; ++k, kPos += blockSize_)
      {
        // data caching
        for (int i = 0, iRow = lhs.beginRows(); i<nbBlocks; ++i, iRow += blockSize_)
        { CopyLhsImpl::arrayToBlock( lhs, tabBlock[i], iRow, kPos);}
        CopyLhsImpl::arrayToBlock( lhs, tabBlock[nbBlocks], iLastRow, kPos, bSize);
        for (int j = 0, jCol = rhs.beginCols(); j<nbPanels; ++j, jCol += panelSize_)
        { CopyRhsImpl::arrayToPanel( rhs, tabPanel[j], kPos, jCol);}
        CopyRhsImpl::arrayToPanel( rhs, tabPanel[nbPanels], kPos, jLastCol, pSize);
        // block by panel products
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i<nbBlocks; ++i)
        {
          int iRow = lhs.beginRows()+ i * blockSize_;
          for (int j = 0, jCol = rhs.beginCols(); j<nbPanels; ++j, jCol += panelSize_)
          { multBlockByPanel( tabBlock[i], tabPanel[j], res, iRow, jCol);}
        }
        for (int i = 0, iRow = lhs.beginRows(); i<nbBlocks; ++i,iRow += blockSize_)
        { multBlockByPanel( tabBlock[i], tabPanel[nbPanels], res, iRow, jLastCol, pSize);}
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j<nbPanels; ++j)
        { multBlockByPanel( tabBlock[nbBlocks], tabPanel[j], res, iLastRow, rhs.beginCols() + j * panelSize_, panelSize_, bSize);}
        multBlockByPanel( tabBlock[nbBlocks], tabPanel[nbPanels], res, iLastRow, jLastCol, pSize, bSize);
      } // InnerLoop
      delete[] tabBlock;
      delete[] tabPanel;
    } // if IneerLoop
    // treat the remaining rows, columns
    switch (tSize)
    {
      case 0: break;
      case 1: MultCoeff::mult1Outer(lhs, rhs, res, kLastPos); break;
      case 2: MultCoeff::mult2Outer(lhs, rhs, res, kLastPos); break;
      case 3: MultCoeff::mult3Outer(lhs, rhs, res, kLastPos); break;
      default:break;
    }
  }
  /** Default dimension */
  static void multBlockByPanel( Block<Type> const& block, Panel<Type> const& panel, Result& res
                              , int iRow, int jCol)
  {
    for (int j=0; j<panelSize_; ++j)
    {
      res.elt(iRow  ,jCol+j) += panel[j*blockSize_]    * block[0]
                              + panel[j*blockSize_+ 1] * block[1]
                              + panel[j*blockSize_+ 2] * block[2]
                              + panel[j*blockSize_+ 3] * block[3];
      res.elt(iRow+1,jCol+j) += panel[j*blockSize_]    * block[4]
                              + panel[j*blockSize_+ 1] * block[5]
                              + panel[j*blockSize_+ 2] * block[6]
                              + panel[j*blockSize_+ 3] * block[7];
      res.elt(iRow+2,jCol+j) += panel[j*blockSize_]    * block[8]
                              + panel[j*blockSize_+ 1] * block[9]
                              + panel[j*blockSize_+ 2] * block[10]
                              + panel[j*blockSize_+ 3] * block[11];
      res.elt(iRow+3,jCol+j) += panel[j*blockSize_]    * block[12]
                              + panel[j*blockSize_+ 1] * block[13]
                              + panel[j*blockSize_+ 2] * block[14]
                              + panel[j*blockSize_+ 3] * block[15];
    }
  }
  /** with panel size given */
  static void multBlockByPanel( Block<Type> const& block, Panel<Type> const& panel, Result& res
                              , int iRow, int jCol, int pSize)
  {
    for (int j=0; j<pSize; ++j)
    {
      res.elt(iRow  ,jCol+j) += panel[j*blockSize_]    * block[0]
                              + panel[j*blockSize_+ 1] * block[1]
                              + panel[j*blockSize_+ 2] * block[2]
                              + panel[j*blockSize_+ 3] * block[3];
      res.elt(iRow+1,jCol+j) += panel[j*blockSize_]    * block[4]
                              + panel[j*blockSize_+ 1] * block[5]
                              + panel[j*blockSize_+ 2] * block[6]
                              + panel[j*blockSize_+ 3] * block[7];
      res.elt(iRow+2,jCol+j) += panel[j*blockSize_]    * block[8]
                              + panel[j*blockSize_+ 1] * block[9]
                              + panel[j*blockSize_+ 2] * block[10]
                              + panel[j*blockSize_+ 3] * block[11];
      res.elt(iRow+3,jCol+j) += panel[j*blockSize_]    * block[12]
                              + panel[j*blockSize_+ 1] * block[13]
                              + panel[j*blockSize_+ 2] * block[14]
                              + panel[j*blockSize_+ 3] * block[15];
    }
  }
  /** with panel size given */
  static void multBlockByPanel( Block<Type> const& block, Panel<Type> const& panel, Result& res
                              , int iRow, int jCol, int pSize, int bSize)
  {
    for (int i=0; i<bSize; ++i)
      for (int j=0; j<pSize; ++j)
      { res.elt(iRow+i,jCol+j) += panel[j*blockSize_]   * block[i*blockSize_]
                                + panel[j*blockSize_+1] * block[i*blockSize_+1]
                                + panel[j*blockSize_+2] * block[i*blockSize_+2]
                                + panel[j*blockSize_+3] * block[i*blockSize_+3];
      }
  }
}; // struct BlockByPanel


/** @ingroup hidden
 *  @brief Methods to use for C=AB with A divided in panels and B divided in blocks.
 *  The structure PanelByBlock use data cache and contains only static method and typedef
 *  and should normally not be used directly.
 *  @sa PanelBlockProduct
 **/
template<typename Lhs, typename Rhs, typename Result>
struct PanelByBlock
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  typedef hidden::CopySubArrayImpl<Lhs, Type> CopyLhsImpl;
  typedef hidden::CopySubArrayImpl<Rhs, Type> CopyRhsImpl;
  /** Main method for Matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
#ifdef STK_ARRAYS_DEBUG
     stk_cout << _T("Entering PanelByBlock::run()\n");
#endif
    // compute dimensions
    int nbInnerLoop = lhs.sizeCols()/blockSize_; // = rhs.sizeRows()/blockSize_;
    int nbBlocks    = rhs.sizeCols()/blockSize_;
    int nbPanels    = lhs.sizeRows()/panelSize_;

    // remaining sizes in the matrices
    int pSize = lhs.sizeRows() - panelSize_ * nbPanels;
    int bSize = rhs.sizeCols() - blockSize_ * nbBlocks;
    int tSize = lhs.sizeCols() - blockSize_ * nbInnerLoop;

    // index of the remaining positions
    int jLastCol = rhs.beginCols() + blockSize_ * nbBlocks;
    int iLastRow = lhs.beginRows() + panelSize_ * nbPanels;
    int kLastPos = rhs.beginRows() + blockSize_ * nbInnerLoop;

    if (nbInnerLoop)
    {
      // create panels
      Panel<Type>* tabPanel = new Panel<Type>[nbPanels+1];
      Block<Type>* tabBlock = new Block<Type>[nbBlocks+1];
      // start blocks by panel
      for (int k = 0, kPos = rhs.beginRows(); k<nbInnerLoop; ++k,kPos += blockSize_)
      {
        // data caching
        for (int i = 0, iRow= lhs.beginRows(); i<nbPanels; ++i, iRow+= panelSize_)
        { CopyLhsImpl::arrayToPanelByCol( lhs, tabPanel[i], iRow, kPos);}
        CopyLhsImpl::arrayToPanelByCol( lhs, tabPanel[nbPanels], iLastRow, kPos, pSize);
        // get blocks
        for (int j = 0, jCol = rhs.beginCols(); j<nbBlocks; ++j, jCol+=blockSize_)
        { CopyRhsImpl::arrayToBlockByCol( rhs, tabBlock[j], kPos, jCol);}
        CopyRhsImpl::arrayToBlockByCol( rhs, tabBlock[nbBlocks], kPos, jLastCol, bSize);
        // perform the products panels * blocks
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j<nbBlocks; ++j)
        {
          int jCol = rhs.beginCols() + j* blockSize_;
          for (int i = 0,iRow = lhs.beginRows(); i<nbPanels; ++i, iRow += panelSize_)
          { multPanelByBlock( tabPanel[i], tabBlock[j], res, iRow, jCol);}
          multPanelByBlock( tabPanel[nbPanels], tabBlock[j], res, iLastRow, jCol, pSize);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i<nbPanels; ++i)
        { multPanelByBlock( tabPanel[i],  tabBlock[nbBlocks], res, lhs.beginRows() + i * panelSize_, jLastCol, panelSize_, bSize);}
        multPanelByBlock( tabPanel[nbPanels],  tabBlock[nbBlocks], res, iLastRow, jLastCol, pSize, bSize);
      } // k loop
      delete[] tabPanel;
      delete[] tabBlock;
    }
    // treat the remaining rows, columns
    switch (tSize)
    {
      case 0: break;
      case 1: MultCoeff::mult1Outer(lhs, rhs, res, kLastPos); break;
      case 2: MultCoeff::mult2Outer(lhs, rhs, res, kLastPos); break;
      case 3: MultCoeff::mult3Outer(lhs, rhs, res, kLastPos); break;
      default:break;
    }
  }
  /** Default dimension */
  static void multPanelByBlock( Panel<Type> const& panel, Block<Type> const& block, Result& res
                              , int iRow, int jCol)
  {
    for (int i=0; i<panelSize_; ++i)
    {
      res.elt(iRow+i,jCol)   += panel[i*blockSize_]    * block[0]
                              + panel[i*blockSize_+ 1] * block[1]
                              + panel[i*blockSize_+ 2] * block[2]
                              + panel[i*blockSize_+ 3] * block[3];
      res.elt(iRow+i,jCol+1) += panel[i*blockSize_]    * block[4]
                              + panel[i*blockSize_+ 1] * block[5]
                              + panel[i*blockSize_+ 2] * block[6]
                              + panel[i*blockSize_+ 3] * block[7];
      res.elt(iRow+i,jCol+2) += panel[i*blockSize_]    * block[8]
                              + panel[i*blockSize_+ 1] * block[9]
                              + panel[i*blockSize_+ 2] * block[10]
                              + panel[i*blockSize_+ 3] * block[11];
      res.elt(iRow+i,jCol+3) += panel[i*blockSize_]    * block[12]
                              + panel[i*blockSize_+ 1] * block[13]
                              + panel[i*blockSize_+ 2] * block[14]
                              + panel[i*blockSize_+ 3] * block[15];
    }
  }
  static void multPanelByBlock( Panel<Type> const& panel, Block<Type> const& block, Result& res
                              , int iRow, int jCol, int pSize)
  {
    for (int i=0; i<pSize; ++i)
    {
      res.elt(iRow+i,jCol)   += panel[i*blockSize_]    * block[0]
                              + panel[i*blockSize_+ 1] * block[1]
                              + panel[i*blockSize_+ 2] * block[2]
                              + panel[i*blockSize_+ 3] * block[3];
      res.elt(iRow+i,jCol+1) += panel[i*blockSize_]    * block[4]
                              + panel[i*blockSize_+ 1] * block[5]
                              + panel[i*blockSize_+ 2] * block[6]
                              + panel[i*blockSize_+ 3] * block[7];
      res.elt(iRow+i,jCol+2) += panel[i*blockSize_]    * block[8]
                              + panel[i*blockSize_+ 1] * block[9]
                              + panel[i*blockSize_+ 2] * block[10]
                              + panel[i*blockSize_+ 3] * block[11];
      res.elt(iRow+i,jCol+3) += panel[i*blockSize_]    * block[12]
                              + panel[i*blockSize_+ 1] * block[13]
                              + panel[i*blockSize_+ 2] * block[14]
                              + panel[i*blockSize_+ 3] * block[15];
    }
  }
  /** with panel size and block size dimension given */
  static void multPanelByBlock( Panel<Type> const& panel, Block<Type> const&  block, Result& res
                              , int iRow, int jCol, int pSize, int bSize)
  {
    for (int i=0; i<pSize; ++i)
      for (int j=0; j<bSize; ++j)
        res.elt(iRow+i,jCol+j) += panel[i*blockSize_]   * block[j*blockSize_]
                                + panel[i*blockSize_+1] * block[j*blockSize_+1]
                                + panel[i*blockSize_+2] * block[j*blockSize_+2]
                                + panel[i*blockSize_+3] * block[j*blockSize_+3];
  }
}; // struct PanelByBlock


/** @ingroup hidden
 *  @brief Methods to use for C=AB with A divided in blocks and B divided in panels.
 *  The structure BlockPanelProduct contains only static methods and typedef
 *  and should normally not be used directly.
 *  @sa BlockByPanel
 **/
template<typename Lhs, typename Rhs, typename Result>
struct BlockPanelProduct
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;

  /** Main method for Matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
#ifdef STK_DEBUG
     stk_cout << _T("Entering BlockPanelProduct::run()\n");
#endif
    // compute number of loops
    int nbInnerLoop = lhs.sizeCols()/blockSize_;
    int nbBlocks    = lhs.sizeRows()/blockSize_;
    int nbPanels    = rhs.sizeCols()/panelSize_;

    // remaining sizes in the matrices
    int bSize = lhs.sizeRows() - blockSize_ * nbBlocks;
    int pSize = rhs.sizeCols() - panelSize_ * nbPanels;
    int tSize = lhs.sizeCols() - blockSize_ * nbInnerLoop;

    // remaining positions in the matrices
    int kLastPos = lhs.beginCols() + blockSize_ * nbInnerLoop;

    // start blocks by panel
    for (int k = 0; k<nbInnerLoop; ++k)
    {
      TRange<blockSize_> innerRange(lhs.beginCols() + k * blockSize_, blockSize_);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i<nbBlocks; ++i)
      {
        TRange<blockSize_> rowRange(lhs.beginRows() + i * blockSize_, blockSize_);
        for (int j = 0; j<nbPanels; ++j)
        {
          TRange<panelSize_> colRange(rhs.beginCols() + j * panelSize_, panelSize_);
          multBlockByPanel( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
        }
        // remaining incomplete panels
        Range colRange(rhs.beginCols() + nbPanels * panelSize_, pSize);
        multBlockByPanel( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
      }
      Range rowRange(lhs.beginRows() + nbBlocks * blockSize_, bSize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int j = 0; j<nbPanels; ++j)
      {
        TRange<panelSize_> colRange(res.beginCols() + j * panelSize_, panelSize_);
        multBlockPartByPanel( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
      }
      // remaining incomplete panels
      Range colRange(res.beginCols() + nbPanels * panelSize_, pSize);
      multBlockPartByPanel( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
    } // InnerLoop
    // treat the remaining rows/columns not in outer loop k
    switch (tSize)
    {
      case 1: MultCoeff::mult1Outer(lhs, rhs, res, kLastPos); break;
      case 2: MultCoeff::mult2Outer(lhs, rhs, res, kLastPos); break;
      case 3: MultCoeff::mult3Outer(lhs, rhs, res, kLastPos); break;
      default:break;
    }
  }
  /** matrix product between a lhs block and a rhs panel */
  template<class SubLhs, class SubRhs>
  static void multBlockByPanel( SubLhs const& lhs, SubRhs const& rhs, Result& res)
  {
    const int k= rhs.beginRows();
    for (int j=rhs.beginCols(); j<rhs.endCols(); ++j)
    {
      int i = lhs.beginRows();
      res.elt(i,j) += lhs(i  ,k)   * rhs(k  ,j) + lhs(i  ,k+1) * rhs(k+1,j)
                    + lhs(i  ,k+2) * rhs(k+2,j) + lhs(i  ,k+3) * rhs(k+3,j);
      ++i;
      res.elt(i,j) += lhs(i,k)   * rhs(k  ,j) + lhs(i,k+1) * rhs(k+1,j)
                    + lhs(i,k+2) * rhs(k+2,j) + lhs(i,k+3) * rhs(k+3,j);
      ++i;
      res.elt(i,j) += lhs(i,k)   * rhs(k  ,j) + lhs(i,k+1) * rhs(k+1,j)
                    + lhs(i,k+2) * rhs(k+2,j) + lhs(i,k+3) * rhs(k+3,j);
      ++i;
      res.elt(i,j) += lhs(i,k)   * rhs(k  ,j) + lhs(i,k+1) * rhs(k+1,j)
                    + lhs(i,k+2) * rhs(k+2,j) + lhs(i,k+3) * rhs(k+3,j);
    }
  }
  /** matrix product between a lhs block and a rhs panel */
  template<class SubLhs, class SubRhs>
  static void multBlockPartByPanel( SubLhs const& lhs, SubRhs const& rhs, Result& res)
  {
    const int k= rhs.beginRows();
    for (int j=rhs.beginCols(); j<rhs.endCols(); ++j)
      for(int i=lhs.beginRows(); i<lhs.endRows(); ++i)
        res.elt(i,j) += lhs(i,k  ) * rhs(k  ,j) + lhs(i,k+1) * rhs(k+1,j)
                      + lhs(i,k+2) * rhs(k+2,j) + lhs(i,k+3) * rhs(k+3,j);
  }
}; // struct BlockPanelProduct

/** @ingroup hidden
 *  @brief Methods to use for C=AB with A divided in panels and B divided in blocks.
 *  The structure PanelBlockProduct contains only static methods and typedef
 *  and should normally not be used directly.
 *  @sa PanelByBlock
 *
 **/
template<typename Lhs, typename Rhs, typename Result>
struct PanelBlockProduct
{
  typedef typename Result::Type Type;
  typedef hidden::MultImpl<Type> Cmult;
  typedef hidden::MultCoefImpl<Lhs, Rhs, Result> MultCoeff;

  /** Main method for Matrices multiplication implementation.
   *  @note res have been resized and initialized to zero outside this method.
   **/
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res)
  {
#ifdef STK_ARRAYS_DEBUG
     stk_cout << _T("Entering PanelBlockProduct::run()\n");
#endif
    // compute dimensions
    int nbInnerLoop = lhs.sizeCols()/blockSize_; // = rhs.sizeRows()/blockSize_;
    int nbBlocks    = rhs.sizeCols()/blockSize_;
    int nbPanels    = lhs.sizeRows()/panelSize_;

    // remaining sizes in the matrices
    int bSize = rhs.sizeCols() - blockSize_ * nbBlocks;
    int pSize = lhs.sizeRows() - panelSize_ * nbPanels;
    int tSize = lhs.sizeCols() - blockSize_ * nbInnerLoop;

    // index of the remaining positions
    int kLastPos = lhs.beginCols() + blockSize_ * nbInnerLoop;

    // start blocks by panel
    for (int k = 0; k<nbInnerLoop; ++k)
    {
      TRange<blockSize_> innerRange(lhs.beginCols() + k * blockSize_, blockSize_);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int j = 0; j<nbBlocks; ++j)
      {
        TRange<blockSize_> colRange(rhs.beginCols() + j * blockSize_, blockSize_);
        for (int i = 0; i<nbPanels; ++i)
        {
          TRange<panelSize_> rowRange(lhs.beginRows() + i * panelSize_, panelSize_);
          multPanelByBlock( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
        }
        // remaining incomplete
        Range rowRange(lhs.beginRows() + nbPanels * panelSize_, pSize);
        multPanelByBlock( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
      }
      Range colRange(rhs.beginCols() + nbBlocks * blockSize_, bSize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i<nbPanels; ++i)
      {
        TRange<panelSize_> rowRange(lhs.beginRows() + i * panelSize_, panelSize_);
        multPanelByBlockPart( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
      }
      // remaining incomplete panels
      Range rowRange(lhs.beginRows() + nbPanels * panelSize_, pSize);
      multPanelByBlockPart( lhs.sub(rowRange, innerRange), rhs.sub(innerRange, colRange), res);
    } // InnerLoop
    // treat the remaining rows/columns not in outer loop k
    switch (tSize)
    {
      case 0: break;
      case 1: MultCoeff::mult1Outer(lhs, rhs, res, kLastPos); break;
      case 2: MultCoeff::mult2Outer(lhs, rhs, res, kLastPos); break;
      case 3: MultCoeff::mult3Outer(lhs, rhs, res, kLastPos); break;
      default:break;
    }
  }
  /** Default dimension */
  template<class SubLhs, class SubRhs>
  static void multPanelByBlock( SubLhs const& lhs, SubRhs const& rhs, Result& res)
  {
    int const k= lhs.beginCols();
    for (int i=lhs.beginRows(); i<lhs.endRows(); ++i)
    {
      int j = rhs.beginCols();
      res.elt(i,j) += lhs(i, k  ) * rhs(k  , j) + lhs(i, k+1) * rhs(k+1, j)
                    + lhs(i, k+2) * rhs(k+2, j) + lhs(i, k+3) * rhs(k+3, j);
      ++j;
      res.elt(i,j) += lhs(i, k  ) * rhs(k  , j) + lhs(i, k+1) * rhs(k+1, j)
                    + lhs(i, k+2) * rhs(k+2, j) + lhs(i, k+3) * rhs(k+3, j);
      ++j;
      res.elt(i,j) += lhs(i, k  ) * rhs(k  , j) + lhs(i, k+1) * rhs(k+1, j)
                    + lhs(i, k+2) * rhs(k+2, j) + lhs(i, k+3) * rhs(k+3, j);
      ++j;
      res.elt(i,j) += lhs(i, k  ) * rhs(k  , j) + lhs(i, k+1) * rhs(k+1, j)
                    + lhs(i, k+2) * rhs(k+2, j) + lhs(i, k+3) * rhs(k+3, j);
    }
  }
  /** Default dimension */
  template<class SubLhs, class SubRhs>
  static void multPanelByBlockPart( SubLhs const& lhs, SubRhs const& rhs, Result& res)
  {
    const int k= lhs.beginCols();
    for (int i=lhs.beginRows(); i<lhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j<rhs.endCols(); ++j)
        res.elt(i,j) += lhs(i, k  ) * rhs(k  , j) + lhs(i, k+1) * rhs(k+1, j)
                      + lhs(i, k+2) * rhs(k+2, j) + lhs(i, k+3) * rhs(k+3, j);
  }
}; // struct PanelBlockProduct

} // namespace hidden

} // namespace STK

#endif /* STK_ARRAYBYARRAYPRODUCT_H */
