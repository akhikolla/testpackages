/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2015  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
 */


/** @file ContinuousLBModel.h
 *  @brief Declares concrete model class ContinuousLBModel derived from ICoClustModel.
 **/

#ifndef CONTINUOUSLBMODEL_H_
#define CONTINUOUSLBMODEL_H_

#include "ICoClustModel.h"

/** @brief Concrete model class for continuous data.
 *  This class does not presume equal variance among various  co-clusters.
 */
class ContinuousLBModel: public ICoClustModel
{
  public:
    ContinuousLBModel( MatrixReal const& m_Dataij
                     , ModelParameters const& Mparam
                     );
    ContinuousLBModel( MatrixReal const& m_Dataij
                     , VectorInt const & rowlabels
                     , VectorInt const & collabels
                     , ModelParameters const& Mparam
                     );
    virtual ~ContinuousLBModel(){};
    virtual ContinuousLBModel* clone(){return new ContinuousLBModel(*this);}

    virtual void logSumRows(MatrixReal & m_sum);
    virtual void logSumCols(MatrixReal & m_sum);

  public:
    virtual void mStepFull();
    virtual bool emRows();
    virtual bool cemRows();
    virtual bool emCols();
    virtual bool cemCols();
    virtual bool semRows();
    virtual bool semCols();
    virtual bool GibbsRows();
    virtual bool GibbsCols();

    virtual STK::Real computeLnLikelihood();
    virtual bool initStopCriteria();
    virtual void parameterStopCriteria();
    virtual void consoleOut();

    virtual void saveThetaInit();
    virtual void modifyTheta();
    virtual void copyTheta();

    MatrixReal const& arrangedDataClusters();
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;

    /**Return various co-clusters mean ContinuousLBModel::m_Mukl_*/
    MatrixReal const& mean() const;
    /**Return various co-clusters Sigma ContinuousLBModel::m_Sigma2kl_*/
    MatrixReal const& sigma2() const;

  protected:
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    MatrixReal m_Dataij2_;
    MatrixReal m_Mukl_;
    MatrixReal m_Sigma2kl_, m_Sigma2kltemp_;
    MatrixReal m_Muklold1_,m_Muklold2_, m_Mukltemp_;
    MatrixReal m_Vjk2_, m_Uil2_;

    //M-steps
    virtual void mStepRows();
    virtual void mStepCols();
    /** Compute m_Vjk_ array for all models */
    virtual void computeVjk();
    /** Compute m_Uil_ array for all models */
    virtual void computeUil();
};

inline MatrixReal const& ContinuousLBModel::mean() const
{ return m_Mukl_;}

inline MatrixReal const& ContinuousLBModel::sigma2() const
{ return m_Sigma2kl_;}

inline void ContinuousLBModel::mStepRows()
{
  mSteplogPiek();
  MatrixReal m_trkl = v_Tk_*v_Rl_.transpose();
  m_Mukl_           = (m_Tik_.transpose()*m_Uil_)/m_trkl;
  m_Sigma2kl_       = (m_Tik_.transpose()*m_Uil2_)/m_trkl - m_Mukl_.square();
}

inline void ContinuousLBModel::mStepCols()
{
  mSteplogRhol();
  MatrixReal m_trkl = v_Tk_*v_Rl_.transpose();
  m_Mukl_           = (m_Vjk_.transpose()*m_Rjl_)/m_trkl;
  m_Sigma2kl_       = (m_Vjk2_.transpose()*m_Rjl_)/m_trkl - m_Mukl_.square();
}

inline void ContinuousLBModel::computeUil()
{
  m_Uil_ = m_Dataij_*m_Rjl_;
  m_Uil2_ = m_Dataij2_*m_Rjl_;
}

inline void ContinuousLBModel::computeVjk()
{
  m_Vjk_ = m_Dataij_.transpose()*m_Tik_;
  m_Vjk2_ = m_Dataij2_.transpose()*m_Tik_;
}

#endif /* CONTINUOUSLBMODEL_H_ */
