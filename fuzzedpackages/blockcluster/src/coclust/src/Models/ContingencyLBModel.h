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

/** @file ContingencyLBModel.h
 *  @brief Declares concrete model class ContingencyLBModel_mu_i_nu_j for Contingency Data.
 **/

#ifndef CONTINGENCYLBMODEL_H_
#define CONTINGENCYLBMODEL_H_

/** @brief Concrete model class for contingency data-sets.
 * This model assumes that the row and column effects are unknown.
 *
 */
#include "ICoClustModel.h"

class ContingencyLBModel: public ICoClustModel
{
  public:
    ContingencyLBModel( MatrixReal const& m_Dataij
                      , ModelParameters const& Mparam
                      );
    ContingencyLBModel( MatrixReal const& m_Dataij
                      , VectorInt const & rowlabels
                      , VectorInt const & collabels
                      , ModelParameters const& Mparam
                      );
    virtual ~ContingencyLBModel(){};
    virtual ContingencyLBModel* clone(){return new ContingencyLBModel(*this);}

    virtual void logSumRows(MatrixReal & m_sum);
    virtual void logSumCols(MatrixReal & m_sum);

//    virtual bool cemInitStep();
//    virtual bool emInitStep();
//    virtual bool randomInitStep();

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

    /**Return Poisson Parameters ContingencyLBModel::m_Gammakl_*/
    MatrixReal const& gamma() const;

  protected:
    //Variables involved in Poisson model
    MatrixReal const& m_Dataij_;
    MatrixReal m_ClusterDataij_;
    STK::Real DataSum_;
    MatrixReal m_Gammakl_, m_Gammaklold_;
    MatrixReal m_Gammakl1_, m_Gammakl1old_;
    MatrixReal m_Gammakltemp_;
    VectorReal v_Ui_,v_Vj_;
    MatrixReal m_Ykl_;

    //M-steps
    virtual void mStepRows();
    virtual void mStepCols();

    /** Compute m_Vjk_ array for all models */
    virtual void computeVjk();
    /** Compute m_Uil_ array for all models */
    virtual void computeUil();
};

inline MatrixReal const& ContingencyLBModel::gamma() const
{ return m_Gammakl_;}

inline void ContingencyLBModel::mStepRows()
{
  mSteplogPiek();
  m_Gammakl_ = (m_Tik_.transpose()*m_Uil_)/(v_Tk_* v_Rl_.transpose());
}

inline void ContingencyLBModel::mStepCols()
{
  mSteplogRhol();
  m_Gammakl_ = (m_Vjk_.transpose()*m_Rjl_)/(v_Tk_* v_Rl_.transpose());
}

inline void ContingencyLBModel::computeUil()
{ m_Uil_= m_Dataij_*m_Rjl_;}

inline void ContingencyLBModel::computeVjk()
{ m_Vjk_ = m_Dataij_.transpose()*m_Tik_;}


#endif /* CONTINGENCYLBMODEL_H_ */
