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

#ifndef CATEGORICALLBMODEL_H_
#define CATEGORICALLBMODEL_H_

/**@file CategoricalLBModel.h
 * @brief
 */
#include "ICoClustModel.h"

class CategoricalLBModel:public ICoClustModel
{
  public:
    CategoricalLBModel( MatrixInt const& m_Dataij
                      , ModelParameters const& Mparam
                      , STK::Real a=1, STK::Real b=1
                      );
    CategoricalLBModel( MatrixInt const& m_Dataij
                      , VectorInt const & rowlabels
                      , VectorInt const & collabels
                      , ModelParameters const& Mparam
                      , STK::Real a=1, STK::Real b=1
                      );
    virtual CategoricalLBModel* clone(){return new CategoricalLBModel(*this);}

    virtual void logSumRows(MatrixReal & m_sum);
    virtual void logSumCols(MatrixReal & m_sum);

//    virtual bool cemInitStep();
//    virtual bool emInitStep();
//    virtual bool randomInitStep();

    virtual void mStepFull(){};
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
    virtual STK::Real iclCriteriaValue();
    virtual void consoleOut();

    virtual void saveThetaInit();
    virtual void modifyTheta();
    virtual void copyTheta();

    MatrixInt const& arrangedDataClusters();
    inline const std::vector<MatrixReal>& mean(){return m3_Alphahkl_;}
    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;

    virtual ~CategoricalLBModel();

  protected:
    /// hyper-parameters
    STK::Real a_,b_;
    MatrixInt const& m_Dataij_;
    MatrixInt m_ClusterDataij_;
    /// Vector of the summed matrices
    VectorReal v_Ui_,v_Vj_;
    int r_; //number of categories
    /// parameters sets
    std::vector<MatrixReal> m3_Alphahkl_,m3_Alphahklold_;
    std::vector<MatrixReal> m3_Alphahkl1_, m3_Alphahkl1old_;
    std::vector<MatrixReal> m3_Alphahkltemp_, m3_logAlphahkl_;
    // TODO replace binary matrices by sparse matrices
    /// binary matrices
    std::vector<MatrixBinary> m3_Yhij_,m3_Yijh_,m3_Yjih_;//different ways to store data

    virtual void mStepRows();
    virtual void mStepCols();
    /** compute logRhol during the m-step */
    virtual void mSteplogRhol();
    /** compute logPiek during the m-step */
    virtual void mSteplogPiek();

    virtual void mGibbsStepRows();
    virtual void mGibbsStepCols();

    /** Compute m_Vjk_ array for all models */
    virtual void computeVjk();
    /** Compute m_Uil_ array for all models */
    virtual void computeUil();

  private:
    void initializeStorages();
};

inline void CategoricalLBModel::computeUil(){}
inline void CategoricalLBModel::computeVjk(){}

#endif /* CATEGORICALLBMODEL_H_ */
