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


/** @file BinaryLBModel.h
 *  @brief Declares concrete model class BinaryLBModel derived from ICoClustModel.
 **/

#ifndef BinaryLBModel_H_
#define BinaryLBModel_H_

/** @brief Concrete model class for binary data.
 * This class assumes unequal dispersion among various co-clusters.
 *
 */
#include <limits.h>
#include <iostream>
#include "../typedefs/typedef.h"
#include "../Models/ICoClustModel.h"

#ifndef RPACKAGE
#include "../../CImg/CImg.h"
#endif

class BinaryLBModel : public ICoClustModel
{
  public:
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b bayesian hyperparameters
     * */
    BinaryLBModel( MatrixBinary const&  m_Dataij
                 , ModelParameters const& Mparam
                 , STK::Real a=1, STK::Real b=1
                 );
    /**Constructor for unsupervised co-clustering
     * @param m_Dataij a constant reference on the data set.
     * @param rowlabels various labels for rows (-1  for unknown row label)
     * @param collabels various labels for columns (-1 for unknown column label)
     * @param Mparam A constant reference to various ModelParameters.
     * @param a,b Bayesian hyperparameters
     * */
    BinaryLBModel( MatrixBinary const&  m_Dataij
                 , VectorInt const& rowlabels
                 , VectorInt const& collabels
                 , ModelParameters const& Mparam
                 , STK::Real a=1, STK::Real b=1
                 );

    /** cloning */
    inline virtual BinaryLBModel* clone(){return new BinaryLBModel(*this);}

    virtual void logSumRows(MatrixReal & _m_sik);
    virtual void logSumCols(MatrixReal & _m_sjl);

//    virtual bool cemInitStep();
//    virtual bool emInitStep();
//    virtual bool randomInitStep();

    virtual void mStepFull();

    virtual bool emRows();
    virtual bool cemRows();
    virtual bool semRows();
    virtual bool GibbsRows();

    virtual bool emCols();
    virtual bool cemCols();
    virtual bool semCols();
    virtual bool GibbsCols();

    virtual STK::Real computeLnLikelihood();
    virtual bool initStopCriteria();
    virtual void parameterStopCriteria();
    virtual STK::Real iclCriteriaValue();
    virtual void finalizeOutput();
    virtual void consoleOut();

    virtual void saveThetaInit();
    virtual void copyTheta();
    virtual void modifyTheta();

    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const;

    MatrixBinary const&  arrangedDataClusters();
    /**Return class mean BinaryLBModel::m_akl_ for all the blocks (co-clusters)*/
    MatrixBinary const&  mean() const;
    /**Return Class despersion BinaryLBModel::m_epsilonkl_ for all the blocks (co-clusters) */
    MatrixReal const& dispersion() const;
    /**Destructor*/
    inline virtual ~BinaryLBModel(){};

#ifndef RPACKAGE
    /**This function display co-clusters for Binary data.*/
    void displayCluster();
#endif

  protected:
    //Variables involved in Bernouilli model
    STK::Real a_,b_;//hyper-parameters
    MatrixBinary const&  m_Dataij_;
    MatrixBinary m_ClusterDataij_;
    MatrixReal m_Alphakl_, m_Alphaklold_;
    MatrixReal m_Alphakl1_, m_Alphakl1old_;
    MatrixReal m_Alphakltemp_;
    MatrixBinary m_akl_;
    MatrixReal m_epsilonkl_, m_epsilonkltemp_;

    virtual void mStepRows();
    virtual void mStepCols();
    /** compute logRhol during the m-step */
    virtual void mSteplogRhol();
    /** compute logPiek during the m-step */
    virtual void mSteplogPiek();

    void mGibbsStepRows();
    void mGibbsStepCols();

    /** Compute m_Vjk_ array for all models */
    virtual void computeVjk();
    /** Compute m_Uil_ array for all models */
    virtual void computeUil();
};

inline MatrixBinary const&  BinaryLBModel::mean() const
{ return m_akl_;}

inline MatrixReal const& BinaryLBModel::dispersion() const
{ return m_epsilonkl_;}

inline void BinaryLBModel::computeUil()
{ m_Uil_ = m_Dataij_.cast<STK::Real>()*m_Rjl_;}

inline void BinaryLBModel::computeVjk()
{ m_Vjk_ = m_Dataij_.cast<STK::Real>().transpose()*m_Tik_;}

#endif /* BinaryLBModel_H_ */
