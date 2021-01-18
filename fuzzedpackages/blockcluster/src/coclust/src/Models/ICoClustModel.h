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


#ifndef ICOCLUSTMODEL_H_
#define ICOCLUSTMODEL_H_

/** @file ICoClustModel.h
 *  @brief This file declares the ICoClustModel abstract model class. All the concrete models
 *  classes are derived from this abstract class.
 **/

#include "../InputParameters/InputParameters.h"

/** @brief This is an abstract class which provides interface for various models.
 *  It provides interfaces for the most common functions in EM, CEM  and SEM algorithms as well
 *  as provide implementation for common functionalities.
 */
class ICoClustModel
{
  protected:
    /** Constructor
     *  @param Mparam ModelParameters
     **/
    ICoClustModel(ModelParameters const& Mparam);
    /** Constructor
     * @param Mparam ModelParameters
     * @param rowlabels Row clusters for each row (-1 for unknown cluster for each row)
     * @param collabels Column clusters for each column (-1 unknown cluster for each column)
     * */
    ICoClustModel( ModelParameters const& Mparam
                 , VectorInt const & rowlabels
                 , VectorInt const & collabels
                 );
  public:
    /** Destructor*/
    inline virtual ~ICoClustModel(){};
    // getters
    /** Get the likelihood value*/
    STK::Real likelihood() const;
    /** This function will return the row classification vector.
     * @return Row classification vector
     */
    VectorInt const& rowClassificationVector() const;
    /** This function will return the column classification vector.
     * @return COlumn classification vector
     */
    VectorInt const& columnClassificationVector() const;
    /** This function will return the row proportions of mixtures models.
     * @return Row proportions of mixing models.
     */
    VectorReal const& rowProportions() const;
    /** This function will return the column proportions of mixtures models.
     * @return column proportions of mixing models.
     */
    const VectorReal & colProportions() const;
    /** This function will return the posterior probabilities for each row.
     * @return Row posterior probabilities (Rows represent original rows and columns represent class probabilities)
     */
    const MatrixReal & rowPosteriorProb() const;
    /** This function will return the posterior probabilities for each column.
     * @return Column posterior probabilities ((Rows represent original columns and columns represent class probabilities))
     */
    const MatrixReal & colPosteriorProb() const;

    // pure virtual
    /** Cloning interface*/
    virtual ICoClustModel* clone() = 0;

    /** Interface function for CEM initialization . It will initialize model parameters
     * using CEM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool cemInitStep();
    /** Interface function for EM initialization . It will initialize model parameters
     * using EM algorithm.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool emInitStep();
    /** Interface function for Random Initialization. It will initialize model parameters
     * using Random initialization.
     * @return boolean value representing success or failure. By default, it will return
     * false if the derived class does not overwrite this method.
     */
    virtual bool randomInitStep();

    /** FUll M-step interface*/
    virtual  void mStepFull() = 0;
    /** Interface for calculating log sum for rows*/
    virtual void logSumRows(MatrixReal&) = 0;
    /** Interface for calculating log sum for columns*/
    virtual void logSumCols(MatrixReal&) = 0;
    /**Interface for EM Algorithm for rows*/
    virtual bool emRows() = 0;
    /**Interface for EM Algorithm for Columns*/
    virtual bool emCols() = 0;
    /**Interface for CEM Algorithm for rows*/
    virtual bool cemRows() = 0;
    /**Interface for CEM Algorithm for Columns*/
    virtual bool cemCols() = 0;
    /**Interface for SEM Algorithm for rows*/
    virtual bool semRows() = 0;
    /**Interface for SEM Algorithm for columns*/
    virtual bool semCols() = 0;
    /**Interface for Gibbs Algorithm for rows*/
    virtual bool GibbsRows() = 0;
    /**Interface for Gibbs Algorithm for columns*/
    virtual bool GibbsCols() = 0;
    /**Interface for maximization step for rows*/
    virtual void mStepRows() = 0;
    /**Interface for maximization step for columns*/
    virtual void mStepCols() = 0;

    /** @return the number of free parameters of the distribution of a block.*/
    virtual int nbFreeParameters() const = 0;
    /** Interface for calculating Stopping condition using percentage Change in Parameter values.
     *  This function will set the ICoClustModel::stopAlgo_ parameter to either true or false
     *  depending on whether the change in parameters is less than Mparam_.epsilon_ or not
     *  respectively.
     */
    virtual void parameterStopCriteria() = 0;
    /** Interface for calculating Stopping condition using percentage Change in Parameter values
     *  during initialization step.
     *  This function will return true or false depending on whether the change in parameters
     *  is less than Mparam_.epsilon_ or not respectively.
     */
    virtual bool initStopCriteria() = 0;
    /** Interface for output of model parameters on console. */
    virtual void consoleOut() = 0;

    /** Save current value of parameters during init step. */
    virtual void saveThetaInit() = 0;
    /** Modify the value of parameters during init/xem/XEM steps if better parameters are found.*/
    virtual void modifyTheta() = 0;
    /** Copy the value of parameters after init/xem/XEM steps.*/
    virtual void copyTheta() = 0;

    /** Estimate likelihood value*/
    virtual STK::Real computeLnLikelihood() = 0;

    /** Interface for finalizing the output. This function will allow the model to finalize
     * all the output parameters after the algorithm terminates. Do nothing by default.
     */
    virtual void finalizeOutput();
    /** This function will provide the current status of ICoClustModel::stopAlgo_ parameter.*/
    bool stopAlgo();
    /** Interface function for calculating ICL criteria value.
     * @return ICL criteria.
     */
    virtual STK::Real iclCriteriaValue();
    /** Interface for calculating Stopping condition using percentage Change in Likelihood. This function will set the
     * ICoClustModel::stopAlgo_ parameter to either true or false depending on whether the change
     * in Likelihood is less than Mparam_.epsilon_ or not respectively.
     */
    virtual void likelihoodStopCriteria();

    //common computation getter/setter functions
    /** @return ModelParameters */
    ModelParameters const& modelParameters() const;
    /** @param Mparam ModelParameters */
    void setModelParameters(ModelParameters& Mparam);
    /** @param Epsilon Value to be set as Mparam_.epsilon_ */
    void setEpsilon(STK::Real Epsilon);

    /**Get Error message*/
    inline std::string errorMsg(){ return Error_msg_;}
    /**Set Error message*/
    inline void setMsg(const std::string& s){ Error_msg_ = s;}
    /**Get status of empty cluster*/
    inline bool isEmptyCluster() const { return empty_cluster_;}
    /**set the value for Empty Cluster*/
    inline void setEmptyCluster(bool val){ empty_cluster_ = val;}

    /** Interface function for containers initialization. This method is called
     *  at creation of the object.
     **/
    void initializeStep();
    /** E-step for rows*/
    bool eStepRows();
    /** CE-step for rows*/
    bool ceStepRows();
    /** SE-step for rows*/
    bool seStepRows();

    /** E-step for columns*/
    bool eStepCols();
    /** CE-step for columns*/
    bool ceStepCols();
    /** SE-step for columns*/
    bool seStepCols();


    template<class T>
    void arrangedDataCluster(T&,const T&);

  protected:
    std::string  Error_msg_;

    //variables use in case of semi-supervised co-clustering
    STK::Array1D<int> UnknownLabelsRows_, UnknownLabelsCols_ ;
    STK::Array1D<std::pair<int,int> > knownLabelsRows_, knownLabelsCols_;

    /** Models parameters: dimensions, tolerances, etc. */
    ModelParameters Mparam_;
    /** current ln-likelihood value and maximal ln-likelihood value */
    STK::Real likelihood_,Lmax_;
    /** does current model has empty cluster ? */
    bool empty_cluster_;
    // row and columns posterior probabilities
    MatrixReal m_Tik_, m_Rjl_;
    MatrixReal m_Tiktemp_, m_Rjltemp_;
    // sum by column and by row of the posterior probabilities
    VectorReal v_Tk_,v_Rl_;
    // proportions and log proportions
    VectorReal v_Piek_, v_Rhol_;
    VectorReal v_logPiek_, v_logRhol_;
    VectorReal v_logPiektemp_, v_logRholtemp_;
    /** Row and column classification matrices respectively*/
    MatrixInt m_Zik_,m_Wjl_;
    /**Row and column classification vector*/
    VectorInt v_Zi_, v_Wj_;
    /** Intermediary array for computing rows/columns eStep, ceStep and semStep */
    MatrixReal m_Vjk_, m_Uil_;
    /** Boolean Variable used to identify whether the stopping condition
     *  is reached or not*/
    bool stopAlgo_;
    /** number of elements in dataij_ (Mparam.nbRow_ * Mparam.nbCol_) */
    STK::Real dimprod_;

    /** Compute m_Vjk_ array for all models */
    virtual void computeVjk() = 0;
    /** Compute m_Uil_ array for all models */
    virtual void computeUil() = 0;
    /** compute logRhol during the m-step */
    virtual void mSteplogRhol();
    /** compute logPiek during the m-step */
    virtual void mSteplogPiek();

    /** Parameter finalization common for all models*/
    void commonFinalizeOutput();
    /** compute the vector v_Tk_ and check if the size block is not too small
     *  @return true if the size block is under the threshold, false otherwise
     **/
    bool finalizeStepRows();
    /** compute the vector v_Rl_ and check if the size block is not too small
     *  @return true if the size block is under the threshold, false otherwise
     **/
    bool finalizeStepCols();

    /**Set the known and unknown row labels for semi-supervised coclustering*/
    void setRowLabels(VectorInt const&);
    /**Set the known  and unknown column labels for semi-supervised coclustering*/
    void setColLabels(VectorInt const&);

  private:
    //make assignment operator private
    ICoClustModel& operator=(const ICoClustModel&);
};

inline ModelParameters const& ICoClustModel::modelParameters() const
{ return Mparam_;}

inline void ICoClustModel::setModelParameters(ModelParameters& Mparam)
{ Mparam_ = Mparam;}

inline void ICoClustModel::setEpsilon(STK::Real Epsilon)
{ Mparam_.epsilon_ = Epsilon;}

inline bool ICoClustModel::stopAlgo()
{ return stopAlgo_;}

inline  VectorInt const& ICoClustModel::rowClassificationVector() const
{ return v_Zi_;}

inline  VectorInt const& ICoClustModel::columnClassificationVector() const
{ return v_Wj_;}

inline VectorReal const& ICoClustModel::rowProportions() const
{ return v_Piek_;}

inline VectorReal const& ICoClustModel::colProportions() const
{ return v_Rhol_;}

inline MatrixReal const& ICoClustModel::rowPosteriorProb() const
{return m_Tik_;}

inline MatrixReal const& ICoClustModel::colPosteriorProb() const
{return m_Rjl_;}

inline STK::Real ICoClustModel::likelihood() const
{ return likelihood_;}

inline void ICoClustModel::mSteplogRhol()
{ if(!Mparam_.fixedproportions_) { v_logRhol_=(v_Rl_/Mparam_.nbCol_).log();} }

inline void ICoClustModel::mSteplogPiek()
{ if(!Mparam_.fixedproportions_) { v_logPiek_=(v_Tk_/Mparam_.nbRow_).log();} }

template<class T>
void ICoClustModel::arrangedDataCluster(T& m_ClusterDataij_, const T& m_Dataij_)
{
  STK::VectorXi v_Zi = v_Zi_;
  STK::VectorXi v_Wj = v_Wj_;
  m_ClusterDataij_.resize(Mparam_.nbRow_,Mparam_.nbCol_);
  m_ClusterDataij_ = 0;
  //Rearrange data into clusters
  VectorInt rowincrement(Mparam_.nbrowclust_, 0);
  VectorInt nbindrows(Mparam_.nbrowclust_+1, 0);

  for ( int k = 1; k < Mparam_.nbrowclust_; ++k)
  { nbindrows[k] = (v_Zi==(k-1)).count()+nbindrows[k-1];}

  VectorInt colincrement(Mparam_.nbcolclust_, 0);
  VectorInt nbindcols(Mparam_.nbcolclust_+1, 0);

  for ( int l = 1; l < Mparam_.nbcolclust_; ++l)
  { nbindcols[l]= (v_Wj==(l-1)).count()+nbindcols[l-1];}

  for ( int j = 0; j < Mparam_.nbCol_; ++j)
  {
    m_ClusterDataij_.col(colincrement[v_Wj[j]] + nbindcols[v_Wj[j]]) = m_Dataij_.col(j);
    colincrement[v_Wj[j]]+=1;
  }
  T temp = m_ClusterDataij_;

  for ( int i = 0; i <Mparam_.nbRow_; ++i)
  {
    m_ClusterDataij_.row( rowincrement[v_Zi[i]] + nbindrows[v_Zi[i]]) = temp.row(i);
    rowincrement[v_Zi[i]]+=1;
  }
}

#endif /* ICOCLUSTMODEL_H_ */
