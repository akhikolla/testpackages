#ifndef RANKCLUSTER_H_
#define RANKCLUSTER_H_

/**@file RankCluster.h
 * @brief Definition of the class @c RankCluster and struct PartialRank, SEMparameters and OutParameters
 */
//#include <RcppEigen.h>
#include <vector>
#include <set>
#include <utility>
#include "Eigen/Dense"

#include "functions.h"


// create your own data structures
struct PartialRank
{
    ///rank
    std::vector<int> rank;
    ///order of presentation
    std::vector<int> y;
    ///if true, the rank contains partial or ties
    bool isNotFull;
    /// missing element of the rank or ties
    std::vector<std::vector<int> > missingData;
    ///index of the 0 or ties
    std::vector<std::vector<int> > missingIndex;

};

struct SEMparameters
{
  ///number of iteration in the gibbs of SE step for each dimension of rank
  std::vector<int> nGibbsSE;
  ///number of iteration in the gibbs of M step for each dimension of rank
  std::vector<int> nGibbsM;
  ///maximum number of iteration of the SEM algorithm
  int maxIt;
  ///burn-in period of SEM algorithm
  int burnAlgo;
  /// number of iteration in the gibbs of the likelihood computation
  int nGibbsL;
  ///burn-in period of the likelihood computation
  int burnL;
  ///maximum number of try of the SEM
  int maxTry;
  ///if true print details
  bool detail;
};

struct OutParameters
{
  ///loglikelihood
  double L;
  ///bic criterion
  double bic;
  ///icl criterion
  double icl;
  ///
  Eigen::ArrayXXd tik;
  ///
  Eigen::ArrayXd entropy;
  ///
  Eigen::ArrayXXd probabilities;
  ///percentage of confidence in final estimation of missing data
  std::vector<std::vector<std::vector<double> > > partialRankScore;

  //algorithm initialization
  std::vector<std::vector<std::vector<int> > > initialPartialRank;
  std::vector<std::vector<double> > initialP;
  std::vector<int> initialZ;
  std::vector<double> initialProportion;
  std::vector<std::vector<std::vector<int> > > initialMu;

  //distance between parameters
  std::vector<std::vector<double> > distProp;
  std::vector<std::vector<std::vector<double> > > distP;
  std::vector<std::vector<std::vector<int> > > distMu;
  std::vector<double> distZ;
  std::vector<std::vector<std::vector<int> > > distPartialRank;

};

class RankCluster
{
  public:
        ///defaut constructor
    RankCluster();
    /**
    * @param X data one row= a multi dimensionnal rank
    * @param g number of clusters
    * @param m size of rank of each dimension
    * @param param parameters of SEM algorithm
    */
    RankCluster(std::vector<std::vector<int> > const& X,int g, std::vector<int> const& m, SEMparameters const& param);
    //constructor with given initialization of parameters
    RankCluster(std::vector<std::vector<int> > const& X, std::vector<int> const& m, SEMparameters const& param,
                    std::vector<double> const& proportion, std::vector<std::vector<double> > const& p,
                    std::vector<std::vector<std::vector<int> > > const& mu);

    ///copy constructor
    RankCluster(RankCluster& rankClusterObject);

    ///destructor
    virtual ~RankCluster();

    /// run the SEM algorithm
    void run();

    //getters
    inline int d() const {return d_;}
    inline int n() const {return n_;}
    inline int g() const {return g_;}
    inline std::vector<int> m() const {return m_;}
    inline std::vector<int> z() const {return z_;}
    inline std::vector<std::vector<double> > p() const {return p_;}
    inline std::vector<std::vector<std::vector<int> > >  mu() const {return mu_;}
    inline std::vector<double> proportion()  const  {return proportion_;}
    inline std::vector<std::vector<int> > indexPartialData() const {return indexPartialData_;}
    inline std::vector<std::vector<PartialRank> > data() const {return data_;}
    inline std::vector<int> rank(int dim, int index) const {return data_[dim][index].rank;}
    inline bool dataOk() const {return dataOk_;}
    inline bool convergence() const {return convergence_;}
    inline bool partial() const {return partial_;}
    inline std::vector<std::vector<int> > indexPb() const {return indexPb_;}
    inline SEMparameters parameter() const {return parameter_;}
    
    //output getters
    inline Eigen::ArrayXXd tik() const {return output_.tik;}
    inline Eigen::ArrayXd entropy() const {return output_.entropy;}
    inline Eigen::ArrayXXd probabilities() const {return output_.probabilities;}
    Eigen::ArrayXd probability() const;
    inline double bic() const {return output_.bic;}
    inline double icl() const {return output_.icl;}
    inline double L() const {return output_.L;}
    inline std::vector<std::vector<std::vector<int> > > initialPartialRank() const {return output_.initialPartialRank;}
    inline std::vector<std::vector<double> > initialP() const {return output_.initialP;}
    inline std::vector<int> initialZ() const {return output_.initialZ;}
    inline std::vector<std::vector<std::vector<int> > > initialMu() const {return output_.initialMu;}
    inline std::vector<double> initialProportion() const {return output_.initialProportion;}
    inline std::vector<std::vector<double> > distProp() const {return output_.distProp;}
    inline std::vector<std::vector<std::vector<double> > > distP() const {return output_.distP;}
    inline std::vector<std::vector<std::vector<int> > > distMu() const {return output_.distMu;}
    inline std::vector<double> distZ() const {return output_.distZ;}
    inline std::vector<std::vector<std::vector<int> > > distPartialRank() const {return output_.distPartialRank;}
    inline std::vector<std::vector<std::vector<double> > > partialRankScore() const {return output_.partialRankScore;}
  
    ///reestimation of criterion
    void estimateCriterion(double &L,double &bic,double &icl);

  protected:
        /**convert X in vector<vector<PartialRank>>
        * @param X raw data one row= a multi dimensionnal rank
        */
    void conversion2data(std::vector<std::vector<int> > const& X);
    /** read rank. used in conversion2data
        * @param X raw data one row= a multi dimensionnal rank
        * @param dim actual dimension
        * @param j actual index of the sample
        * @param indM transformation of m_
        */
    void readRankingRank(std::vector<std::vector<int> > const& X, int const& dim, int const& j, std::vector<int> const& indM);

    ///initialization of parameters
    void initialization();
    /** SE step */
    void SEstep();
    /**unidimensionnal gibbs sampler for y estimation
    * @param indexDim index of the dimension (<d_)
    */
    void gibbsY(int indexDim);
    ///simulation of z_
    void zSimulation();
    /**unidimensionnal gibbs sampler for partial rank estimation
    * @param indexDim index of the dimension (<d_)
    */
    void gibbsX(int indexDim);
    ///M step
    void Mstep();
    /** simulation of mu
    * @param indexDim index of the dimension
    * @param indCl index of teh cluster
    */
    void simuM(int indexDim,int indCl);
    /** likelihood step
     * @param listeMu mu for every iterations of the SEM
     * @param resP pi for every iterations of the SEM
     * @param resProp proportion for every iterations of the SEM
     */
    void likelihood(std::vector<std::vector<std::vector<std::vector<int> > > > &listeMu,std::vector<std::vector<std::vector<double> > > &resP,
            std::vector<std::vector<double> > &resProp);

    /** compute the log likelihood for a set of parameter
     * @param mu estimated central rank
     * @param p estimated p (dispersion parameter)
     * @param proportion estimated proportion
     * @param tik used for store tik
     * @param Y used for store estimated y
     * @param xTemp used for store estimated partial rank
     * @param score used for confidence in estimated partial rank
     * @param iterproba probability of each individuals at each iteration
     */
    double computeLikelihood(std::vector<std::vector<std::vector<int> > > const& mu, std::vector<std::vector<double> > const& p,
        std::vector<double> const& proportion, Eigen::ArrayXXd &tik, std::vector<std::vector<std::vector<int> > > &Y,
        std::vector<std::vector<std::vector<int> > > &xTemp, Eigen::ArrayXXd &probabilities,
        std::vector<std::vector<std::vector<double> > > &score);
    ///compute the final z_
    void computePartition();
    ///compute distance between final parameters and each iteration parameters
    void computeDistance(std::vector<std::vector<double> > const& resProp,std::vector<std::vector<std::vector<double> > > const& resP,
        std::vector<std::vector<std::vector<std::vector<int> > > > const& resMu,std::vector<std::vector<int> > const& resZ,
        std::vector<std::vector<std::vector<std::vector<int> > > > const& resDonneesPartiel);


  private:
    ///contains the size of rank for each dim
    std::vector<int> m_;
    ///number of individuals
    int n_;
    ///number of dimension
    int d_;
    ///number of cluster
    int g_;
    ///data of teh form data[dimension[index]]
    std::vector<std::vector<PartialRank> > data_;
    ///estimated cluster of each individual
    std::vector<int> z_;
    /// estimated rank parameter of each cluster :  mu_[dimension][cluster][indice]
    std::vector<std::vector<std::vector<int> > > mu_;
    /// estimated probability parameter of each cluster :  p_[dimension][cluster]
    std::vector<std::vector<double> > p_;
    /// estimated proportion of the mixture model
    std::vector<double> proportion_;
    ///algorithm parameters
    SEMparameters parameter_;
    ///distance and initialization of the algorithm
    OutParameters output_;
    ///true if there is partial rank in the data
    bool partial_;
    ///index of partial data
    std::vector<std::vector<int> > indexPartialData_;
    /// if true, SEM has converged
    bool convergence_;
    /// if true, good data
    bool dataOk_;
    ///index of rank with problem for each dimension
    std::vector<std::vector<int> > indexPb_;
};

#endif /* RANKCLUSTER_H_ */
