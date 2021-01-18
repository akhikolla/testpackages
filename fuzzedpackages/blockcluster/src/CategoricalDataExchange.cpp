#include "CategoricalDataExchange.h"
#include "coclust/src/Models/CategoricalLBModel.h"

void CategoricalDataExchange::dataOutput(Rcpp::S4& obj,ICoClustModel* model,bool successful)
{
  if(!successful)
  {
    obj.slot("successful") = false;
    obj.slot("message")    = std::string("Co-Clustering Failed! ") + model->errorMsg();
  }
  else
  {
    obj.slot("successful") = true;
    obj.slot("message")    = "Co-Clustering successfully terminated!";
    CategoricalLBModel* ptrLBM;
    MatrixReal dispersion;

    ptrLBM = dynamic_cast<CategoricalLBModel*>(model);
    const std::vector<MatrixReal> mean = ptrLBM->mean();
    std::vector<std::vector<std::vector<double> > > tempmean(mean.size(),std::vector<std::vector<double> >(Mparam_.nbrowclust_,(std::vector<double>(Mparam_.nbcolclust_))));
    for (int h = 0; h < int(mean.size()); ++h)
    {
      for (int k = 0; k < Mparam_.nbrowclust_; ++k)
      {
        for (int l = 0; l < Mparam_.nbcolclust_; ++l)
        { tempmean[h][k][l] = mean[h](k,l);}
      }
    }

    obj.slot("classmean") = Rcpp::wrap(tempmean);
    obj.slot("coclusterdata") =STK::wrap(ptrLBM->arrangedDataClusters());
    obj.slot("rowclass") = STK::wrap(model->rowClassificationVector());
    obj.slot("colclass") = STK::wrap(model->columnClassificationVector());
    obj.slot("rowproportions") = STK::wrap(model->rowProportions());
    obj.slot("columnproportions") = STK::wrap(model->colProportions());
    obj.slot("rowposteriorprob") = STK::wrap(model->rowPosteriorProb());
    obj.slot("colposteriorprob") = STK::wrap(model->colPosteriorProb());
    obj.slot("likelihood") = model->likelihood();
    obj.slot("ICLvalue") = model->iclCriteriaValue();
  }
}

void CategoricalDataExchange::dataInput(Rcpp::S4& obj)
{
  //Rcpp::IntegerMatrix data(SEXP(obj.slot("data")));
  STK::RMatrix<int> data(SEXP(obj.slot("data")));
  m_Dataij_ = data;
  Mparam_.nbRow_ = m_Dataij_.sizeRows();
  Mparam_.nbCol_ = m_Dataij_.sizeCols();

  //Get Strategy
  //Rcpp::S4 strategy(obj.slot("strategy"));
  //get hyper-parameters
  Rcpp::NumericVector hyperparam(obj.slot("hyperparam"));
  a_ = hyperparam(0);
  b_ = hyperparam(1);
}

void CategoricalDataExchange::instantiateModel(ICoClustModel*& model)
{
  if(!strategy_.SemiSupervised)
  {
    switch (strategy_.Model_)
    {
      case pi_rho_multi_:
        Mparam_.fixedproportions_ = true;
        model = new CategoricalLBModel(m_Dataij_,Mparam_,a_,b_);
        break;
      case pik_rhol_multi_:
        Mparam_.fixedproportions_ = false;
        model = new CategoricalLBModel(m_Dataij_,Mparam_,a_,b_);
        break;
      default:
        Rcpp::stop("Wrong Model in CategoricalDataExchange. Please report Bug.");
        break;
    }
  }
  else
  {
    switch (strategy_.Model_)
    {
      case pi_rho_multi_:
        Mparam_.fixedproportions_ = true;
        model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_,a_,b_);
        break;
      case pik_rhol_multi_:
        Mparam_.fixedproportions_ = false;
        model = new CategoricalLBModel(m_Dataij_,v_rowlabels_,v_collabels_,Mparam_,a_,b_);
        break;
      default:
        Rcpp::stop("Wrong Model in CategoricalDataExchange. Please report Bug.");
        break;
    }
  }
}

