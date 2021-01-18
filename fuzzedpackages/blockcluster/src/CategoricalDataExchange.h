#ifndef CATEGORICALDATAEXCHANGE_H_
#define CATEGORICALDATAEXCHANGE_H_
/**@file CategoricalDataExchange.h
 * @brief
 */
#include "IDataExchange.h"
#include "coclust/src/Models/CategoricalLBModel.h"
class CategoricalDataExchange: public IDataExchange
{
  public:
    CategoricalDataExchange():a_(1), b_(1) {};
    virtual void dataOutput(Rcpp::S4& obj,ICoClustModel*,bool);
    virtual void dataInput(Rcpp::S4 & obj);
    virtual void instantiateModel(ICoClustModel*& model);
    inline const MatrixInt& GetData() const {return m_Dataij_;}
    virtual ~CategoricalDataExchange(){};
  protected:
    MatrixInt m_Dataij_;
    STK::Real a_,b_;
};

#endif /* CATEGORICALDATAEXCHANGE_H_ */
