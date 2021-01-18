#ifndef CONSTRAINTH_
#define CONSTRAINTH_
#include <vector>
#include "Offset.h"
#include "limits.h"

namespace lolog {




/*!
 * Constraints are just an offset which is log(0) if satisfied, and
 * log(-Inf) if not.
 */
template<class NetworkEngine, class OffsetEngine>
class Constraint : public Offset<NetworkEngine,OffsetEngine>{
public:
    Constraint() : Offset<NetworkEngine,OffsetEngine>(){}

    Constraint(Rcpp::List params) : Offset<NetworkEngine,OffsetEngine>(params){}

    virtual ~Constraint() {}

    virtual AbstractOffset<NetworkEngine>* vCreateUnsafe(Rcpp::List params) const{
        return createUnsafe(params);
    }

    inline AbstractOffset<NetworkEngine>* createUnsafe(Rcpp::List params) const{
        return new Constraint(params);
    }

    /*!
     * \return an identical un-aliased version of the Stat
     */
    virtual boost::shared_ptr< AbstractOffset<NetworkEngine> > vClone(){
        return clone();
    }

    inline boost::shared_ptr< AbstractOffset<NetworkEngine> > clone(){
        return boost::shared_ptr< AbstractOffset<NetworkEngine> >(cloneUnsafe());
    }

    virtual AbstractOffset<NetworkEngine>* vCloneUnsafe(){
        return cloneUnsafe();
    }
    inline AbstractOffset<NetworkEngine>* cloneUnsafe(){
        return new Constraint<NetworkEngine,OffsetEngine>(*this);
    }

    virtual void vCalculate(const BinaryNet<NetworkEngine>& net){
        calculate(net);
    }

    inline void calculate(const BinaryNet<NetworkEngine>& net){
        this->off.updateLogLik(this->off.initialize(net));
    }

    virtual void vDyadUpdate(const BinaryNet<NetworkEngine>& net,const int& from, int& to,const std::vector<int> &order,const int &actorIndex){
        dyadUpdate(net,from,to, order, actorIndex);
    }

    inline void dyadUpdate(const BinaryNet<NetworkEngine>& net,const int &from,const int &to,const std::vector<int> &order, const int &actorIndex){
        this->off.updateLogLik(this->off.dyadUpdateDistance(net, from, to));
    }

    virtual void vDiscreteVertexUpdate(const BinaryNet<NetworkEngine>& net, const  int& vert,
            const int& variable, const  int& newValue, const  std::vector<int> &order, const  int &actorIndex){
        discreteVertexUpdate(net,vert,variable,newValue,order,actorIndex);
    }

    inline void discreteVertexUpdate(const BinaryNet<NetworkEngine>& net, const  int& vert,
            const int& variable, const  int& newValue, const  std::vector<int> &order, const  int &actorIndex){
        this->off.updateLogLik(this->off.discreteVertexUpdateDistance(net, vert, variable,newValue));
    }

    virtual void vContinVertexUpdate(const BinaryNet<NetworkEngine>& net,  const int& vert,
            const int& variable, const  double& newValue, const  std::vector<int> &order, const  int &actorIndex){
        continVertexUpdate(net,vert,variable,newValue,order,actorIndex);
    }

    inline void continVertexUpdate(const BinaryNet<NetworkEngine>& net,  const int& vert,
            const int& variable, const  double& newValue, const  std::vector<int> &order, const  int &actorIndex){
        this->off.updateLogLik(this->off.continVertexUpdateDistance(net, vert, variable, newValue));
    }

    virtual void vRollback(const BinaryNet<NetworkEngine>& net){
        rollback(net);
    }

    inline void rollback(const BinaryNet<NetworkEngine>& net){
        this->off.rollback(net);
    }

};


/*!
 * Constraints restrict the sample space of the model by setting it's log likelihood to
 * -infinity when outside the restricted space. If the network does not satisfy the
 * constraint, it will be snapped to it when MCMC is run.
 */
template<class Engine>
class BaseConstraint{
private:
    double logValue;
    double lastValue;
protected:

public:
    virtual ~BaseConstraint(){};

    /*!
     * calculate how many steps away the constraint is from being satisfied
     */
    double initialize(const BinaryNet<Engine>& net){
        Rf_error("initialize must be implemented");
    }

    /*!
     * dyad update for how many steps away the constraint is from being satisfied
     */
    double dyadUpdateDistance(const BinaryNet<Engine>& net,const int& from, const int& to){
        Rf_error("distanceFromSatisfaction (dyad) must be implemented");
    }

    /*!
     * disc vertex update for how many steps away the constraint is from being satisfied
     */
    double discreteVertexUpdateDistance(const BinaryNet<Engine>& net,
            const int& vert, const int& variable, const int& newValue){
        Rf_error("distanceFromSatisfaction (vertex) must be implemented");
    }

    /*!
     * cont vertex update for how many steps away the constraint is from being satisfied
     */
    double continVertexUpdateDistance(const BinaryNet<Engine>& net, const int& vert,
            const int& variable, const double& newValue){
        Rf_error("distanceFromSatisfaction (vertex cont) must be implemented");
    }


    /*!
     * takes the distance from constraint satisfaction, and sets the value of the
     * constraint offset o(T).
     *
     */
    void updateLogLik(double dist){
        lastValue = logValue;
        if(near(dist,0.0)){
            logValue = 0.0;
        }else{
            logValue = -10000000000.0 - 100000.0*dist;
        }
    }

    void rollback(const BinaryNet<Engine>& net){
        logValue = lastValue;
    }

    /*!
     * value of offset
     */
    double logLik(){
        return logValue;
    };

    std::vector<double> values(){
        return std::vector<double>(1,logValue);
    }

    /*!
     * number of terms
     */
    int size(){
        return 1;
    }

    void calculate(const BinaryNet<Engine>& net){
        Rf_error("BaseConstraint calculate should not be called");
    }
    void dyadUpdate(const BinaryNet<Engine>& net, int from, int to, const  std::vector<int> &order, const  int &actorIndex){
        Rf_error("BaseConstraint dyadUpdate should not be called");
    }
    void discreteVertexUpdate(const BinaryNet<Engine>& net, int vert,
            int variable, int newValue, const  std::vector<int> &order, const  int &actorIndex){
        Rf_error("BaseConstraint discreteVertexUpdate should not be called");
    }
    void continVertexUpdate(const BinaryNet<Engine>& net, int vert,
            int variable, double newValue, const  std::vector<int> &order, const  int &actorIndex){
        Rf_error("BaseConstraint continVertexUpdate should not be called");
    }

    bool isOrderIndependent(){
        return false;
    }

    bool isDyadIndependent(){
        return false;
    }
};


} /* namespace lolog */
#endif /* CONSTRAINTH_ */
