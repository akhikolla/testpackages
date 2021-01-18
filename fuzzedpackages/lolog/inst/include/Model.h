#ifndef MODELH_
#define MODELH_

#include "Stat.h"
#include "Offset.h"
#include "StatController.h"
#include <vector>
#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <Rcpp.h>
#include <RcppCommon.h>
#include "ShallowCopyable.h"

namespace lolog{


/*!
 * a representation of an lolog model
 */
template<class Engine>
class Model : public ShallowCopyable{
protected:
    typedef boost::shared_ptr< AbstractStat<Engine > > StatPtr;
    typedef std::vector< StatPtr >  StatVector;

    typedef boost::shared_ptr< AbstractOffset<Engine > > OffsetPtr;
    typedef std::vector< OffsetPtr >  OffsetVector;

    typedef boost::shared_ptr< std::vector<int> > VectorPtr;

    StatVector stats;						//!statistics
    OffsetVector offsets;
    boost::shared_ptr< BinaryNet<Engine> > net;			//!the relevant network

    /**
     * A vector giving the (partial) ordering for vertex inclusion
     */
    VectorPtr vertexOrder;

public:
    Model(){
        //std::cout << "m1";
        boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>());
        net=n;
        vertexOrder = VectorPtr(new std::vector<int>());
    }

    Model(BinaryNet<Engine>& network){
        //std::cout << "m2";
        boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
        net = n;
        vertexOrder = VectorPtr(new std::vector<int>());
    }

    Model(const Model& mod){
        //std::cout << "m3";
        stats = mod.stats;
        offsets = mod.offsets;
        net = mod.net;
        vertexOrder = mod.vertexOrder;
    }

    /*!
     * if deep, then the model statistics are de-aliased
     */
    Model(const Model& mod, bool deep){
        //std::cout << "m4";
        stats = mod.stats;
        offsets = mod.offsets;
        net = mod.net;
        vertexOrder = mod.vertexOrder;
        if(deep){
            for(int i=0;i<stats.size();i++)
                stats[i] = stats[i]->vClone();
            for(int i=0;i<offsets.size();i++)
                offsets[i] = offsets[i]->vClone();
            vertexOrder = boost::shared_ptr< std::vector<int> >(new std::vector<int>);
            *vertexOrder = *mod.vertexOrder;
        }
    }

    virtual ~Model(){}

    /*!
     * R constructor for RCPP
     *
     */
    Model(SEXP sexp){
        //std::cout << "m5";
        boost::shared_ptr<Model> xp = unwrapRobject< Model<Engine> >(sexp);
        stats = xp->stats;
        offsets = xp->offsets;
        net = xp->net;
        vertexOrder = xp->vertexOrder;
    }

    virtual ShallowCopyable* vShallowCopyUnsafe() const{
        return new Model(*this);
    }

    /*!
     * coerce to R object. for RCPP
     */
    operator SEXP() const{
        //std::cout << "mWrap";
        return wrapInReferenceClass(*this,Engine::engineName() + "Model");
    }

    /*!
     * clone the model (but not the network assosiated with it)
     */
    boost::shared_ptr< Model<Engine> > clone() const{
        return boost::shared_ptr< Model<Engine> >(new Model<Engine>(*this, true));
    }

    virtual boost::shared_ptr< Model<Engine> > vClone() const{
        return clone();
    }

    void copy(Model<Engine>& mod){
        stats = mod.stats;
        offsets = mod.offsets;
        net = mod.net;
        vertexOrder = mod.vertexOrder;
    }

    void copy(Model<Engine>& mod,bool deep){
        net = mod.net;
        if(deep){
            stats.resize(mod.stats.size());
            offsets.resize(mod.offsets->vSize());
            for(int i=0;i<stats.size();i++)
                stats[i] = mod.stats[i]->vClone();
            for(int i=0;i<offsets->vSize();i++)
                offsets[i] = mod.offsets[i]->vClone();
            vertexOrder = boost::shared_ptr< std::vector<int> >(new std::vector<int>);
            *vertexOrder = *mod.vertexOrder;
        }else{
            stats = mod.stats;
            offsets = mod.offsets;
            vertexOrder = mod.vertexOrder;
        }
    }

    /*!
     * the model terms
     */
    std::vector<double> terms(){
        int n=0;
        for(int i=0;i<stats.size();i++){
            n += stats.at(i)->vSize();
        }
        for(int i=0;i<offsets.size();i++){
            n += offsets.at(i)->vSize();
        }
        std::vector<double> v(n,0.0);
        int c=0;
        for(int i=0;i<stats.size();i++){
            std::vector<double> vals = stats.at(i)->vValues();
            for(int j=0;j<vals.size();j++){
                v[c] = vals[j];
                c++;
            }
        }
        for(int i=0;i<offsets.size();i++){
            std::vector<double> vals = offsets.at(i)->vValues();
            for(int j=0;j<vals.size();j++){
                v[c] = vals[j];
                c++;
            }
        }
        return v;
    }

    /*!
     * the model parameters
     */
    std::vector<double> thetas(){
        int n=0;
        for(int i=0;i<stats.size();i++){
            n += stats.at(i)->vTheta().size();
        }
        std::vector<double> v(n,0.0);
        int c=0;
        for(int i=0;i<stats.size();i++){
            std::vector<double> vals = stats.at(i)->vTheta();
            for(int j=0;j<vals.size();j++){
                v[c] = vals[j];
                //cout << stats.at(i)->theta()[j];
                c++;
            }
        }
        return v;
    }

    /*!
     * set the model paramters
     */
    void  setThetas(std::vector<double> newThetas){
        int n=0;
        for(int i=0;i<stats.size();i++){
            n += stats.at(i)->vTheta().size();
        }
        if(newThetas.size()!= n){
            //Rcpp::Rcout << n  << " " << newThetas.size() << " ";
            ::Rf_error("Model.setThetas: size mismatch:");
        }
        int c=0;
        for(int i=0;i<stats.size();i++){
            std::vector<double>* vals = &stats.at(i)->vTheta();
            for(int j=0;j<vals->size();j++){
                (*vals)[j] = newThetas[c];
                //cout << stats.at(i)->theta()[j];
                c++;
            }
        }
    }


    /*!
     * the model statistics
     */
    std::vector<double> statistics(){
        int n=0;
        for(int i=0;i<stats.size();i++){
            n += stats.at(i)->vSize();
        }
        std::vector<double> v(n,0.0);
        int c=0;
        for(int i=0;i<stats.size();i++){
            //std::vector<double> vals = stats.at(i)->vStatistics();
            for(int j=0;j<stats.at(i)->vStatistics().size();j++){
                v[c] = stats.at(i)->vStatistics()[j];
                c++;
            }
        }
        return v;
    }
    /*!
     * Copy model statistics into v
     */
    void statistics(std::vector<double>& v){
        int c=0;
        for(int i=0;i<stats.size();i++){
            //std::vector<double> vals = stats.at(i)->vStatistics();
            for(int j=0;j<stats[i]->vStatistics().size();j++){
                v[c] = stats[i]->vStatistics()[j];
                c++;
            }
        }
    }

    /*!
     * returns statistics with names for R
     */
    NumericVector statisticsR(){
        NumericVector res = wrap(statistics());
        res.attr("names") = wrap(names());
        return res;
    }




    /*!
     * returns thetas with names for R
     */
    NumericVector thetasR(){
        NumericVector res = wrap(thetas());
        res.attr("names") = wrap(names());
        return res;
    }

    /*!
     * the model statistic names
     */
    std::vector<std::string> names(){
        int n=0;
        for(int i=0;i<stats.size();i++){
            n += stats.at(i)->vSize();
        }
        std::vector<std::string> v(n,"??");
        int c=0;
        for(int i=0;i<stats.size();i++){
            std::vector<std::string> vals = stats.at(i)->vStatNames();
            for(int j=0;j<vals.size();j++){
                v[c] = vals[j];
                c++;
            }
        }
        return v;
    }

    /*!
     * the model offsets
     */
    std::vector<double> offset(){
        int n=0;
        for(int i=0;i<offsets.size();i++){
            n += offsets.at(i)->vSize();
        }
        std::vector<double> v(n,0.0);
        int c=0;
        for(int i=0;i<offsets.size();i++){
            std::vector<double> vals = offsets.at(i)->vValues();
            for(int j=0;j<vals.size();j++){
                v[c] = vals[j];
                c++;
            }
        }
        return v;
    }

    /*!
     * the log likelihood of the model
     */
    double logLik(){
        double ll = 0.0;
        for(int i=0;i<stats.size();i++){
            ll += stats[i]->vLogLik();
        }
        for(int i=0;i<offsets.size();i++){
            ll += offsets[i]->vLogLik();
        }
        return ll;
    }

    virtual double vLogLik(){
        return logLik();
    }

    /*!
     * add a statistic
     */
    void addStatPtr(StatPtr  s){
        stats.push_back(s);
        s->vCalculate(*net);
    }

    /*!
     * add a statistic to the model
     */
    void addStat(const AbstractStat<Engine>&  s){
        StatPtr ps((&s)->clone());
        ps->vCalculate(*net);
        stats.push_back(ps);
    }

    /*!
     * add a offset
     */
    void addOffsetPtr(OffsetPtr  o){
        offsets.push_back(o);
        o->vCalculate(*net);
    }

    /*!
     * add a offset to the model
     */
    void addOff(const AbstractOffset<Engine>&  o){
        OffsetPtr ps((&o)->vClone());
        ps->vCalculate(*net);
        offsets.push_back(ps);
    }

    /*!
     * add a statistic by name. uses StatController
     */
    void addStatistic(const std::string name, Rcpp::List params){
        AbstractStat<Engine>* ps = StatController<Engine>::getStat(name, params);
        if(ps==NULL){
            ::Rf_error("Invalid stat");
            return;
        }
        ps->vCalculate(*net);
        stats.push_back(StatPtr(ps));
    }

    /*!
     * add a offset by name. uses StatController
     */
    void addOffset(const std::string name, Rcpp::List params){
        AbstractOffset<Engine>* ps = StatController<Engine>::getOffset(name, params);
        if(ps==NULL){
            ::Rf_error("Invalid offset");
            return;
        }
        ps->vCalculate(*net);
        offsets.push_back(OffsetPtr(ps));
    }

    /*!
     * calculates statistics and offsets
     */
    void calculate(){
        calculateStatistics();
        calculateOffsets();
    }

    /*!
     * calculate the statistics
     */
    void calculateStatistics(){
        for(int i=0;i<stats.size();i++){
            stats[i]->vCalculate(*net);
        }
    }

    /*!
     * calculate the statistics
     */
    void calculateOffsets(){
        for(int i=0;i<offsets.size();i++){
            offsets[i]->vCalculate(*net);
        }
    }

    void dyadUpdate(int &from, int &to, std::vector<int> &order, int &actorIndex){
        for(int k=0;k<stats.size();k++){
            stats[k]->vDyadUpdate(*net, from, to, order, actorIndex);
        }
        for(int k=0;k<offsets.size();k++){
            offsets[k]->vDyadUpdate(*net, from, to, order, actorIndex);
        }
    }


    void discreteVertexUpdate(int vertex, int variable, int newValue, std::vector<int> &order, int &actorIndex){
        for(int k=0;k<stats.size();k++)
            stats[k]->vDiscreteVertexUpdate(*net,vertex, variable, newValue);
        for(int k=0;k<offsets.size();k++)
            offsets[k]->vDiscreteVertexUpdate(*net,vertex, variable, newValue);
    }

    void continVertexUpdate(int vertex, int variable, double newValue, std::vector<int> &order, int &actorIndex){
        for(int k=0;k<stats.size();k++)
            stats[k]->vContinVertexUpdate(*net,vertex, variable, newValue, order, actorIndex);
        for(int k=0;k<offsets.size();k++)
            offsets[k]->vContinVertexUpdate(*net,vertex, variable, newValue, order, actorIndex);
    }

    void rollback(){
        for(int k=0;k<stats.size();k++)
            stats[k]->vRollback(*net);
        for(int k=0;k<offsets.size();k++)
            offsets[k]->vRollback(*net);
    }

    /*!
     * get the network
     */
    boost::shared_ptr< BinaryNet<Engine> > network() const{
        return net;
    }

    Rcpp::RObject getNetworkR() const{
        return wrap(*net);
    }

    /*!
     * set the network
     */
    void setNetwork(const BinaryNet<Engine>& network){
        boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
        net = n;
    }
    void setNetworkR(const BinaryNet<Engine>& network){
        boost::shared_ptr< BinaryNet<Engine> > n(new BinaryNet<Engine>(network));
        net = n;
    }
    void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > network){
        net = network;
    }

    const VectorPtr& getVertexOrder() const {
        return vertexOrder;
    }

    void setVertexOrder(const VectorPtr& vertexOrder) {
        this->vertexOrder = vertexOrder;
    }

    std::vector<int> getVertexOrderVector() const {
        if(!vertexOrder)
            return std::vector<int>();
        return *vertexOrder;
    }

    void setVertexOrderVector(std::vector<int> vertexOrder) {
        if(net){
            if(vertexOrder.size() != 0 && net->size() != vertexOrder.size())
                Rf_error("Vertex ordering does not have the same number of elements as there are vertices in the network.");
        }
        if(!this->vertexOrder)
            this->vertexOrder = boost::shared_ptr< std::vector<int> >(new std::vector<int>);
        *this->vertexOrder = vertexOrder;
    }

    bool hasVertexOrder(){
        return this->vertexOrder->size() != 0;
    }

    /*!
     * The independence type of each term
     *
     * \param dyad if true, dyad independence is evaluated, otherwise, order independence is evaluated
     * \param statistic. If true, independence is evaluated for each statistic, otherwise it is evaulated for each offset.
     *
     */
    std::vector<bool> isIndependent(bool dyad, bool statistic){
        bool di;
        if(statistic){
            int n=0;
            for(int i=0;i<stats.size();i++){
                n += stats.at(i)->vSize();
            }
            std::vector<bool> v(n,0.0);
            int c=0;
            for(int i=0;i<stats.size();i++){

                if(dyad)
                    di = stats.at(i)->vIsDyadIndependent();
                else
                    di = stats.at(i)->vIsOrderIndependent();
                //std::vector<double> vals = stats.at(i)->vStatistics();
                for(int j=0;j<stats.at(i)->vStatistics().size();j++){
                    v[c] = di;
                    c++;
                }
            }
            return v;
        }else{
            int n=0;
            for(int i=0;i<offsets.size();i++){
                n += offsets.at(i)->vSize();
            }
            std::vector<bool> v(n,0.0);
            int c=0;
            for(int i=0;i<offsets.size();i++){
                if(dyad)
                    di = offsets.at(i)->vIsDyadIndependent();
                else
                    di = offsets.at(i)->vIsOrderIndependent();
                //std::vector<double> vals = stats.at(i)->vStatistics();
                for(int j=0;j<offsets.at(i)->vSize();j++){
                    v[c] = di;
                    c++;
                }
            }
            return v;
        }
    }

};


#include <Rcpp.h>

}

#endif /* MODELH_ */
