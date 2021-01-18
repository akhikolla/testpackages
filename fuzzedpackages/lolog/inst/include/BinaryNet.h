



#ifndef NETH_
#define NETH_

#include <vector>
#include <map>
#include <utility>
#include <Rcpp.h>
#include "DirectedVertex.h"
#include "UndirectedVertex.h"
#include "VarAttrib.h"
#include "util.h"
#include "ShallowCopyable.h"
#include <memory>
#include <boost/shared_ptr.hpp>
#include <stdint.h>

namespace lolog {

using namespace Rcpp;


#ifdef _WIN32
typedef unsigned __int64	unsigned64_t;
#else

typedef uint64_t			unsigned64_t;

#endif

/*!
 * The network. the fundamental structure of this package. Takes Engine as
 * a template parameter, which controls the underlying representation of the network.
 * currently Directed and Undirected are implemented.
 *
 * To implement a different type of network
 * or network data structure, create the appropriate Engine and optionally subclass BinaryNet
 * with new features. This compile time inheritance structure maximizes modualrity and
 * performance.
 *
 */
template <class Engine>
class BinaryNet : public ShallowCopyable{
protected:
    Engine engine; /*!< the driver of the interface */

public:

    typedef typename Engine::NeighborIterator NeighborIterator;
    /*!
     * constructor
     */
    BinaryNet() : engine(){}

    /*!
     * shallow copy constructor
     * \param net net to copy
     */
    BinaryNet(const BinaryNet& net) : engine(net.engine){}

    /*!
     * copy constructor
     * \param net network to copy
     * \param deep should a full clone be made
     */
    BinaryNet(const BinaryNet& net,bool deep) : engine(net.engine,deep){}

    /*!
     * constructor from R S4 object or external pointer. For RCPP
     */
    BinaryNet(SEXP sexp){
        boost::shared_ptr<BinaryNet> ptr = unwrapRobject< BinaryNet<Engine> >(sexp);
        engine = Engine(ptr->engine);
    }

    /*!
     * destructor
     */
    ~BinaryNet(){}

    /*!
     * coerce to R Reference class. For RCPP
     */
    operator SEXP() const{
        std::string rClassName = Engine::engineName() + "Net";
        return wrapInReferenceClass(*this,rClassName);
    }

    static std::string engineName(){
        return Engine::engineName();
    }

    /*!
     * shallow copy
     */
    BinaryNet& operator=(const BinaryNet& rhs){
        if(this != &rhs){
            engine = rhs.engine;
        }
        return *this;
    }

    /*!
     * deep copy
     *  \returns a deep copy of the network
     */
    boost::shared_ptr<BinaryNet<Engine> >  clone() const{
        return boost::shared_ptr<BinaryNet<Engine> >(new BinaryNet(*this,true));
    }

    virtual ShallowCopyable* vShallowCopyUnsafe() const{
        return new BinaryNet(*this, false);
    }

    /*!
     * add a vertex to the network
     */
    void addVertex(){
        engine.addVertex();
    }

    /*!
     * remove a vertex
     */
    void removeVertex(int which){
        engine.removeVertex();
    }

    /*!
     * permute the order of the nodes
     */
    void reorderVertices(std::vector<int> order){
        engine.reorderVertices();
    }

    /*!
     * test if an edge exists
     *
     * \param from the id of the from node
     * \param to the id of the to node
     * \returns whether the edge exists
     */
    bool hasEdge(int from,int to) const{
        return engine.hasEdge(from,to);
    }

    /*!
     * set the value of a dyad
     *
     * \param from the id of the from node
     * \param to the id of the to node
     * \param value whether the dyad is to be set to an edge or not
     */
    void setDyad(int from,int to,bool value){
        if(value)
            engine.addEdge(from,to);
        else
            engine.removeEdge(from,to);
    }

    /*!
     * adds an edge to the network
     * \param from the id of the from node
     * \param to the id of the to node
     */
    void addEdge(int from,int to){
        engine.addEdge(from,to);
    }



    /*!
     * removes an edge
     * \param from the id of the from node
     * \param to the id of the to node
     * \returns whether the edge was present
     */
    bool removeEdge(int from,int to){
        return engine.removeEdge(from,to);
    }

    /*!
     * removes all edges
     */
    void emptyGraph(){
        engine.emptyGraph();
    }

    /*!
     * adds an edge to a dyad if none exists, otherwise removes the edge.
     * \param from the id of the from node
     * \param to the id of the to node
     */
    void toggle(int from,int to){
        bool hadEdge = removeEdge(from,to);
        if(!hadEdge)
            addEdge(from,to);
    }

    /*!
     * network size
     * \returns the number of vertices in the network
     */
    int size() const{
        return engine.size();
    }

    /*!
     * is a dyad missing
     * \param from the id of the from node
     * \param to the id of the to node
     * \returns true if dyad is unobserved
     */
    bool isMissing(int from, int to) const{
        return engine.isMissing(from,to);
    }

    /*!
     * The number of missing (out-)dyads connected to node from.
     * \param from id of node
     */
    int nMissing(int from) const{
        return engine.nMissing(from);
    }

    /*!
     * set a  dyad's missingness
     * \param from the id of the from node
     * \param to the id of the to node
     * \param missing if true dyad is set to missing, otherwise it is set to observed.
     * \returns true if dyad was unobserved
     */
    bool setMissing(int from, int to, bool missing){
        return engine.setMissing(from,to,missing);
    }

    /*!
     * sets all dyads to missing
     *
     */
    void setAllDyadsMissing(){
        engine.setAllDyadsMissing();
    }

    /*!
     * sets all (out-)dyads to missing
     *
     * \param nodes the node ids
     * \param missing should be set to missing or observed
     *
     */
    void setAllDyadsMissing(std::vector<int> nodes,bool missing){
        engine.setAllDyadsMissing(nodes,missing);
    }

    /*!
     * sets all dyads to observed
     *
     */
    void setAllDyadsObserved(){
        engine.setAllDyadsObserved();
    }

    /*!
     * missingness mask
     * \returns a smart pointer to a vector of dyads that are unobserved
     */
    boost::shared_ptr< std::vector< std::pair<int,int> > > missingDyads() const{
        return engine.missingDyads();
    }

    /*!
     * number of edges
     * \returns the number of edges in the network
     */
    int nEdges() const{
        return engine.nEdges();
    }

    /*!
     * are edges directed.
     */
    bool isDirected() const{
        return engine.isDirected();
    }
    /*!
     * \returns the maximum number of edges possible
     */
    unsigned64_t maxEdges() const{
        return engine.maxEdges();
    }

    /*!
     * nodal in-degree
     * \param which the id of the node
     * \returns the number of in-edges
     */
    int indegree(int which) const{
        return engine.indegree(which);
    }

    /*!
     * in-neighbors
     * \param which the id of the node
     * \returns a collection containing the ids of the
     * 			nodes that have an edge pointing to the node
     */
    template<class Collection>
    Collection inneighbors(int which) const{
        return engine.template inneighbors<Collection>(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator inBegin(int which) const{
        return engine.inBegin(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator inEnd(int which) const{
        return engine.inEnd(which);
    }

    /*
	const Set& inneighbors(int which) const{
		return engine.inneighbors(which);
	}
     */
    /*!
     * nodal out-degree
     * \param which the id of the node
     * \returns the number of out-edges
     */
    int outdegree(int which) const{
        return engine.outdegree(which);
    }

    /*!
     * out-neighbors
     * \param which the id of the node
     * \returns a collection containing the ids of the
     * 			nodes that have an edge pointing from the node
     */
    template<class Collection>
    Collection outneighbors(int which) const{
        return engine.template outneighbors< Collection >(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator outBegin(int which) const{
        return engine.outBegin(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator outEnd(int which) const{
        return engine.outEnd(which);
    }
    /*
	const Set& outneighbors(int which) const{
		return engine.outneighbors(which);
	}
     */
    /*!
     * nodal degree
     * throws error for directed nets
     * \param which the id of the node
     * \returns the number of edges
     */
    int degree(int which) const{
        return engine.degree(which);
    }

    /*!
     * neighbors
     * throws error for directed nets
     * \param which the id of the node
     * \returns a collection containing the ids of the
     * 			nodes that share an edge  the node
     */
    template<class Collection>
    Collection neighbors(int which) const{
        return engine.template neighbors<Collection>(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator begin(int which) const{
        return engine.begin(which);
    }

    /*!
     * Access to a bidirectional iterator traversing the neighbors
     * in increasing order
     * \param which id of node
     */
    NeighborIterator end(int which) const{
        return engine.end(which);
    }
    /*
	const Set& neighbors(int which) const{
		return engine.neighbors(which);
	}
     */
    /*!
     * select a random dyad
     * \returns a pair of ids representing the dyad
     */
    std::pair<int,int> randomDyad() const{
        return engine.randomDyad();
    }

    /*!
     * select a random dyad
     * \param toggle a pair of ids representing the dyad
     */
    void randomDyad(std::pair<int,int>& toggle){
        engine.randomDyad(toggle);
    }

    /*!
     * select a random edge
     * \returns a pair of ids representing the edge
     */
    std::pair<int,int> randomEdge() const{
        return engine.randomEdge();
    }

    /*!
     * choose a random dyad from a node.
     * \param from the node. Only outedges from this node are considered.
     * \param missing if true, only missing dyads are considered
     */
    int randomDyad(int from,bool missing){
        return engine.randomDyad(from,missing);
    }

    /*!
     * select a random non-edge
     * \returns a pair of ids representing the dyad
     */
    std::pair<int,int> randomNonEdge() const{
        return engine.randomNonEdge();
    }

    /*!
     * get the edge list of the network
     */
    boost::shared_ptr< std::vector<std::pair<int,int> > > edgelist() const{
        return engine.edgelist();
    }



    /*!
     * the names of the continuous variables
     */
    std::vector<std::string> continVarNames() const{
        return engine.continVarNames();
    }

    /*!
     * get a value of a continuous variable at a vertex
     * \param which the variable index
     * \param at the id of the vertex
     * \returns the value
     */
    double continVariableValue(int which,int at) const{
        return engine.continVariableValue(which,at);
    }

    /*!
     * set a value of a continuous variable at a vertex
     * \param which the variable index
     * \param at the id of the vertex
     * \param newValue the new value
     */
    void setContinVariableValue(int which,int at,double newValue){
        engine.setContinVariableValue(which,at,newValue);
    }

    /*!
     * the variable attributes for a continuous variable
     * \param which the index of the variable
     */
    ContinAttrib continVariableAttributes(int which){
        return engine.continVariableAttributes(which);
    }

    /*!
     * removes a continuous variable
     * \param which the index of the variable to remove
     */
    void removeContinVariable(int which){
        engine.removeContinVariable(which);
    }

    /*!
     * add a continuous variable
     * \param vals a vector of values
     * \param attribs the variable attributes
     */
    void addContinVariable(const std::vector<double>& vals,ContinAttrib& attribs){
        engine.addContinVariable(vals,attribs);
    }

    /*!
     * is variable observed
     * \param which vertex id
     * \param at variable id
     */
    bool continVariableObserved(int which,int at){
        return engine.continVariableObserved(which,at);
    }

    /*!
     * is variable observed
     * \param which variable id
     */
    std::vector<bool> continVariableObserved(int which){
        return engine.continVariableObserved(which);
    }

    /*!
     * sets missing mask
     * \param which variable id
     * \param at vertex id
     * \param observed true if observed
     */
    void setContinVariableObserved(int which,int at,bool observed){
        engine.setContinVariableObserved(which,at,observed);
    }

    /*!
     * the names of the discrete variables
     */
    std::vector<std::string> discreteVarNames() const{
        return engine.discreteVarNames();
    }

    /*!
     * get the values for a discrete variable
     * \param which the index of the variable
     * \returns its values
     */
    std::vector<int> discreteVariableValues(int which) const{
        return engine.discreteVariableValues(which);
    }

    /*!
     * get a value of a discrete variable at a vertex
     * \param which the variable index
     * \param at the id of the vertex
     * \returns the value
     */
    int discreteVariableValue(int which,int at) const{
        return engine.discreteVariableValue(which,at);
    }

    /*!
     * set a value of a discrete variable at a vertex
     * \param which the variable index
     * \param at the id of the vertex
     * \param newValue the new value
     */
    void setDiscreteVariableValue(int which,int at,int newValue){
        engine.setDiscreteVariableValue(which,at,newValue);
    }

    /*!
     * get the label of a discrete variable at a vertex
     * \param which the variable index
     * \param at the id of the vertex
     * \returns the label
     */
    std::string discreteVariableLabel(int which,int at) const{
        return engine.discreteVariableLabel(which,at);
    }

    /*!
     * get the labels of a discrete variable
     * \param which the index of the variable
     * \returns the variable
     */
    std::vector<std::string> discreteVariable(int which) const{
        return engine.discreteVariable(which);
    }

    /*!
     * the variable attributes for a discrete variable
     * \param which the index of the variable
     */
    DiscreteAttrib discreteVariableAttributes(int which) const{
        return engine.discreteVariableAttributes(which);
    }

    /*!
     * removes a discrete variable
     * \param which the index of the variable to remove
     */
    void removeDiscreteVariable(int which){
        engine.removeDiscreteVariable(which);
    }

    /*!
     * is variable observed
     * \param which variable id
     * \param at vertex id
     */
    bool discreteVariableObserved(int which,int at){
        return engine.discreteVariableObserved(which,at);
    }

    /*!
     * is variable observed
     * \param which variable id
     */
    std::vector<bool> discreteVariableObserved(int which){
        return engine.discreteVariableObserved(which);
    }

    /*!
     * sets missing mask
     * \param which variable id
     * \param at vertex id
     * \param observed true if observed
     */
    void setDiscreteVariableObserved(int which,int at,bool observed){
        engine.setDiscreteVariableObserved(which,at,observed);
    }


    /*!
     * add a discrete variable
     * \param vals a vector of values
     * \param attribs the variable attributes
     */
    void addDiscreteVariable(const std::vector<int>& vals,DiscreteAttrib& attribs){
        engine.addDiscreteVariable(vals,attribs);
    }

    /*!
     * adds a discrete variable from an R object
     * \param robj a vector coersable into a factor with the same length as the number of nodes
     * \param name the name of the variable
     */
    void addDiscreteVariableR(RObject robj,std::string name){
        engine.addDiscreteVariableR(robj,name);
    }



    /***************************************************************************/
    //
    //								R Interface
    //
    /***************************************************************************/
    //
    //	Note: 	indices coming from and going to R start at 1,
    //			where as the internal representation starts at 0.
    //
    //			Currently this is not as modular as it should be. these should all
    //			be pass though functions to Engine
    /*!
     * construct from edgelist.
     * \param edgeList the edgelist
     * \param numNodes the number of nodes
     */
    BinaryNet(Rcpp::IntegerMatrix edgeList,int numNodes) : engine(edgeList,numNodes){}

    /*!
     * deep copy
     *  \returns an R Reference Class deep copy of the network
     */
    Rcpp::RObject cloneR() const{
        return wrap(BinaryNet(*this,true));
    }

    /*!
     * sets the specified dyads to the supplied values
     *
     * \param from ids from
     * \param to ids to
     * \param value logical, may be NA
     */
    void setDyadsR(IntegerVector from,IntegerVector to, LogicalVector values){
        if(from.size()!=to.size() || from.size() != values.size())
            ::Rf_error("setDyadsR: vectors must be of the same length");
        if(!is_true(all(from>0)) || !is_true(all(from <= size())))
            ::Rf_error("setDyadsR: range check");
        if(!is_true(all(to>0)) || !is_true(all(to <= size())))
            ::Rf_error("setDyadsR: range check");
        IntegerVector::iterator fit = from.begin();
        IntegerVector::iterator fend = from.end();
        IntegerVector::iterator tit = to.begin();
        IntegerVector::iterator vit = values.begin();
        while(fit!=fend){
            setMissing(*fit-1,*tit-1,(*vit) == NA_LOGICAL);
            if((*vit) != NA_LOGICAL)
                setDyad(*fit-1,*tit-1,(*vit)==true);
            fit++;
            tit++;
            vit++;
        }
    }

    /*!
     * gets dyad values
     *
     * \param from ids from
     * \param to ids to
     */
    LogicalVector getDyadsR(IntegerVector from,IntegerVector to){
        if(from.size()!=to.size())
            ::Rf_error("setDyadsR: vectors must be of the same length");
        if(!is_true(all(from>0)) || !is_true(all(from <= size())))
            ::Rf_error("setDyadsR: range check");
        if(!is_true(all(to>0)) || !is_true(all(to <= size())))
            ::Rf_error("setDyadsR: range check");
        IntegerVector::iterator fit = from.begin();
        IntegerVector::iterator fend = from.end();
        IntegerVector::iterator tit = to.begin();

        LogicalVector values(from.size());
        IntegerVector::iterator vit = values.begin();
        while(fit!=fend){
            *vit = hasEdge(*fit-1,*tit-1);
            if(isMissing(*fit-1,*tit-1))
                *vit = NA_LOGICAL;
            fit++;
            tit++;
            vit++;
        }
        return values;
    }


    /*!
     *  get the edge list of the network as an IntegerMatrix
     */
    Rcpp::IntegerMatrix edgelistR1() const{
        return engine.edgelistR(false);
    }

    /*!
     * get the edge list of the network as an IntegerMatrix
     */
    Rcpp::IntegerMatrix edgelistR2(bool includeMissing) const{
        return engine.edgelistR(includeMissing);
    }

    /*!
     * for the [ operator
     */
    Rcpp::LogicalMatrix getDyadMatrixR(Rcpp::IntegerVector from,Rcpp::IntegerVector to,
            bool maskMissing){
        if(!is_true(all(from>0)) || !is_true(all(from <= size())))
            ::Rf_error("getDyadMatrixR: range check");
        if(!is_true(all(to>0)) || !is_true(all(to <= size())))
            ::Rf_error("getDyadMatrixR: range check");
        LogicalMatrix result(from.size(),to.size());
        for(int i=0;i<from.size();i++){
            for(int j=0;j<to.size();j++){
                if(maskMissing && isMissing(from[i]-1,to[j]-1))
                    result(i,j) = NA_LOGICAL;
                else
                    result(i,j) = hasEdge(from[i]-1,to[j]-1);
            }
        }
        return result;
    }

    /*!
     * for the [<- operator
     */
    void setDyadMatrixR(Rcpp::IntegerVector from,Rcpp::IntegerVector to,
            Rcpp::LogicalMatrix mat){
        if(!is_true(all(from>0)) || !is_true(all(from <= size())))
            ::Rf_error("setDyadMatrixR: range check");
        if(!is_true(all(to>0)) || !is_true(all(to <= size())))
            ::Rf_error("setDyadMatrixR: range check");
        if(mat.nrow()!=from.size())
            ::Rf_error("setDyadMatrixR: number of rows in matrix does not match assignment");
        if(mat.ncol()!=to.size())
            ::Rf_error("setDyadMatrixR: number of column in matrix does not match assignment");
        for(int i=0;i<from.size();i++){
            for(int j=0;j<to.size();j++){
                int value = mat(i,j);
                if(value==NA_LOGICAL)
                    setMissing(from[i]-1,to[j]-1,true);
                else
                    setDyad(from[i]-1,to[j]-1,value);
            }
        }
    }

    /*!
     * A list of character vectors representing the discrete and continuous variables
     */
    Rcpp::RObject getVariableNamesR1(bool simplify){
        Rcpp::List result;
        std::vector<std::string> tmp;
        tmp.push_back("discrete");
        tmp.push_back("continuous");
        result.push_back(wrap(discreteVarNames()));
        result.push_back(wrap(continVarNames()));
        result.attr("names") = wrap(tmp);
        if(!simplify){
            return result;
        }else{
            Language call("unlist",result);
            return call.eval();
        }
    }

    Rcpp::RObject getVariableNamesR2(){
        return getVariableNamesR1(false);
    }

    /*!
     * get a variable by name
     */
    Rcpp::RObject getVariableR(std::string name){
        return engine.getVariableR(name,true);
    }

    /*!
     * get a variable by name
     */
    Rcpp::RObject getVariableR1(std::string name,bool maskMissing){
        return engine.getVariableR(name,maskMissing);
    }

    /*!
     * set a variable
     */
    void setVariableR(std::string name, SEXP variable){
        engine.setVariableR(variable, name);
    }

    /*!
     * the number of missing (out-)dyads around from
     */
    IntegerVector nMissingR(IntegerVector from){
        if(!is_true(all(from>0)) || !is_true(all(from <= size())))
            ::Rf_error("nMissingR: range check");
        IntegerVector vec(from.size());
        IntegerVector::iterator toIt = vec.begin();
        IntegerVector::iterator it = from.begin();
        IntegerVector::iterator end = from.end();
        while(it!=end){
            *toIt = nMissing((*it)-1);
            it++;
            toIt++;
        }
        return vec;
    }

    int nEdgesR1(bool includeMissing){
        if(!includeMissing)
            return nEdges();
        else{
            return edgelistR1().nrow();
        }
    }

    int nEdgesR2(){
        return nEdgesR1(false);
    }

    /*!
     * in degree
     */
    IntegerVector indegreeR(IntegerVector node){
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("indegreeR: range check");
        IntegerVector vec(node.size());
        IntegerVector::iterator vecIt = vec.begin();
        IntegerVector::iterator it = node.begin();
        IntegerVector::iterator end = node.end();
        while(it!=end){
            NeighborIterator nbrIt = this->inBegin((*it)-1);
            NeighborIterator nbrEnd = this->inEnd((*it)-1);
            int deg = 0;
            while(nbrIt!=nbrEnd){
                if(!isMissing(*nbrIt,(*it)-1))
                    deg++;
                nbrIt++;
            }
            *vecIt = deg;
            it++;
            vecIt++;
        }
        return vec;
    }

    /*!
     * out degree
     */
    IntegerVector outdegreeR(IntegerVector node){
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("outdegreeR: range check");
        IntegerVector vec(node.size());
        IntegerVector::iterator vecIt = vec.begin();
        IntegerVector::iterator it = node.begin();
        IntegerVector::iterator end = node.end();
        while(it!=end){
            NeighborIterator nbrIt = this->outBegin((*it)-1);
            NeighborIterator nbrEnd = this->outEnd((*it)-1);
            int deg = 0;
            while(nbrIt!=nbrEnd){
                if(!isMissing((*it)-1,*nbrIt))
                    deg++;
                nbrIt++;
            }
            *vecIt = deg;
            it++;
            vecIt++;
        }
        return vec;
    }
    /*!
     * degree
     */
    IntegerVector degreeR(IntegerVector node){

        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("degreeR: range check");
        IntegerVector vec(node.size());
        IntegerVector::iterator vecIt = vec.begin();
        IntegerVector::iterator it = node.begin();
        IntegerVector::iterator end = node.end();
        while(it!=end){

            NeighborIterator nbrIt = this->begin((*it)-1);
            NeighborIterator nbrEnd = this->end((*it)-1);
            int deg = 0;
            while(nbrIt!=nbrEnd){
                if(!isMissing((*it)-1,*nbrIt))
                    deg++;
                nbrIt++;
            }
            *vecIt = deg;
            it++;
            vecIt++;
        }
        return vec;
    }

    List inneighborsR(IntegerVector node) const{
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("inneighborsR: range check");
        List l;
        for(int i=0;i<node.size();i++){
            IntegerVector vec = engine.template inneighbors<IntegerVector>(node[i]-1);
            for(int j=0;j<vec.size();j++){
                if(isMissing(vec[j], node[i]-1)){
                    vec.erase(j);
                    j--;
                    continue;
                }
                vec[j]++;
            }
            l.push_back(vec);
        }
        return l;
    }

    List outneighborsR(IntegerVector node) const{
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("outneighborsR: range check");
        List l;
        for(int i=0;i<node.size();i++){
            IntegerVector vec = engine.template outneighbors<IntegerVector>(node[i]-1);
            for(int j=0;j<vec.size();j++){
                if(isMissing(node[i]-1,vec[j])){
                    vec.erase(j);
                    j--;
                    continue;
                }
                vec[j]++;
            }
            l.push_back(vec);
        }
        return l;
    }

    List neighborsR(IntegerVector node) const{
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("inneighborsR: range check");
        List l;
        for(int i=0;i<node.size();i++){
            IntegerVector vec = engine.template neighbors<IntegerVector>(node[i]-1);
            for(int j=0;j<vec.size();j++){
                if(isMissing(vec[j], node[i]-1)){
                    vec.erase(j);
                    j--;
                    continue;
                }
                vec[j]++;
            }
            l.push_back(vec);
        }
        return l;
    }

    /*!
     * Sets the (out-)dyads of the specified nodes to missing
     *
     * \param node the nodes
     * \param missing if true the dyads are set to missing,
     * 					otherwise they are set to observed
     */
    void setAllDyadsMissingR1(IntegerVector node,bool missing){
        if(!is_true(all(node>0)) || !is_true(all(node <= size())))
            ::Rf_error("inneighborsR: range check");
        setAllDyadsMissing(as<std::vector<int> >(node),missing);
    }

    /*!
     * sets missingness for all dyads
     * \param missing if true the dyads are set to missing,
     * 					otherwise they are set to observed
     */
    void setAllDyadsMissingR2(bool missing){
        if(missing)
            setAllDyadsMissing();
        else
            setAllDyadsObserved();
    }

    /*!
     * sets all dyads to missing
     */
    void setAllDyadsMissingR3(){
        setAllDyadsMissingR2(true);
    }
};



class Directed{
protected:
    typedef DirectedVertex VertType;
    typedef boost::shared_ptr<VertType> vertPtr;
    typedef boost::shared_ptr< std::vector<ContinAttrib> > cAttrVecPtr;
    typedef boost::shared_ptr< std::vector<DiscreteAttrib> > dAttrVecPtr;
    std::vector< vertPtr > verts;
    cAttrVecPtr contMeta;
    dAttrVecPtr disMeta;
    boost::shared_ptr<double> numEdges;
    void refreshIds(){
        for(int i=0;i<verts.size();i++)
            verts[i]->setId(i);
    }

public:

    typedef Set::const_iterator NeighborIterator;


    Directed(){
        cAttrVecPtr cm(new std::vector<ContinAttrib>());
        dAttrVecPtr dm(new std::vector<DiscreteAttrib>());
        contMeta=cm;
        disMeta=dm;
        numEdges = boost::shared_ptr<double>(new double);
        (*numEdges)=0.0;
    }

    Directed(const Directed& net){
        verts = net.verts;
        contMeta = net.contMeta;
        disMeta = net.disMeta;
        numEdges = net.numEdges;
    }

    Directed(const Directed& net,bool deepCopy){
        if(!deepCopy){
            verts = net.verts;
            contMeta = net.contMeta;
            disMeta = net.disMeta;
            numEdges = net.numEdges;
        }else{
            verts.clear();
            verts.resize(net.verts.size());
            //cout << verts.size()<<"\n"<<net.contMeta<<" "<<net.disMeta;
            //cout << net.verts.size()<<"\n";
            for(int i=0;i<net.verts.size();i++){
                vertPtr v(new VertType(*(net.verts.at(i))));
                verts[i] = v;
            }
            cAttrVecPtr cm(new std::vector<ContinAttrib>(*net.contMeta));
            dAttrVecPtr dm(new std::vector<DiscreteAttrib>(*net.disMeta));
            contMeta = cm;
            disMeta = dm;
            numEdges = boost::shared_ptr<double>(new double);
            (*numEdges) =(*net.numEdges);
        }
    }

    Directed(Rcpp::IntegerMatrix edgeList,int numNodes){
        for(int i=0;i<numNodes;i++){
            vertPtr ver(new DirectedVertex(numNodes));
            verts.push_back(ver);
        }
        numEdges = boost::shared_ptr<double>(new double);
        (*numEdges) = 0.0;
        refreshIds();
        for(int i=0;i<edgeList.nrow();i++){
            //cout<< *numEdges <<" \n "<<this->hasEdge(edgeList(i,0)-1,edgeList(i,1)-1);
            int from = edgeList(i,0)-1;
            int to = edgeList(i,1)-1;
            if(from < 0 || from >= size() || to<0 || to >= size())
                Rf_error("Edgelist indices out of range");
            this->addEdge(from,to);
            //cout <<  (edgeList(i,0)-1) << " " << (edgeList(i,1)-1)<<"\n";
        }
        cAttrVecPtr cm(new std::vector<ContinAttrib>());
        dAttrVecPtr dm(new std::vector<DiscreteAttrib>());
        contMeta = cm;
        disMeta = dm;
    }


    static std::string engineName(){
        return "Directed";
    }

    void addVertex(){
        vertPtr pv(new VertType(size()+1));
        pv->setId(verts.size());
        //TODO: make vetex variable length match
        verts.push_back(pv);
        for(int i=0;i<size();i++){
            verts[i]->setNetworkSize(size());
        }
    }


    void removeVertex(int pos){
        vertPtr pV = verts.at(pos);
        verts.erase(verts.begin() + pos);
        refreshIds();
        (*numEdges) -= pV->indegree();
        (*numEdges) -= pV->outdegree();
        //TODO: correct sets to reflect new ids
        for(int i=0;i<size();i++){
            verts[i]->setNetworkSize(size());
        }
    }

    void reorderVertices(std::vector<int> order){
        std::vector<vertPtr> tmp = verts;
        for(int i=0;i<verts.size();i++)
            verts[i] = tmp[order[i]];
        //TODO correct sets to reflect new ids
    }

    int size() const{
        return verts.size();
    }

    bool hasEdge(int from, int to) const{
        return verts[from]->hasOutedge(to);
    }

    bool removeEdge(int from,int to){
        bool has = false;
        has = verts[from]->removeOutedge(to);
        if(has){
            verts[to]->removeInedge(from);
            (*numEdges)--;
        }
        return has;
    }

    void emptyGraph(){
        for(int i=0;i < verts.size();i++){
            verts[i]->clearOutedges();
            verts[i]->clearInedges();
        }
        (*numEdges) = 0;
    }

    void addEdge(int from,int to){
        if(from==to)
            return;
        if(verts[from]->addOutedge(to)){
            verts[to]->addInedge(from);
            (*numEdges)++;
        }

    }

    bool isMissing(int from,int to) const{
        return verts[from]->isOutmissing(to);
    }

    int nMissing(int from) const{
        return verts[from]->nMissing();
    }


    bool setMissing(int from,int to, bool value){
        if(from==to)
            return false;
        return verts[from]->setOutmissing(to,value);

    }

    void setAllDyadsMissing(){
        for(int i=0;i<size();i++){
            verts[i]->setAllMissing();
        }
    }
    void setAllDyadsObserved(){
        for(int i=0;i<size();i++){
            verts[i]->setAllObserved();
        }
    }

    void setAllDyadsMissing(std::vector<int> nodes,bool missing){
        if(missing){
            for(int i=0;i<nodes.size();i++){
                verts[nodes[i]]->setAllMissing();
            }
        }else{
            for(int i=0;i<nodes.size();i++){
                verts[nodes[i]]->setAllObserved();
            }
        }
    }

    boost::shared_ptr< std::vector< std::pair<int,int> > > missingDyads() const{
        Set s;
        boost::shared_ptr< std::vector< std::pair<int,int> > > vec(new std::vector< std::pair<int,int> >);
        for(int i=0;i<size();i++){
            s = verts[i]->outmissing();
            for(Set::iterator it = s.begin();it!=s.end();it++){
                vec->push_back(std::make_pair(i,*it));
            }
        }
        return vec;
    }

    int nEdges() const{
        return (*numEdges);
    }

    bool isDirected() const{
        return true;
    }

    int indegree(int which) const{
        return verts[which]->indegree();
    }

    template<class Collection>
    Collection inneighbors(int which) const{
        //Set s = verts[which]->inedges();
        return Collection(verts[which]->inedges().begin(),verts[which]->inedges().end());
    }

    NeighborIterator inBegin(int which) const{
        return verts[which]->inedges().begin();
    }

    NeighborIterator inEnd(int which) const{
        return verts[which]->inedges().end();
    }

    const Set& inneighbors(int which) const{
        return verts[which]->inedges();
    }

    int outdegree(int which) const{
        return verts[which]->outdegree();
    }

    template<class Collection>
    Collection outneighbors(int which) const{
        //Set s = verts[which]->outedges();
        return Collection(verts[which]->outedges().begin(),verts[which]->outedges().end());
    }

    NeighborIterator outBegin(int which) const{
        return verts[which]->outedges().begin();
    }

    NeighborIterator outEnd(int which) const{
        return verts[which]->outedges().end();
    }

    const Set& outneighbors(int which) const{
        return verts[which]->outedges();
    }

    int degree(int which) const{
        ::Rf_error("degree not meaningful for directed networks");
        return -1;
    }

    NeighborIterator begin(int which) const{
        ::Rf_error("begin not meaningful for directed networks");
    }

    NeighborIterator end(int which) const{
        ::Rf_error("end not meaningful for directed networks");
    }

    template<class Collection>
    Collection neighbors(int which) const{
        ::Rf_error("neighbors not meaningful for directed networks");
    }

    const Set& neighbors(int which) const{
        ::Rf_error("neighbors not meaningful for directed networks");
        return verts[which]->inedges();
    }

    std::pair<int,int> randomDyad() const{
        int n = size();
        double d1 = Rf_runif(0,(double)n);
        int i1 = floor(d1);
        double d2 = Rf_runif(0,(double)(n-1));
        int i2 = floor(d2);
        if(i2>=i1)
            i2++;
        return std::make_pair(i1,i2);
    }

    void randomDyad(std::pair<int,int>& toggle) const{
        int n = size();
        double d1 = Rf_runif(0,(double)n);
        int i1 = floor(d1);
        double d2 = Rf_runif(0,(double)(n-1));
        int i2 = floor(d2);
        if(i2>=i1)
            i2++;
        toggle.first = i1;
        toggle.second = i2;
    }

    std::pair<int,int> randomEdge() const{
        int n= this->nEdges();
        if(n==0)
            ::Rf_error("randomEdge: network has no edges");
        int c=0;
        int edgeNumber = floor(Rf_runif(0,(double)n));
        int degree;
        for(int i=0;i<verts.size();i++){
            degree = verts[i]->outdegree();
            if(c+degree > edgeNumber){
                const Set* edges = &verts[i]->outedges();
                Set::const_iterator it = edges->begin();
                std::advance(it,edgeNumber-c);
                int to = *(it);
                return std::make_pair(i,to);
            }
            c += degree;
        }
        //cout << "edge: "<<edgeNumber<<" total edges: "<<n<<" c:"<<c<<"\n";
        ::Rf_error("randomEdge: edge not found. ");
        return std::pair<int,int>(-1,-1);
    }

    int randomDyad(int from,bool missing){
        if(!missing){
            int index = floor(Rf_runif(0,size()-1.0));
            if(index>=from)
                index++;
            return index;
        }else
            return verts[from]->randomMissingDyad();
    }

    std::pair<int,int> randomNonEdge() const{
        ::Rf_error("randomNonEdge unimplemented ");

        //the following only works if set is used, not hashSet
        /*		int n= this->nEdges();
		int nVerts = this->size();
		int m= this->maxEdges();
		int c=0;
		int nonEdgeNumber = floor(Rf_runif(0,(double)(m-n)));
		int degree;
		for(int i=0;i<verts.size();i++){
			degree = verts[i]->outdegree();
			if(c+(nVerts-degree) < nonEdgeNumber){
				//find the nth (which) non-edge of vertex i
				int which = nonEdgeNumber - c;
				Set::iterator currentEdge = verts[i]->outedges().begin();
				Set::iterator end = verts[i]->outedges().end();
				while(currentEdge!=end){
					if(which<(*currentEdge)){
						return std::make_pair(i,which);
					}else{
						which++;
						currentEdge++;
					}
				}
				return std::make_pair(i,which);
			}
			c += (nVerts-degree);
		}
		::Rf_error("randomNonEdge: dyad not found");
         */return std::pair<int,int>(-1,-1);
    }


    boost::shared_ptr< std::vector< std::pair<int,int> > > edgelist() const{
        boost::shared_ptr< std::vector< std::pair<int,int> > > v(new std::vector<std::pair<int,int> >());
        v->reserve(nEdges());
        for(int i=0;i<verts.size();i++){
            const Set* out = &(verts[i]->outedges());
            for(Set::const_iterator it = out->begin();it!=out->end();it++){
                std::pair<int,int> p = std::make_pair(i,*it);
                //cout << p.first << " " << p.second<<"\n";
                v->push_back(p);
            }
        }
        return v;
    }

    Rcpp::IntegerMatrix edgelistR(bool includeMissing) const{
        boost::shared_ptr< std::vector< std::pair<int,int> > > v = edgelist();
        int s = 0;
        if(!includeMissing){
            for(int i=0;i<v->size();i++){
                if(!isMissing((*v)[i].first,(*v)[i].second))
                    s++;
            }
        }else
            s = v->size();
        Rcpp::IntegerMatrix rV(s,2);
        int ind = 0;
        for(int i=0;i<v->size();i++){
            if(!includeMissing && isMissing((*v)[i].first,(*v)[i].second))
                continue;
            rV(ind,0) = (*v)[i].first+1;
            rV(ind,1) = (*v)[i].second+1;
            ind++;
        }
        return rV;
    }

    std::vector<std::string> continVarNames() const{
        std::vector<std::string> v(contMeta->size());
        for(int i=0;i<contMeta->size();i++)
            v[i] =contMeta->at(i).getName();
        return v;
    }

    ContinAttrib continVariableAttributes(int which){
        ContinAttrib va = contMeta->at(which);
        return va;
    }

    double continVariableValue(int which,int at) const{
        return verts[at]->continVariable(which);
    }

    void setContinVariableValue(int which,int at,double newValue){
        verts[at]->setContinVariable(newValue, which);
    }


    bool continVariableObserved(int which,int at){
        return verts[at]->continObserved(which);
    }

    std::vector<bool> continVariableObserved(int which){
        std::vector<bool> obs(size(),false);
        for(int i=0;i<size();i++){
            obs[i] = verts[i]->continObserved(which);
        }
        return obs;
    }


    void setContinVariableObserved(int which,int at,bool observed){
        verts[at]->setContinObserved(which,observed);
    }

    void removeContinVariable(int which){
        contMeta->erase(contMeta->begin()+which);
        for(int i=0;i<verts.size();i++){
            verts[i]->removeContinVariable(which);
        }
    }

    void addContinVariable(const std::vector<double>& vals,ContinAttrib& attribs){
        contMeta->push_back(attribs);
        for(int i=0;i<verts.size();i++){
            verts[i]->addContinVariable(vals[i]);
        }
    }

    std::vector<std::string> discreteVarNames() const{
        std::vector<std::string> v(disMeta->size());
        for(int i=0;i<disMeta->size();i++)
            v[i] = disMeta->at(i).getName();
        return v;
    }

    DiscreteAttrib discreteVariableAttributes(int which) const{
        return disMeta->at(which);
    }

    void removeDiscreteVariable(int which){
        disMeta->erase(disMeta->begin()+which);
        for(int i=0;i<verts.size();i++){
            verts[i]->removeDiscreteVariable(which);
        }
    }

    std::vector<int> discreteVariableValues(int which) const{
        std::vector<int> v(verts.size(),0);
        for(int i=0;i<verts.size();i++)
            v[i] = verts[i]->discreteVariable(which);
        return v;
    }

    int discreteVariableValue(int which,int at) const{
        return verts[at]->discreteVariable(which);
    }

    void setDiscreteVariableValue(int which,int at,int newValue){
        verts[at]->setDiscreteVariable(newValue,which);
    }

    std::string discreteVariableLabel(int which,int at) const{
        return disMeta->at(which).labels().at(verts[at]->discreteVariable(which));
    }

    std::vector<std::string> discreteVariable(int which) const{
        std::vector<std::string> v(verts.size(),"");
        for(int i=0;i<verts.size();i++)
            v[i] = disMeta->at(which).labels().at(verts[i]->discreteVariable(which)-1);
        return v;
    }


    void addDiscreteVariable(const std::vector<int>& vals,DiscreteAttrib& attribs){
        disMeta->push_back(attribs);
        for(int i=0;i<verts.size();i++){
            verts[i]->addDiscreteVariable(vals[i]);
        }
    }

    bool discreteVariableObserved(int which,int at){
        return verts[at]->discreteObserved(which);
    }


    std::vector<bool> discreteVariableObserved(int which){
        std::vector<bool> obs(size(),false);
        for(int i=0;i<size();i++){
            obs[i] = verts[i]->discreteObserved(which);
        }
        return obs;
    }


    void setDiscreteVariableObserved(int which,int at,bool observed){
        verts[at]->setDiscreteObserved(which,observed);
    }

    void addDiscreteVariableR(SEXP robj,std::string name){
        std::vector<int> vals;
        std::vector<std::string> levels;
        try{
            Language call("as.factor",robj);
            Rcpp::RObject facSexp = call.eval();
            Language call1("as.integer",facSexp);
            Rcpp::RObject intSexp = call1.eval();
            Language call2("levels",facSexp);
            Rcpp::RObject levelsSexp = call2.eval();
            vals = as< std::vector<int> >(intSexp);
            levels = as< std::vector<std::string> >(levelsSexp);
        } catch(std::exception &ex) {	
            forward_exception_to_r(ex);
        } catch(...){
            ::Rf_error("error, invalid object addDiscreteVariableR");
        }
        if(vals.size() != this->size())
            ::Rf_error("vertex variable size does not match network size");
        std::vector<bool> missing(size(),false);
        for(int i=0;i<size();i++){
            if(vals[i] == NA_INTEGER){
                vals[i] = 1;
                missing[i]=true;
            }
        }
        DiscreteAttrib attrib;
        attrib.setLabels(levels);
        attrib.setName(name);
        addDiscreteVariable(vals,attrib);
        int index = indexOf(name,discreteVarNames());
        for(int i=0;i<missing.size();i++){
            setDiscreteVariableObserved(index,i,!missing[i]);
        }
    }

    Rcpp::RObject getVariableR(std::string name,bool maskMissing){
        std::vector<std::string> names = discreteVarNames();
        int index = -1;
        for(int i=0;i<names.size();i++)
            if(names[i]==name)
                index = i;
        if(index>-1){
            CharacterVector labels = wrap(disMeta->at(index).labels());
            IntegerVector result(size());
            for(int i=0;i<size();i++){
                if(maskMissing && !discreteVariableObserved(index,i)){
                    result[i] = NA_INTEGER;
                }else{
                    result[i] = discreteVariableValue(index,i);
                }
            }
            if(labels.size() > 0){
                result.attr("levels")=labels;
                result.attr("class")="factor";
            }
            return result;
        }
        names = continVarNames();
        for(int i=0;i<names.size();i++)
            if(names[i]==name)
                index = i;
        if(index<0)
            return RObject();
        NumericVector res(size());
        for(int i=0;i<size();i++){
            if(maskMissing && !continVariableObserved(index,i)){
                res[i] = NA_REAL;
            }else{
                res[i] = continVariableValue(index,i);
            }
        }
        if(contMeta->at(index).hasLowerBound()){
            res.attr("lowerBound") = contMeta->at(index).lowerBound();
        }
        if(contMeta->at(index).hasUpperBound()){
            res.attr("upperBound") = contMeta->at(index).upperBound();
        }
        return res;
    }

    void setVariableR(SEXP var, std::string name){
        if(Rf_isNull(var)){
            int index = indexOf(name,discreteVarNames());
            if(index>=0){
                removeDiscreteVariable(index);
            }else{
                index = indexOf(name,continVarNames());
                if(index>=0)
                    removeContinVariable(index);
            }
            return;
        }
        if(Rf_isNumeric(var) && !Rf_isLogical(var)){
            std::vector<bool> missing(size(),false);
            NumericVector vec = as<NumericVector>(var);
            if(vec.size()!=size())
                ::Rf_error("invalid assignment");
            std::vector<double> d(vec.size());
            for(int i=0;i<size();i++){
                if(R_IsNA(vec[i])){
                    d[i]=0.0; //TODO: this should be more intel.
                    missing[i] = true;
                }else
                    d[i]=vec[i];
            }
            int index = indexOf(name,continVarNames());
            if(index>-1)
                removeContinVariable(index);
            index = indexOf(name,discreteVarNames());
            if(index>-1)
                removeDiscreteVariable(index);
            ContinAttrib atr;
            atr.setName(name);
            Rcpp::RObject lb = vec.attr("lowerBound");
            if(!Rf_isNull(lb))
                atr.setLowerBound(as<double>(lb));
            Rcpp::RObject ub = vec.attr("upperBound");
            if(!Rf_isNull(ub))
                atr.setUpperBound(as<double>(ub));
            addContinVariable(d,atr);
            index = indexOf(name,continVarNames());
            for(int i=0;i<missing.size();i++){
                setContinVariableObserved(index,i,!missing[i]);
            }
        }else{
            int index = indexOf(name,continVarNames());
            if(index>-1)
                removeContinVariable(index);
            index = indexOf(name,discreteVarNames());
            if(index>-1)
                removeDiscreteVariable(index);
            addDiscreteVariableR(var,name);
        }
    }

    unsigned64_t maxEdges() const{
        unsigned64_t n = verts.size();
        return n*(n-1);
    }

};


typedef BinaryNet<Directed> DirectedNet;




class Undirected{
protected:
    typedef UndirectedVertex VertType;
    typedef boost::shared_ptr<VertType> vertPtr;
    typedef boost::shared_ptr< std::vector<ContinAttrib> > cAttrVecPtr;
    typedef boost::shared_ptr< std::vector<DiscreteAttrib> > dAttrVecPtr;
    std::vector< vertPtr > verts;
    cAttrVecPtr contMeta;
    dAttrVecPtr disMeta;
    boost::shared_ptr<double> numEdges;
    void refreshIds(){
        for(int i=0;i<verts.size();i++)
            verts[i]->setId(i);
    }

public:

    typedef Set::const_iterator NeighborIterator;

    NeighborIterator getBeginIterator(int node){
        return verts[node]->edges().begin();
    }

    Undirected(){
        cAttrVecPtr cm(new std::vector<ContinAttrib>());
        dAttrVecPtr dm(new std::vector<DiscreteAttrib>());
        contMeta=cm;
        disMeta=dm;
        numEdges = boost::shared_ptr<double>(new double);
        (*numEdges)=0.0;
    }

    Undirected(const Undirected& net){
        verts = net.verts;
        contMeta = net.contMeta;
        disMeta = net.disMeta;
        numEdges = net.numEdges;
    }

    Undirected(const Undirected& net,bool deepCopy){
        if(!deepCopy){
            verts = net.verts;
            contMeta = net.contMeta;
            disMeta = net.disMeta;
            numEdges = net.numEdges;
        }else{
            verts.clear();
            verts.resize(net.verts.size());
            //cout << verts.size()<<"\n"<<net.contMeta<<" "<<net.disMeta;
            //cout << net.verts.size()<<"\n";
            for(int i=0;i<net.verts.size();i++){
                vertPtr v(new VertType(*(net.verts.at(i))));
                verts[i] = v;
            }
            cAttrVecPtr cm(new std::vector<ContinAttrib>(*net.contMeta));
            dAttrVecPtr dm(new std::vector<DiscreteAttrib>(*net.disMeta));
            contMeta = cm;
            disMeta = dm;
            numEdges = boost::shared_ptr<double>(new double);
            (*numEdges) =(*net.numEdges);
        }
    }

    Undirected(Rcpp::IntegerMatrix edgeList,int numNodes){
        for(int i=0;i<numNodes;i++){
            vertPtr ver(new UndirectedVertex(numNodes));
            verts.push_back(ver);
        }
        numEdges = boost::shared_ptr<double>(new double);
        (*numEdges) = 0.0;
        refreshIds();
        for(int i=0;i<edgeList.nrow();i++){
            //cout<< *numEdges <<" \n "<<this->hasEdge(edgeList(i,0)-1,edgeList(i,1)-1);
            int from = edgeList(i,0)-1;
            int to = edgeList(i,1)-1;
            if(from < 0 || from >= size() || to<0 || to >= size())
                Rf_error("Edgelist indices out of range");
            this->addEdge(from,to);
            //cout <<  (edgeList(i,0)-1) << " " << (edgeList(i,1)-1)<<"\n";
        }
        cAttrVecPtr cm(new std::vector<ContinAttrib>());
        dAttrVecPtr dm(new std::vector<DiscreteAttrib>());
        contMeta = cm;
        disMeta = dm;
    }


    static std::string engineName(){
        return "Undirected";
    }

    void addVertex(){
        vertPtr pv(new VertType(size()+1));
        pv->setId(verts.size());
        //TODO: make vetex variable length match
        verts.push_back(pv);
        for(int i=0;i<size();i++){
            verts[i]->setNetworkSize(size());
        }
    }


    void removeVertex(int pos){
        vertPtr pV = verts.at(pos);
        verts.erase(verts.begin() + pos);
        refreshIds();
        (*numEdges) -= pV->degree();
        //TODO: correct sets to reflect new ids
        for(int i=0;i<size();i++){
            verts[i]->setNetworkSize(size());
        }
    }

    void reorderVertices(std::vector<int> order){
        std::vector<vertPtr> tmp = verts;
        for(int i=0;i<verts.size();i++)
            verts[i] = tmp[order[i]];
    }

    int size() const{
        return verts.size();
    }

    bool hasEdge(int from, int to) const{
        return verts[from]->hasEdge(to);
    }

    bool removeEdge(int from,int to){
        bool has = false;
        has = verts[from]->removeEdge(to);
        if(has){
            verts[to]->removeEdge(from);
            (*numEdges)--;
        }
        return has;
    }

    void emptyGraph(){
        for(int i=0;i < verts.size();i++){
            verts[i]->clearEdges();
        }
        (*numEdges) = 0;
    }

    void addEdge(int from,int to){
        if(from==to)
            return;
        if(verts[from]->addEdge(to)){
            verts[to]->addEdge(from);
            (*numEdges)++;
        }

    }

    bool isMissing(int from,int to) const{
        return verts[from]->isMissing(to);
    }

    int nMissing(int from) const{
        return verts[from]->nMissing();
    }

    bool setMissing(int from,int to, bool value){
        if(from==to)
            return false;
        bool wasMissing;
        if(value){
            wasMissing = verts[from]->setMissing(to,true);
            if(!wasMissing)
                verts[to]->setMissing(from,true);
        }else{
            wasMissing = verts[from]->setMissing(to,false);
            if(wasMissing)
                verts[to]->setMissing(from,false);
        }
        return wasMissing;
    }


    void setAllDyadsMissing(){
        for(int i=0;i<size();i++){
            verts[i]->setAllMissing();
        }
    }

    void setAllDyadsObserved(){
        for(int i=0;i<size();i++){
            verts[i]->setAllObserved();
        }
    }

    void setAllDyadsMissing(std::vector<int> nodes,bool missing){
        if(missing){
            for(int i=0;i<nodes.size();i++){
                verts[nodes[i]]->setAllMissing();
                for(int j=0;j<size();j++)
                    if(j!=nodes[i])
                        verts[j]->setMissing(nodes[i],true);
            }
        }else{
            for(int i=0;i<nodes.size();i++){
                verts[nodes[i]]->setAllObserved();
                for(int j=0;j<size();j++)
                    if(j!=nodes[i])
                        verts[j]->setMissing(nodes[i],false);
            }
        }
    }

    boost::shared_ptr< std::vector< std::pair<int,int> > > missingDyads() const{
        Set s;
        boost::shared_ptr< std::vector< std::pair<int,int> > > vec(new std::vector< std::pair<int,int> >);
        for(int i=0;i<size();i++){
            s = verts[i]->missing();
            for(Set::iterator it = s.begin();it!=s.end();it++){
                if(*it>i)
                    vec->push_back(std::make_pair(i,*it));
            }
        }
        return vec;
    }

    int nEdges() const{
        return (*numEdges);
    }

    bool isDirected() const{
        return false;
    }

    int indegree(int which) const{
        ::Rf_error("indegree not meaningful for undirected networks");
        return -1;
    }

    template<class Collection>
    Collection inneighbors(int which) const{
        ::Rf_error("inneighbors not meaningful for undirected networks");
    }

    NeighborIterator inBegin(int which) const{
        ::Rf_error("degree not meaningful for directed networks");
    }

    NeighborIterator inEnd(int which) const{
        ::Rf_error("degree not meaningful for directed networks");
    }

    const Set& inneighbors(int which) const{
        ::Rf_error("inneighbors not meaningful for undirected networks");

    }

    int outdegree(int which) const{
        ::Rf_error("outdegree not meaningful for undirected networks");
        return -1;
    }

    template<class Collection>
    Collection outneighbors(int which) const{
        ::Rf_error("outneighbors not meaningful for undirected networks");
    }

    NeighborIterator outBegin(int which) const{
        ::Rf_error("degree not meaningful for directed networks");
    }

    NeighborIterator outEnd(int which) const{
        ::Rf_error("degree not meaningful for directed networks");
    }


    const Set& outneighbors(int which) const{
        ::Rf_error("outneighbors not meaningful for undirected networks");

    }

    int degree(int which) const{
        return verts[which]->degree();
    }

    template<class Collection>
    Collection neighbors(int which) const{
        return Collection(verts[which]->edges().begin(),verts[which]->edges().end());
    }

    NeighborIterator begin(int which) const{
        return verts[which]->edges().begin();
    }

    NeighborIterator end(int which) const{
        return verts[which]->edges().end();
    }

    const Set& neighbors(int which) const{
        return verts[which]->edges();
    }

    std::pair<int,int> randomDyad() const{
        int n = size();
        double d1 = Rf_runif(0,(double)n);
        int i1 = floor(d1);
        double d2 = Rf_runif(0,(double)(n-1));
        int i2 = floor(d2);
        if(i2>=i1)
            i2++;
        return std::make_pair(i1,i2);
    }

    void randomDyad(std::pair<int,int>& toggle) const{
        int n = size();
        double d1 = Rf_runif(0,(double)n);
        int i1 = floor(d1);
        double d2 = Rf_runif(0,(double)(n-1));
        int i2 = floor(d2);
        if(i2>=i1)
            i2++;
        toggle.first = i1;
        toggle.second = i2;
    }

    std::pair<int,int> randomEdge() const{
        int n= this->nEdges()*2.0;
        if(n==0)
            ::Rf_error("randomEdge: network has no edges");
        int c=0;
        int edgeNumber = floor(Rf_runif(0,(double)n));
        int degree;
        for(int i=0;i<verts.size();i++){
            degree = verts[i]->degree();
            if(c+degree > edgeNumber){
                const Set* edges = &verts[i]->edges();
                Set::const_iterator it = edges->begin();
                std::advance(it,edgeNumber-c);
                int to = *(it);
                return std::make_pair(i,to);
            }
            c += degree;
        }
        //cout << "edge: "<<edgeNumber<<" total edges: "<<n<<" c:"<<c<<"\n";
        ::Rf_error("randomEdge: edge not found. ");
        return std::pair<int,int>(-1,-1);
    }

    int randomDyad(int from,bool missing){
        if(!missing)
            return floor(Rf_runif(0,(double)size()));
        else
            return verts[from]->randomMissingDyad();
    }

    std::pair<int,int> randomNonEdge() const{
        ::Rf_error("randomNonEdge unimplemented ");

        //the following only works if set is used, not hashSet
        /*		int n= this->nEdges();
		int nVerts = this->size();
		int m= this->maxEdges();
		int c=0;
		int nonEdgeNumber = floor(Rf_runif(0,(double)(m-n)));
		int degree;
		for(int i=0;i<verts.size();i++){
			degree = verts[i]->outdegree();
			if(c+(nVerts-degree) < nonEdgeNumber){
				//find the nth (which) non-edge of vertex i
				int which = nonEdgeNumber - c;
				Set::iterator currentEdge = verts[i]->outedges().begin();
				Set::iterator end = verts[i]->outedges().end();
				while(currentEdge!=end){
					if(which<(*currentEdge)){
						return std::make_pair(i,which);
					}else{
						which++;
						currentEdge++;
					}
				}
				return std::make_pair(i,which);
			}
			c += (nVerts-degree);
		}
		::Rf_error("randomNonEdge: dyad not found");
         */return std::pair<int,int>(-1,-1);
    }


    boost::shared_ptr< std::vector< std::pair<int,int> > > edgelist() const{
        boost::shared_ptr< std::vector< std::pair<int,int> > > v(new std::vector<std::pair<int,int> >());
        v->reserve(nEdges());
        for(int i=0;i<verts.size();i++){
            const Set* out = &(verts[i]->edges());
            for(Set::const_iterator it = out->begin();it!=out->end();it++){
                if(*it<i)
                    continue;
                std::pair<int,int> p = std::make_pair(i,*it);
                //cout << p.first << " " << p.second<<"\n";
                v->push_back(p);
            }
        }
        return v;
    }

    Rcpp::IntegerMatrix edgelistR(bool includeMissing) const{
        boost::shared_ptr< std::vector< std::pair<int,int> > > v = edgelist();
        int s = 0;
        if(!includeMissing){
            for(int i=0;i<v->size();i++){
                if(!isMissing((*v)[i].first,(*v)[i].second))
                    s++;
            }
        }else
            s = v->size();
        Rcpp::IntegerMatrix rV(s,2);
        int ind = 0;
        for(int i=0;i<v->size();i++){
            if(!includeMissing && isMissing((*v)[i].first,(*v)[i].second))
                continue;
            rV(ind,0) = (*v)[i].first+1;
            rV(ind,1) = (*v)[i].second+1;
            ind++;
        }
        return rV;
    }

    std::vector<std::string> continVarNames() const{
        std::vector<std::string> v(contMeta->size());
        for(int i=0;i<contMeta->size();i++)
            v[i] =contMeta->at(i).getName();
        return v;
    }

    ContinAttrib continVariableAttributes(int which){
        ContinAttrib va = contMeta->at(which);
        return va;
    }

    double continVariableValue(int which,int at) const{
        return verts[at]->continVariable(which);
    }


    void setContinVariableValue(int which,int at,double newValue){
        verts[at]->setContinVariable(newValue, which);
    }


    bool continVariableObserved(int which,int at){
        return verts[at]->continObserved(which);
    }

    std::vector<bool> continVariableObserved(int which){
        std::vector<bool> obs(size(),false);
        for(int i=0;i<size();i++){
            obs[i] = verts[i]->continObserved(which);
        }
        return obs;
    }


    void setContinVariableObserved(int which,int at,bool observed){
        verts[at]->setContinObserved(which,observed);
    }

    void removeContinVariable(int which){
        contMeta->erase(contMeta->begin()+which);
        for(int i=0;i<verts.size();i++){
            verts[i]->removeContinVariable(which);
        }
    }

    void addContinVariable(const std::vector<double>& vals,ContinAttrib& attribs){
        contMeta->push_back(attribs);
        for(int i=0;i<verts.size();i++){
            verts[i]->addContinVariable(vals[i]);
        }
    }

    std::vector<std::string> discreteVarNames() const{
        std::vector<std::string> v(disMeta->size());
        for(int i=0;i<disMeta->size();i++)
            v[i] = disMeta->at(i).getName();
        return v;
    }

    DiscreteAttrib discreteVariableAttributes(int which) const{
        return disMeta->at(which);
    }

    void removeDiscreteVariable(int which){
        disMeta->erase(disMeta->begin()+which);
        for(int i=0;i<verts.size();i++){
            verts[i]->removeDiscreteVariable(which);
        }
    }

    std::vector<int> discreteVariableValues(int which) const{
        std::vector<int> v(verts.size(),0);
        for(int i=0;i<verts.size();i++)
            v[i] = verts[i]->discreteVariable(which);
        return v;
    }

    int discreteVariableValue(int which,int at) const{
        return verts[at]->discreteVariable(which);
    }

    void setDiscreteVariableValue(int which,int at,int newValue){
        verts[at]->setDiscreteVariable(newValue,which);
    }

    std::string discreteVariableLabel(int which,int at) const{
        return disMeta->at(which).labels().at(verts[at]->discreteVariable(which));
    }

    std::vector<std::string> discreteVariable(int which) const{
        std::vector<std::string> v(verts.size(),"");
        for(int i=0;i<verts.size();i++)
            v[i] = disMeta->at(which).labels().at(verts[i]->discreteVariable(which)-1);
        return v;
    }


    void addDiscreteVariable(const std::vector<int>& vals,DiscreteAttrib& attribs){
        disMeta->push_back(attribs);
        for(int i=0;i<verts.size();i++){
            verts[i]->addDiscreteVariable(vals[i]);
        }
    }

    bool discreteVariableObserved(int which,int at){
        return verts[at]->discreteObserved(which);
    }


    std::vector<bool> discreteVariableObserved(int which){
        std::vector<bool> obs(size(),false);
        for(int i=0;i<size();i++){
            obs[i] = verts[i]->discreteObserved(which);
        }
        return obs;
    }


    void setDiscreteVariableObserved(int which,int at,bool observed){
        verts[at]->setDiscreteObserved(which,observed);
    }

    void addDiscreteVariableR(RObject robj,std::string name){
        std::vector<int> vals;
        std::vector<std::string> levels;
        try{
            Language call("as.factor",robj);
            Rcpp::RObject facSexp = call.eval();
            Language call1("as.integer",facSexp);
            Rcpp::RObject intSexp = call1.eval();
            Language call2("levels",facSexp);
            Rcpp::RObject levelsSexp = call2.eval();
            vals = as< std::vector<int> >(intSexp);
            levels = as< std::vector<std::string> >(levelsSexp);
        }catch(...){
            ::Rf_error("error, invalid object addDiscreteVariableR");
        }
        if(vals.size() != this->size())
            ::Rf_error("vertex variable size does not match network size");

        std::vector<bool> missing(size(),false);
        for(int i=0;i<size();i++){
            if(vals[i] == NA_INTEGER){
                vals[i] = 1;
                missing[i]=true;
            }
        }
        DiscreteAttrib attrib;
        attrib.setLabels(levels);
        attrib.setName(name);
        addDiscreteVariable(vals,attrib);
        int index = indexOf(name,discreteVarNames());
        for(int i=0;i<missing.size();i++){
            setDiscreteVariableObserved(index,i,!missing[i]);
        }
    }

    Rcpp::RObject getVariableR(std::string name,bool maskMissing){
        std::vector<std::string> names = discreteVarNames();
        int index = -1;
        for(int i=0;i<names.size();i++)
            if(names[i]==name)
                index = i;
        if(index>-1){
            CharacterVector labels = wrap(disMeta->at(index).labels());
            IntegerVector result(size());
            for(int i=0;i<size();i++){
                if(maskMissing && !discreteVariableObserved(index,i)){
                    result[i] = NA_INTEGER;
                }else{
                    result[i] = discreteVariableValue(index,i);
                }
            }
            if(labels.size() > 0){
                result.attr("levels")=labels;
                result.attr("class")="factor";
            }
            return result;
        }
        names = continVarNames();
        for(int i=0;i<names.size();i++)
            if(names[i]==name)
                index = i;
        if(index<0)
            return RObject();
        NumericVector res(size());
        for(int i=0;i<size();i++){
            if(maskMissing && !continVariableObserved(index,i)){
                res[i] = NA_REAL;
            }else{
                res[i] = continVariableValue(index,i);
            }
        }
        if(contMeta->at(index).hasLowerBound()){
            res.attr("lowerBound") = contMeta->at(index).lowerBound();
        }
        if(contMeta->at(index).hasUpperBound()){
            res.attr("upperBound") = contMeta->at(index).upperBound();
        }
        return res;
    }

    void setVariableR(SEXP var, std::string name){
        if(Rf_isNull(var)){
            int index = indexOf(name,discreteVarNames());
            if(index>=0){
                removeDiscreteVariable(index);
            }else{
                index = indexOf(name,continVarNames());
                if(index>=0)
                    removeContinVariable(index);
            }
            return;
        }
        if(Rf_isNumeric(var) && !Rf_isLogical(var)){
            std::vector<bool> missing(size(),false);
            NumericVector vec = as<NumericVector>(var);
            ContinAttrib atr;
            atr.setName(name);
            Rcpp::RObject lb = vec.attr("lowerBound");
            if(!Rf_isNull(lb))
                atr.setLowerBound(as<double>(lb));
            Rcpp::RObject ub = vec.attr("upperBound");
            if(!Rf_isNull(ub))
                atr.setUpperBound(as<double>(ub));
            if(vec.size()!=size())
                ::Rf_error("invalid assignment");
            std::vector<double> d(vec.size());
            for(int i=0;i<size();i++){
                if(R_IsNA(vec[i])){
                    if(atr.hasUpperBound() && atr.hasLowerBound()){
                        d[i] = Rf_runif(atr.lowerBound(),atr.upperBound());
                    }else
                        d[i]=0.0;
                    missing[i] = true;
                }else
                    d[i]=vec[i];
                //Rcpp::Rcout<<vec[i] <<" "<< NA_REAL <<" " << missing[i] <<"\n";
            }
            int index = indexOf(name,continVarNames());
            if(index>-1)
                removeContinVariable(index);
            index = indexOf(name,discreteVarNames());
            if(index>-1)
                removeDiscreteVariable(index);
            addContinVariable(d,atr);
            index = indexOf(name,continVarNames());
            for(int i=0;i<missing.size();i++){
                setContinVariableObserved(index,i,!missing[i]);
            }
        }else{
            int index = indexOf(name,continVarNames());
            if(index>-1)
                removeContinVariable(index);
            index = indexOf(name,discreteVarNames());
            if(index>-1)
                removeDiscreteVariable(index);
            addDiscreteVariableR(var,name);
        }
    }

    unsigned64_t maxEdges() const{
        unsigned64_t n= verts.size();
        return n*(n-1LL)/2LL;
    }

};


typedef BinaryNet<Undirected> UndirectedNet;




}

#endif /* NETH_ */
