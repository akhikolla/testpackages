#ifndef STATS_H_
#define STATS_H_

#include "Stat.h"
#include "ParamParser.h"

#include <map>
#include <set>
#include <vector>
#include <utility>
#include <boost/container/flat_map.hpp>
#include <boost/unordered_map.hpp>
namespace lolog{



/*!
 * Counts shared neighbors of two nodes. type == 4
 */
template<class Engine>
int sharedNbrs(const BinaryNet<Engine>& net, int from, int to){
     return sharedNbrs(net, from, to, 4);
}


/*!
 * Counts shared neighbors of two nodes.
 * type = 1     :   from -> to -> nbr -> from
 * type = 2     :   from -> to <- nbr <- from (homogeneous)
 * type = 3     :   either type 1 or 2
 * type = 4     :   all combinations
 */
template<class Engine>
int sharedNbrs(const BinaryNet<Engine>& net, int from, int to, int type){
    if(net.isDirected()){
        if(type == 4)
            return allDirectedSharedNbrs(net, from, to);
        return directedSharedNbrs(net, from, to, type);
    }
    return undirectedSharedNbrs(net, from, to);
}



template<class Engine>
int undirectedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    NeighborIterator fit = net.begin(from);
    NeighborIterator fend = net.end(from);
    NeighborIterator tit = net.begin(to);
    NeighborIterator tend = net.end(to);
    int shared = 0;
    while(tit!=tend && fit!=fend){
        if(*tit==*fit){
            shared++;
            tit++;
            fit++;
        }else if(*tit<*fit){
            tit++;
        }else
            fit++;
    }
    return shared;
}

template<class Engine>
int allDirectedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    NeighborIterator ifit = net.inBegin(from);
    NeighborIterator ifend = net.inEnd(from);
    NeighborIterator ofit = net.outBegin(from);
    NeighborIterator ofend = net.outEnd(from);
    int shared = 0;
    while(ifit != ifend){
        shared += net.hasEdge(*ifit, to);
        shared += net.hasEdge(to, *ifit);
        ifit++;
    }
    while(ofit != ofend){
        shared += net.hasEdge(*ofit, to);
        shared += net.hasEdge(to, *ofit);
        ofit++;
    }
    return shared;
}


/*!
 * type = 1     :   from -> to -> nbr -> from
 * type = 2     :   from -> to <- nbr <- from (homogeneous)
 * type = 3     :   either type 1 or 2
 */
template<class Engine>
int directedSharedNbrs(const BinaryNet<Engine>& net, int from, int to, int type){
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    NeighborIterator fit, fend, tit, tend;

    int sn = 0;
    if(type == 1 || type == 3){
        fit = net.inBegin(from);
        fend = net.inEnd(from);
        tit = net.outBegin(to);
        tend = net.outEnd(to);
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                sn++;
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
    }
    if(type == 2 || type == 3){
        fit = net.outBegin(from);
        fend = net.outEnd(from);
        tit = net.inBegin(to);
        tend = net.inEnd(to);
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                if(type == 3){
                    bool counted = net.hasEdge(to, *tit) && net.hasEdge(*tit, from);
                    if(!counted)
                        sn++;
                }else
                    sn++;
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
    }
    return sn;
}


/*!
 * the number of edges in the network
 */
template<class Engine>
class Edges : public BaseStat<Engine>{
public:
    Edges(){
    }

    /*!
     * constructor. params is unused
     */
    Edges(List params){
    }


    std::string name(){
        return "edges";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"edges");
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        this->initSingle(net.nEdges());
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        BaseOffset<Engine>::update(net.hasEdge(from,to) ? -1.0 : 1.0, 0);
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }
};

typedef Stat<Directed, Edges<Directed> > DirectedEdges;
typedef Stat<Undirected, Edges<Undirected> > UndirectedEdges;


/*!
 * in/out/k-Stars
 */
template<class Engine>
class Star : public BaseStat< Engine > {
protected:
    std::vector<int> starDegrees; /*!< the star degrees */
    EdgeDirection direction;
public:

    Star(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas = t;
        direction=IN;
    }

    /*!
     * \param params 	a list
     */
    Star(List params){
        ParamParser p(name(), params);
        starDegrees = p.parseNext< std::vector<int> >("k");
        direction = p.parseNextDirection("direction", IN);
        p.end();

        this->init(starDegrees.size());
    }


    std::string name(){
        return "star";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<starDegrees.size();i++){
            int d = starDegrees[i];
            std::string nm = "star."+asString(d);
            if(direction == IN)
                nm = "in-" + nm;
            if(direction == OUT)
                nm = "out-" + nm;
            statnames.push_back(nm);
        }
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        if(!net.isDirected())
            direction = UNDIRECTED;
        std::vector<double> v(starDegrees.size(),0.0);
        this->lastStats = std::vector<double>(starDegrees.size(),0.0);
        for(int i=0; i<net.size();i++){
            double nEd;
            if(!net.isDirected())
                nEd = net.degree(i);
            else{
                if(direction == IN)
                    nEd = net.indegree(i);
                else
                    nEd = net.outdegree(i);
            }
            for(int j=0;j<starDegrees.size();j++){
                v[j] += nchoosek(nEd,(double)starDegrees[j]);
                //cout << "n:"<<nEd<<" s:"<<starDegrees[j]<<" nck:"<<nchoosek(nEd,(double)starDegrees[j])<<"\n" ;
            }
        }
        this->stats=v;
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        int n;
        if(!net.isDirected())
            n = net.degree(to);
        else{
            if(direction==IN)
                n = net.indegree(to);
            else
                n = net.outdegree(from);
        }
        bool edge = net.hasEdge(from,to);
        for(int i=0;i<starDegrees.size();i++){
            if(edge){
                BaseOffset<Engine>::update(-nchoosek(n,starDegrees[i]) +
                        nchoosek(n-1.0,starDegrees[i]),i);

            }else{
                BaseOffset<Engine>::update(nchoosek(n+1.0,starDegrees[i])-nchoosek(n,starDegrees[i]), i);
            }
        }
        if(!net.isDirected()){
            n = net.degree(from);
            edge = net.hasEdge(from,to);
            for(int i=0;i<starDegrees.size();i++){
                if(edge){
                    BaseOffset<Engine>::update(-nchoosek(n,starDegrees[i]) +
                            nchoosek(n-1.0,starDegrees[i]), i);

                }else{
                    BaseOffset<Engine>::update( nchoosek(n+1.0,starDegrees[i])-nchoosek(n,starDegrees[i]), i);
                }
            }
        }
    }

    bool isOrderIndependent(){
        return true;
    }

};

typedef Stat<Directed, Star<Directed> > DirectedStar;
typedef Stat<Undirected, Star<Undirected> > UndirectedStar;


/*!
 * The number of triangles in the network.
 */
template<class Engine>
class Triangles : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
public:


    Triangles(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas = t;
    }
    Triangles(List params){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas = t;
    }

    std::string name(){
        return "triangles";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"triangles");
        return statnames;
    }
    int sharedNbrs(const BinaryNet<Engine>& net, int from, int to){
        if(net.isDirected()){
            return directedSharedNbrs(net, from, to);
        }
        return undirectedSharedNbrs(net, from, to);
    }
    int undirectedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
        NeighborIterator fit = net.begin(from);
        NeighborIterator fend = net.end(from);
        NeighborIterator tit = net.begin(to);
        NeighborIterator tend = net.end(to);
        int shared = 0;
        while(tit!=tend && fit!=fend){
            if(*tit==*fit){
                shared++;
                tit++;
                fit++;
            }else if(*tit<*fit){
                tit++;
            }else
                fit++;
        }
        return shared;
    }

    int directedSharedNbrs(const BinaryNet<Engine>& net, int from, int to){
        NeighborIterator ifit = net.inBegin(from);
        NeighborIterator ifend = net.inEnd(from);
        NeighborIterator ofit = net.outBegin(from);
        NeighborIterator ofend = net.outEnd(from);
        int shared = 0;
        while(ifit != ifend){
            shared += net.hasEdge(*ifit, to);
            shared += net.hasEdge(to, *ifit);
            ifit++;
        }
        while(ofit != ofend){
            shared += net.hasEdge(*ofit, to);
            shared += net.hasEdge(to, *ofit);
            ofit++;
        }
        return shared;
    }


    void calculate(const BinaryNet<Engine>& net){
        this->initSingle(0.0);
        double sumTri = 0.0;

        boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

        std::vector< std::pair<int,int> >::iterator it = edges->begin();
        while(it != edges->end()){
            int shared = sharedNbrs(net, (*it).first,(*it).second);
            sumTri += shared;
            it++;
        }
        sumTri = sumTri/3.0;
        this->stats[0] = sumTri;//sumSqrtTri - sumSqrtExpected;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        int shared = sharedNbrs(net, from, to);
        bool hasEdge = net.hasEdge(from,to);
        if(hasEdge){
            BaseOffset<Engine>::update(-shared,0);
            //sumTri -= shared;
        }else{
            BaseOffset<Engine>::update(shared,0);
            //sumTri += shared;
        }
        //this->stats[0] = sumTri;//sumSqrtTri - sumSqrtExpected;
    }

    bool isOrderIndependent(){
        return true;
    }
};

typedef Stat<Directed, Triangles<Directed> > DirectedTriangles;
typedef Stat<Undirected, Triangles<Undirected> > UndirectedTriangles;


template<class Engine>
class Clustering : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double triangles;
    double twostars;

    double lastTriangles;
    double lastTwostars;

public:

    Clustering(){
        lastTwostars = lastTriangles = twostars = triangles = 0.0;
    }

    /*!
     * \param params
     */
    Clustering(List params){
        lastTwostars = lastTriangles = twostars = triangles = 0.0;
    }

    std::string name(){
        return "clustering";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"clustering");
        return statnames;
    }

    void resetMemory(){
        lastTriangles = triangles;
        lastTwostars = twostars;
    }

    void rollback(const BinaryNet<Engine>& net){
        BaseOffset<Engine>::rollback(net);
        twostars = lastTwostars;
        triangles = lastTriangles;
    }


    void calculate(const BinaryNet<Engine>& net){
        int nstats = 1;

        this->init(nstats);
        triangles = twostars = 0.0;
        boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

        std::vector< std::pair<int,int> >::iterator it = edges->begin();
        while(it != edges->end()){
            int shared = sharedNbrs(net, (*it).first,(*it).second);
            triangles += shared;
            it++;
        }
        triangles = triangles/3.0;

        twostars = 0.0;
        for(int i=0; i<net.size();i++){
            double nEd = net.degree(i);
            twostars += nchoosek(nEd,2.0);
        }

        this->stats[0] = 3.0 * triangles / twostars;
        if(twostars < 0.5 )
            this->stats[0] = 0.0;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        resetMemory();
        int shared = sharedNbrs(net, from, to);
        bool hasEdge = net.hasEdge(from,to);
        if(hasEdge){
            triangles -= shared;
        }else{
            triangles += shared;
        }


        int n = net.degree(to);

        if(hasEdge){
            twostars += -nchoosek(n,2.0) + nchoosek(n-1.0,2.0);

        }else{
            twostars += nchoosek(n+1.0,2.0)-nchoosek(n,2.0);
        }

        if(!net.isDirected()){
            n = net.degree(from);
            if(hasEdge){
                twostars += -nchoosek(n,2.0) + nchoosek(n-1.0,2.0);

            }else{
                twostars += nchoosek(n+1.0,2.0)-nchoosek(n,2.0);
            }
        }
        this->stats[0] = 3.0 * triangles / twostars;
        if(twostars < 0.5 )
            this->stats[0] = 0.0;
    }

    bool isOrderIndependent(){
        return true;
    }

};

//typedef Stat<Directed, Clustering<Directed> > DirectedClustering;
typedef Stat<Undirected, Clustering<Undirected> > UndirectedClustering;



template<class Engine>
class Transitivity : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double triads;
    double nPosTriads;

    double lastTriads;
    double lastNPosTriads;

public:

    Transitivity(){
        triads = lastTriads = nPosTriads = lastNPosTriads = 0.0;
    }

    /*!
     * \param params
     */
    Transitivity(List params){
        triads = lastTriads = nPosTriads = lastNPosTriads = 0.0;
    }

    std::string name(){
        return "transitivity";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"transitivity");
        return statnames;
    }

    void resetMemory(){
        lastTriads = triads;
        lastNPosTriads = nPosTriads;
    }

    void rollback(const BinaryNet<Engine>& net){
        BaseOffset<Engine>::rollback(net);
        nPosTriads = lastNPosTriads;
        triads = lastTriads;
    }


    void calculate(const BinaryNet<Engine>& net){
        int nstats = 1;

        this->init(nstats);
        triads = nPosTriads = 0.0;
        boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

        std::vector< std::pair<int,int> >::iterator it = edges->begin();
        while(it != edges->end()){
            int shared = sharedNbrs(net, (*it).first,(*it).second);
            triads += shared;
            nPosTriads += std::min(net.degree((*it).first), net.degree((*it).second)) - 1.0;
            it++;
        }
        this->stats[0] = (1.0 + triads) / (1.0 + nPosTriads);
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        resetMemory();
        int shared = sharedNbrs(net, from, to);
        bool hasEdge = net.hasEdge(from,to);
        int change = hasEdge ? -1 : 1;
        int fromDeg = net.degree(from);
        int toDeg = net.degree(to);
        triads += change * 3.0 * shared;
        NeighborIterator it = net.begin(from);
        NeighborIterator end = net.end(from);
        while(it != end){

            if(*it != to && net.degree(*it) >= (fromDeg + !hasEdge)){
                nPosTriads += change;
            }
            it++;
        }
        it = net.begin(to);
        end = net.end(to);
        while(it != end){
            if(*it != from && net.degree(*it) >= (toDeg + !hasEdge)){
                nPosTriads += change;
            }
            it++;
        }
        nPosTriads += change * (std::min(fromDeg, toDeg) + !hasEdge - 1.0);

        this->stats[0] =  (1.0 + triads) / (1.0 + nPosTriads);
    }

    bool isOrderIndependent(){
        return true;
    }
};
typedef Stat<Directed, Transitivity<Directed> > DirectedTransitivity;
typedef Stat<Undirected, Transitivity<Undirected> > UndirectedTransitivity;

/*!
 * the number of reciprocal edges in the network
 */
template<class Engine>
class Mutual : public BaseStat< Engine > {
public:
    Mutual(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats = v;
        this->thetas = t;
    }

    /*!
     * constructor. params is unused
     */
    Mutual(List params){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats = v;
        this->thetas = t;
    }

    std::string name(){
        return "mutual";
    }
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"mutual");
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        this->init(1);
        if(!net.isDirected())
            Rf_error("Mutual only make sense for directed networks");

        double rec = 0.0;
        int from, to;
        boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
        for(int i=0;i<edges->size();i++){
            from = (*edges)[i].first;
            to = (*edges)[i].second;
            if(from<to && net.hasEdge(to,from))
                rec++;
        }
        std::vector<double> v(1,rec);
        this->stats=v;
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        bool addingEdge = !net.hasEdge(from,to);
        bool hasReverse = net.hasEdge(to,from);
        double change;
        if(addingEdge && hasReverse)
            change = 1.0;
        else if(!addingEdge && hasReverse)
            change = -1.0;
        else
            change = 0.0;
        BaseOffset<Engine>::update(change,0);
    }

    bool isOrderIndependent(){
        return true;
    }
};

typedef Stat<Directed, Mutual<Directed> > DirectedMutual;


/*!
 * Adds a statistic for each edge in which the 'to' and 'from' nodes
 * match on a variable
 */
template<class Engine>
class NodeMatch : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    std::string variableName; /*!< the name of the matching variable */
    int varIndex; /*!< the index of the variable in the network */
    int nstats; /*!< the number of stats generated (i.e. the number of levels squared) */
    int nlevels; /*!< the number of levels of the variable */
public:
    NodeMatch(){
        variableName="";
        nstats=nlevels=varIndex = -1;
    }

    NodeMatch(std::string name){
        variableName=name;
        nstats=nlevels=varIndex = -1;
    }

    NodeMatch(List params){
        nstats=nlevels=varIndex = -1;
        ParamParser p(name(), params);
        variableName = p.parseNext< std::string >("name");
        p.end();
    }


    std::string name(){
        return "nodeMatch";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"nodematch."+variableName);
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        int from,to;
        int value1, value2;
        std::vector<std::string> vars = net.discreteVarNames();
        int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName){
                variableIndex = i;
            }
        }
        if(variableIndex<0)
            ::Rf_error("NodeMatch::calculate nodal attribute not found in network");
        varIndex = variableIndex;
        //nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
        //nstats = nlevels*nlevels;
        nstats = 1;
        this->init(nstats);
        boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
        for(int i=0;i<edges->size();i++){
            from = (*edges)[i].first;
            to = (*edges)[i].second;
            value1 = net.discreteVariableValue(varIndex,from) - 1;
            value2 = net.discreteVariableValue(varIndex,to) - 1;
            //this->stats[value1 + nlevels*value2]++;
            if(value1==value2)
                this->stats[0]++;
        }

    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        bool addingEdge = !net.hasEdge(from,to);
        int value1 = net.discreteVariableValue(varIndex,from) - 1;
        int value2 = net.discreteVariableValue(varIndex,to) - 1;
        if(value1==value2){
            if(addingEdge)
                BaseOffset<Engine>::update(1.0,0);//this->stats[0]++;
            else
                BaseOffset<Engine>::update(-1.0,0);//this->stats[0]--;
        }
    }

    void discreteVertexUpdate(const BinaryNet<Engine>& net, const  int& vert,
            const int& variable, const  int& newValue, const  std::vector<int> &order, const  int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        if(variable != varIndex)
            return;
        int val = net.discreteVariableValue(varIndex,vert);
        if(net.isDirected()){
            NeighborIterator it = net.outBegin(vert);
            NeighborIterator end = net.outEnd(vert);
            while(it!=end){
                int val2 = net.discreteVariableValue(varIndex,*it);
                if(val2==val)
                    BaseOffset<Engine>::update(-1.0,0);//this->stats[0]--;
                if(val2==newValue)
                    BaseOffset<Engine>::update(1.0,0);//this->stats[0]++;
                it++;
            }
            it = net.inBegin(vert);
            end = net.inEnd(vert);
            while(it!=end){
                int val2 = net.discreteVariableValue(varIndex,*it);
                if(val2==val)
                    BaseOffset<Engine>::update(-1.0,0);//this->stats[0]--;
                if(val2==newValue)
                    BaseOffset<Engine>::update(1.0,0);//this->stats[0]++;
                it++;
            }
        }else{
            NeighborIterator it = net.begin(vert);
            NeighborIterator end = net.end(vert);
            while(it!=end){
                int val2 = net.discreteVariableValue(varIndex,*it);
                if(val2==val)
                    BaseOffset<Engine>::update(-1.0,0);//this->stats[0]--;
                if(val2==newValue)
                    BaseOffset<Engine>::update(1.0,0);//this->stats[0]++;
                it++;
            }
        }
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};

typedef Stat<Directed, NodeMatch<Directed> > DirectedNodeMatch;
typedef Stat<Undirected, NodeMatch<Undirected> > UndirectedNodeMatch;



/*!
 * Adds a statistic for each edge in which the 'to' and 'from' nodes
 * match on a variable
 */
template<class Engine>
class NodeMix : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    std::string variableName; /*!< the name of the matching variable */
    int varIndex; /*!< the index of the variable in the network */
    int nstats; /*!< the number of stats generated (i.e. the number of levels squared) */
    int nlevels; /*!< the number of levels of the variable */
    std::vector<std::string> levels;
public:
    NodeMix(){
        variableName="";
        nstats=nlevels=varIndex = -1;
    }

    NodeMix(std::string name){
        variableName=name;
        nstats=nlevels=varIndex = -1;
    }

    NodeMix(List params){
        nstats=nlevels=varIndex = -1;
        ParamParser p(name(), params);
        variableName = p.parseNext< std::string >("name");
        p.end();
    }


    std::string name(){
        return "nodeMix";
    }

    int getIndex(int i,int j){
        int c;
        if(i>j){
            c=i;
            i=j;
            j=c;
        }
        c = 0;
        for(int k=1;k<=i;k++){
            c +=  nlevels - k;
        }
        return c + j;
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(nstats,"");
        for(int i=0;i<levels.size();i++){
            for(int j=i;j<levels.size();j++){
                std::string name = "nodemix." + levels.at(j) + "." + levels.at(i);
                statnames.at(getIndex(i,j)) = name;
            }
        }
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        int from,to;
        int value1, value2;
        std::vector<std::string> vars = net.discreteVarNames();
        int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName){
                variableIndex = i;
            }
        }
        if(variableIndex<0)
            ::Rf_error("NodeMatch::calculate nodal attribute not found in network");
        varIndex = variableIndex;
        levels = net.discreteVariableAttributes(varIndex).labels();
        nlevels = levels.size();
        nstats = nlevels * (nlevels + 1) / 2;
        this->init(nstats);
        boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
        for(int i=0;i<edges->size();i++){
            from = (*edges)[i].first;
            to = (*edges)[i].second;
            value1 = net.discreteVariableValue(varIndex,from) - 1;
            value2 = net.discreteVariableValue(varIndex,to) - 1;
            //this->stats[value1 + nlevels*value2]++;
            this->stats[getIndex(value1,value2)]++;

        }

    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        bool addingEdge = !net.hasEdge(from,to);
        double change = addingEdge ? 1.0 : -1.0;
        int value1 = net.discreteVariableValue(varIndex,from) - 1;
        int value2 = net.discreteVariableValue(varIndex,to) - 1;
        BaseOffset<Engine>::update(change, getIndex(value1,value2));//this->stats[getIndex(value1,value2)] += change;
    }

    void discreteVertexUpdate(const BinaryNet<Engine>& net, const  int& vert,
            const int& variable, const  int& newValue, const  std::vector<int> &order, const  int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        if(variable != varIndex)
            return;
        Rf_error("NodeMix unimplemented");
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }
};

typedef Stat<Directed, NodeMix<Directed> > DirectedNodeMix;
typedef Stat<Undirected, NodeMix<Undirected> > UndirectedNodeMix;







/*!
 * Adds a stat for the counts of degrees
 */
template<class Engine>
class Degree : public BaseStat< Engine > {
protected:
    EdgeDirection direction;
    std::vector<int> degrees;
    bool lessThanOrEqual;
public:

    Degree(){
        direction = UNDIRECTED;
        lessThanOrEqual = false;
    }

    Degree(std::vector<int> deg){
        direction = UNDIRECTED;
        degrees = deg;
        lessThanOrEqual = false;
    }

    /*!
     * \param params 	a list of length 1, the first element of which is an integer vector
     * 					of degrees
     */
    Degree(List params){
        ParamParser p(name(), params);
        degrees = p.parseNext< std::vector<int> >("d");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        lessThanOrEqual = p.parseNext("lessThanOrEqual", false);
        p.end();
    }


    std::string name(){
        return "degree";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<degrees.size();i++){
            int d = degrees[i];
            std::string nm = "degree."+asString(d);
            if(direction == IN)
                nm = "in-" + nm;
            if(direction == OUT)
                nm = "out-" + nm;
            statnames.push_back(nm);
        }
        return statnames;
    }


    void calculate(const BinaryNet<Engine>& net){
        int nstats = degrees.size();
        this->init(nstats);
        double n = net.size();
        for(int i=0;i<n;i++){
            for(int j=0;j<nstats;j++){
                if(net.isDirected()){
                    if(direction==UNDIRECTED){
                        this->stats[j] += comp((net.outdegree(i) + net.indegree(i)), degrees[j]);
                    }else if(direction==OUT)
                        this->stats[j] += comp(net.outdegree(i), degrees[j]);
                    else if(direction==IN)
                        this->stats[j] += comp(net.indegree(i), degrees[j]);
                }else{
                    this->stats[j] += comp(net.degree(i), degrees[j]);
                }
            }
        }
        //std::cout << lessThanOrEqual << " \n";
    }

    inline bool comp(int d,int d1){
        return lessThanOrEqual ? d <= d1 : d == d1;
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        int change = !net.hasEdge(from,to) ? 1 : -1;
        int fromDegree = 0;
        int fromDegreeNew = 0;
        int toDegree = 0;
        int toDegreeNew = 0;
        if(net.isDirected()){
            if(direction==UNDIRECTED){
                fromDegree = net.outdegree(from) + net.indegree(from);
                toDegree = net.outdegree(to) + net.indegree(to);
                fromDegreeNew = change;
                toDegreeNew = change;
            }else if(direction==OUT){
                fromDegree = net.outdegree(from);
                toDegree = net.outdegree(to);
                fromDegreeNew = change;
            }else if(direction==IN){
                fromDegree += net.indegree(from);
                toDegree += net.indegree(to);
                toDegreeNew = change;
            }
        }else{
            fromDegree = net.degree(from);
            toDegree = net.degree(to);
            fromDegreeNew = change;
            toDegreeNew = change;
        }
        toDegreeNew += toDegree;
        fromDegreeNew += fromDegree;

        for(int j=0;j<degrees.size();j++){
            if(comp(fromDegree, degrees[j]))
                BaseOffset<Engine>::update(-1.0,j);//this->stats[j]--;
            if(comp(toDegree, degrees[j]))
                BaseOffset<Engine>::update(-1.0,j);//this->stats[j]--;
            if(comp(fromDegreeNew, degrees[j]))
                BaseOffset<Engine>::update(1.0,j);//this->stats[j]++;
            if(comp(toDegreeNew, degrees[j]))
                BaseOffset<Engine>::update(1.0,j);//this->stats[j]++;
        }

    }

    bool isOrderIndependent(){
        return true;
    }
};

typedef Stat<Directed, Degree<Directed> > DirectedDegree;
typedef Stat<Undirected, Degree<Undirected> > UndirectedDegree;


template<class Engine>
class DegreeCrossProd : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double nEdges;
    double crossProd;

    double lastNEdges;
    double lastCrossProd;

public:

    DegreeCrossProd(){
        lastNEdges = lastCrossProd = crossProd = nEdges = 0.0;
    }

    /*!
     * \param params
     */
    DegreeCrossProd(List params){
        lastNEdges = lastCrossProd = nEdges = crossProd = 0.0;
    }

    std::string name(){
        return "degreeCrossProd";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"degreeCrossProd");
        return statnames;
    }

    void resetMemory(){
        lastNEdges = nEdges;
        lastCrossProd = crossProd;
    }

    void rollback(const BinaryNet<Engine>& net){
        BaseOffset<Engine>::rollback(net);
        nEdges = lastNEdges;
        crossProd = lastCrossProd;
    }


    void calculate(const BinaryNet<Engine>& net){
        int nstats = 1;

        this->init(nstats);
        nEdges = net.nEdges();
        crossProd = 0.0;
        boost::shared_ptr<std::vector< std::pair<int,int> > > edges = net.edgelist();

        std::vector< std::pair<int,int> >::iterator it = edges->begin();
        while(it != edges->end()){
            crossProd += net.degree((*it).first) * net.degree((*it).second);
            it++;
        }
        if(nEdges==0)
            this->stats[0] = 0;
        else
            this->stats[0] = crossProd / nEdges;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        resetMemory();
        double toDeg;
        double fromDeg;
        bool addingEdge = !net.hasEdge(from,to);
        double edgeChange = 2.0*(addingEdge - 0.5);

        if(addingEdge)
            crossProd += (net.degree(from) + 1.0) * (net.degree(to) + 1.0);
        else
            crossProd -= net.degree(from) * net.degree(to);

        NeighborIterator it = net.begin(from);
        NeighborIterator end = net.end(from);
        double deg = net.degree(from);
        while(it!=end){
            double deg2 = net.degree(*it);
            if(addingEdge)
                crossProd += deg2;//(deg+1.0)*deg2 - deg*deg2
            else if(*it != to)
                crossProd -= deg2;// (deg-1.0)*deg - deg*deg2;

            it++;
        }

        it = net.begin(to);
        end = net.end(to);
        deg = net.degree(to);
        while(it!=end){
            double deg2 = net.degree(*it);
            if(addingEdge)
                crossProd += deg2;//(deg+1.0)*deg2 - deg*deg2
            else if(*it != from)
                crossProd -= deg2;// (deg-1.0)*deg - deg*deg2;

            it++;
        }
        nEdges += edgeChange;
        if(nEdges==0)
            BaseOffset<Engine>::update(-this->stats[0],0);//this->stats[0] = 0;
        else
            BaseOffset<Engine>::update(crossProd / nEdges - this->stats[0], 0);//this->stats[0] = crossProd / nEdges;
    }

    bool isOrderIndependent(){
        return true;
    }

};

typedef Stat<Directed, DegreeCrossProd<Directed> > DirectedDegreeCrossProd;
typedef Stat<Undirected, DegreeCrossProd<Undirected> > UndirectedDegreeCrossProd;


/*!
 * Main effect of a continuous covariate.
 */
template<class Engine>
class NodeCov : public BaseStat< Engine > {
protected:
    EdgeDirection direction;
    std::string variableName;
    int varIndex;
    bool isDiscrete;
public:

    NodeCov(){
        varIndex =  0;
        direction = UNDIRECTED;
        isDiscrete = false;
    }

    NodeCov(std::string name,EdgeDirection d){
        varIndex = 0;
        direction = d;
        variableName = name;
        isDiscrete = false;
    }

    NodeCov(std::string name){
        varIndex = 0;
        direction = UNDIRECTED;
        variableName = name;
        isDiscrete = false;
    }


    NodeCov(List params){
        varIndex = 0;
        isDiscrete=false;

        ParamParser p(name(), params);
        variableName = p.parseNext< std::string >("name");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        p.end();

    }

    std::string name(){
        return "nodeCov";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        statnames.assign(1,"nodecov."+variableName);
        return statnames;
    }

    double getValue(const BinaryNet<Engine>& net, int ind){
        double val;
        if(isDiscrete)
            val = net.discreteVariableValue(varIndex,ind);
        else
            val = net.continVariableValue(varIndex,ind);
        return val;
    }

    void calculate(const BinaryNet<Engine>& net){
        isDiscrete = false;
        std::vector<std::string> vars = net.continVarNames();
        int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName){
                variableIndex = i;
            }
        }
        if(variableIndex == -1){
            isDiscrete = true;
            vars = net.discreteVarNames();
            for(int i=0;i<vars.size();i++){
                if(vars[i] == variableName){
                    variableIndex = i;
                }
            }
        }
        if(variableIndex<0)
            ::Rf_error("nodal attribute not found in network");
        varIndex = variableIndex;
        int nstats = 1;
        this->init(nstats);
        this->stats[0] = 0;
        for(int i=0;i<net.size();i++){
            double val = getValue(net,i);
            if(net.isDirected()){
                if(direction == IN || direction == UNDIRECTED)
                    this->stats[0] += val * net.indegree(i);
                if(direction == OUT || direction == UNDIRECTED)
                    this->stats[0] += val * net.outdegree(i);
            }else{
                this->stats[0] += val * net.degree(i);
            }
        }
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        if(net.isDirected()){
            if(direction == IN || direction == UNDIRECTED)
                this->stats[0] += change * getValue(net,to);
            if(direction == OUT || direction == UNDIRECTED)
                this->stats[0] += change * getValue(net,from);
        }else{
            this->stats[0] += change * (getValue(net,to)+getValue(net,from));
        }
    }

    void discreteVertexUpdate(const BinaryNet<Engine>& net, const  int& vert,
            const int& variable, const  int& newValue, const  std::vector<int> &order, const  int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        if(isDiscrete && variable==varIndex){
            double oldValue = getValue(net,vert);
            int deg = 0;
            if(net.isDirected()){
                if(direction == IN || direction == UNDIRECTED)
                    deg += net.indegree(vert);
                if(direction == OUT || direction == UNDIRECTED)
                    deg += net.outdegree(vert);
            }else
                deg = net.degree(vert);
            this->stats[0] += deg*(newValue - oldValue);
        }
    }

    void continVertexUpdate(const BinaryNet<Engine>& net, const int& vert,
            const int& variable, const double& newValue, const std::vector<int> &order, const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        if(!isDiscrete && variable==varIndex){
            double oldValue = getValue(net,vert);
            int deg = 0;
            if(net.isDirected()){
                if(direction == IN || direction == UNDIRECTED)
                    deg += net.indegree(vert);
                if(direction == OUT || direction == UNDIRECTED)
                    deg += net.outdegree(vert);
            }else
                deg = net.degree(vert);
            this->stats[0] += deg*(newValue - oldValue);
        }
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};



typedef Stat<Directed, NodeCov<Directed> > DirectedNodeCov;
typedef Stat<Undirected, NodeCov<Undirected> > UndirectedNodeCov;



template<class Engine>
class Gwesp : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double alpha;
    double oneexpa;
    double expa;
    std::vector< boost::container::flat_map<int,int> > sharedValues;
    int lastFrom;
    int lastTo;
public:

    Gwesp() : alpha(.5), oneexpa(1.0 - exp(-alpha)), expa(exp(alpha)),
    sharedValues(), lastFrom(0), lastTo(0){
    }

    virtual ~Gwesp(){};

    Gwesp(List params) : sharedValues(), lastFrom(0), lastTo(0){
        ParamParser p(name(), params);
        alpha = p.parseNext< double >("alpha");
        p.end();
        oneexpa = 1.0 - exp(-alpha);
        expa = exp(alpha);
    }

    std::string name(){
        return "gwesp";
    }

    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwesp."+a;
        std::vector<std::string> statnames(1,termname);
        return statnames;

    }


    //counts the number of shared neighbors between f and t.
    //in directed networks this only counts | t --> f --> neighbor --> t | cycles.
    int sharedNbrs(const BinaryNet<Engine>& net, int f, int t){
        if(!net.isDirected()){
            int tmp = f;
            f = std::min(f,t);
            t = std::max(tmp,t);
        }
        boost::container::flat_map<int,int>::iterator it = sharedValues[f].find(t);
        if(it != sharedValues[f].end()){
            return it->second;
        }
        NeighborIterator fit1, fend1, tit1, tend1;
        if(!net.isDirected()){
            fit1 = net.begin(f);
            fend1 = net.end(f);
            tit1 = net.begin(t);
            tend1 = net.end(t);
        }else{
            fit1 = net.inBegin(f);
            fend1 = net.inEnd(f);
            tit1 = net.outBegin(t);
            tend1 = net.outEnd(t);
        }
        int sn = 0;
        while(fit1 != fend1 && tit1 != tend1){
            if(*tit1 == *fit1){
                sn++;
                tit1++;
                fit1++;
            }else if(*tit1 < *fit1)
                tit1 = std::lower_bound(tit1,tend1,*fit1);
            else
                fit1 = std::lower_bound(fit1,fend1,*tit1);
        }
        return sn;
    }

    void setSharedValue(const BinaryNet<Engine>& net, int f, int t, int nbrs){
        if(!net.isDirected()){
            int tmp = f;
            f = std::min(f,t);
            t = std::max(tmp,t);
        }
        sharedValues[f][t] = nbrs;
    }
    void eraseSharedValue(const BinaryNet<Engine>& net, int f, int t){
        if(!net.isDirected()){
            int tmp = f;
            f = std::min(f,t);
            t = std::max(tmp,t);
        }
        sharedValues[f].erase(t);
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        this->init(1);
        double result = 0.0;
        sharedValues = std::vector< boost::container::flat_map<int,int> >();
        for(int i = 0 ; i<net.size();i++)
            sharedValues.push_back(boost::container::flat_map<int,int>());
        boost::shared_ptr<std::vector< std::pair<int,int> > > el = net.edgelist();
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            int sn = sharedNbrs(net, from, to);
            setSharedValue(net,from,to,sn);
            result += 1.0 - pow(oneexpa,sn);
        }
        this->stats[0] = expa * result;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(from);
            fend = net.end(from);
            tit = net.begin(to);
            tend = net.end(to);
        }else{
            fit = net.inBegin(from);
            fend = net.inEnd(from);
            tit = net.outBegin(to);
            tend = net.outEnd(to);
        }
        bool add = !net.hasEdge(from,to);
        double change = 2.0 * (add - 0.5);
        double delta = 0.0;
        int sn = 0;
        double mult = (1.0 - (!add ? 1.0/oneexpa : oneexpa));
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                sn++;
                //tie from to --> shared neighbor
                int tnsn = sharedNbrs(net,to,*tit);
                setSharedValue(net,to,*tit,tnsn + (add ? 1 : -1));
                delta += pow(oneexpa,tnsn) * mult;//pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change);

                //tie from shared neighbor --> from
                int nfsn = sharedNbrs(net,*tit,from);
                setSharedValue(net,*tit,from,nfsn + (add ? 1 : -1));
                delta += pow(oneexpa,nfsn) * mult;
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit = std::lower_bound(tit,tend,*fit);
            else
                fit = std::lower_bound(fit,fend,*tit);
        }
        if(add)
            setSharedValue(net,from,to,sn);
        else
            eraseSharedValue(net,from,to);
        lastFrom=from;
        lastTo=to;
        this->stats[0] += expa * (delta + change * (1.0 - pow(oneexpa,sn)));
    }

    void rollback(const BinaryNet<Engine>& net){
        BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
        pnet->toggle(lastFrom, lastTo);
        this->dyadUpdate(net, lastFrom, lastTo, std::vector<int>(), -1);
        pnet->toggle(lastFrom, lastTo);
    }

    bool isOrderIndependent(){
        return true;
    }
};

typedef Stat<Directed, Gwesp<Directed> > DirectedGwesp;
typedef Stat<Undirected, Gwesp<Undirected> > UndirectedGwesp;






/*!
 * geometrically weighted degree term - jc
 */

template<class Engine>
class GwDegree : public BaseStat< Engine > {
protected:
    double alpha;
    EdgeDirection direction; //same way of dealing with directionality as stars term. second param 1=IN , 2 =OUT
    double oneexpa;
    double expalpha;
public:

    GwDegree() : alpha(.5), direction(), oneexpa(0), expalpha(0){ //jc
    }

    virtual ~GwDegree(){};

    GwDegree(List params) : oneexpa(0), expalpha(0){
        ParamParser p(name(), params);
        alpha = p.parseNext< double >("alpha");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        p.end();
    }

    std::string name(){
        return "gwdegree";
    }

    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwdegree."+a;
        if(direction == IN)
            termname = "in-" + termname;
        if(direction == OUT)
            termname = "out-" + termname;
        std::vector<std::string> statnames(1,termname);
        return statnames;

    }

    virtual void calculate(const BinaryNet<Engine>& net){
        oneexpa = 1.0 - exp(-alpha);
        expalpha = exp(alpha);
        this->init(1);
        double result = 0.0;
        if (!net.isDirected()) {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.degree(i));
            }
        }
        else if (net.isDirected() && direction==IN) {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.indegree(i));
            }
        } else {
            for(int i=0;i<net.size();i++){
                result += 1.0 - pow(oneexpa,net.outdegree(i));
            }
        }


        this->stats[0] = expalpha * result;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        //we'll toggle the dyad betwen from and to
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5); //change in edge value (1, -1)
        double delta1 = 0.0;
        double delta2 = 0.0;

        if(!net.isDirected()){ //does checking every time slow down the toggling?
            delta1 = pow(oneexpa,net.degree(from)) - pow(oneexpa,net.degree(from)+change);
            delta2 = pow(oneexpa,net.degree(to)) - pow(oneexpa,net.degree(to)+change);
        }

        else if(net.isDirected() && direction==IN){
            delta1 = pow(oneexpa,net.indegree(to)) - pow(oneexpa,net.indegree(to)+change);
            delta2 = 0;
        }else {
            delta1 = pow(oneexpa,net.outdegree(from)) - pow(oneexpa,net.outdegree(from)+change);
            delta2 = 0;
        }

        this->stats[0] += expalpha*(delta1 + delta2);
    }

    bool isOrderIndependent(){
        return true;
    }
};

typedef Stat<Directed, GwDegree<Directed> > DirectedGwDegree;
typedef Stat<Undirected, GwDegree<Undirected> > UndirectedGwDegree;


template<class Engine>
class Gwdsp : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double alpha;
public:

    Gwdsp() : alpha(.5){
    }

    virtual ~Gwdsp(){};

    Gwdsp(List params){
        ParamParser p(name(), params);
        alpha = p.parseNext< double >("alpha");
        p.end();
    }

    std::string name(){
        return "gwdsp";
    }

    std::vector<std::string> statNames(){
        std::string a = asString(alpha);
        std::string termname = "gwdsp."+a;
        std::vector<std::string> statnames(1,termname);
        return statnames;

    }

    //counts the number of shared neighbors between f and t.
    //in directed networks this only counts | t --> f --> neighbor --> t | cycles.

    int sharedNbrs(const BinaryNet<Engine>& net, int f, int t){
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(f);
            fend = net.end(f);
            tit = net.begin(t);
            tend = net.end(t);
        }else{
            fit = net.inBegin(f);
            fend = net.inEnd(f);
            tit = net.outBegin(t);
            tend = net.outEnd(t);
        }

        int sn = 0;
        while(fit != fend && tit != tend){
            if(*tit == *fit){
                sn++;
                tit++;
                fit++;
            }else if(*tit < *fit)
                tit++;
            else
                fit++;
        }
        return sn;
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        this->init(1);

        //for each node, how many neighbors does its neighbor have? Where end index is greater than starting index, to avoid duplicates

        double result = 0.0;

        double oneexpa = 1 - exp(-alpha);
        int n = net.size();
        //std::vector<int> dp ;

        for (int f=0;f<n;f++){
            std::set<int> twoaways;
            NeighborIterator fit, fend;
            if(!net.isDirected()){
                fit = net.begin(f);
                fend = net.end(f);
            }else{
                fit = net.inBegin(f);
                fend = net.inEnd(f);
            }
            while(fit != fend){

                NeighborIterator tit, tend;

                if(!net.isDirected()){
                    tit = net.begin(*fit);
                    tend = net.end(*fit);
                }else{
                    tit = net.inBegin(*fit);
                    tend = net.inEnd(*fit);
                }
                while(tit != tend){
                    if(f < *tit){
                        //dp.push_back(sharedNbrs(net,f,*tit));
                        twoaways.insert(*tit);
                        tit++;
                    }else{
                        tit++;
                    }
                }
                fit++;
            }
            std::set<int>::iterator it;//set iterator
            for (it = twoaways.begin() ; it!=twoaways.end(); ++it)                 result += 1.0 - pow(oneexpa,sharedNbrs(net,f,*it));
        }

        this->stats[0] = exp(alpha) * result;
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double oneexpa = 1.0 - exp(-alpha);
        NeighborIterator fit, fend, tit, tend;
        if(!net.isDirected()){
            fit = net.begin(from);
            fend = net.end(from);
            tit = net.begin(to);
            tend = net.end(to);
        }else{
            fit = net.inBegin(from);
            fend = net.inEnd(from);
            tit = net.outBegin(to);
            tend = net.outEnd(to);
        }
        //double add = !net.hasEdge(from,to);
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        double delta = 0.0;
        //Would this be faster if I used the edgelist as in gwesp? -jc
        //boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();

        while(fit != fend){
            if (*fit != to){ //iterate over the neighbors of from, except for to
                int tnsn = sharedNbrs(net,*fit,to);
                delta += pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change);
            }
            fit++;
        }
        while(tit != tend){
            if(*tit != from){
                int tnsn = sharedNbrs(net,from,*tit);
                delta += (pow(oneexpa,tnsn) - pow(oneexpa,tnsn + change));
            }
            tit++;
        }
        this->stats[0] += (exp(alpha)*delta);
    }

    bool isOrderIndependent(){
        return true;
    }

};

typedef Stat<Directed, Gwdsp<Directed> > DirectedGwdsp;
typedef Stat<Undirected, Gwdsp<Undirected> > UndirectedGwdsp;

//Edgewise Shared Parnters. One stat for each user-generated value.
template<class Engine>
class Esp : public BaseStat< Engine > {
    //first part of code based on the code for Degree
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    std::vector<int> esps;
    int type;

public:

    Esp(){
        type = 2;
    }

    virtual ~Esp(){}; //stick with vcalculate?

    Esp(std::vector<int> esps1){
        esps = esps1;
        type = 2;
    }


    /*!
     * \param params 	a list of length 1, the first element of which is an integer vector of edgewise shared partners
     */

    Esp(List params){
        ParamParser p(name(), params);
        esps = p.parseNext< std::vector<int> >("d");
        type = p.parseNext("type", 2);
        if(type < 1 || type > 4)
            ::Rf_error("ESP: type must be 1,2,3, or 4");
        p.end();
    }

    std::string name(){
        return "esp";
    }


    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<esps.size();i++){
            int e = esps[i];
            std::string nm = "esp."+asString(e);
            if(type != 2)
                nm = asString(type) + "-" + nm;
            statnames.push_back(nm);
        }
        return statnames;
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        int nstats = esps.size();
        this->init(nstats);

        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();

        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            int espi = sharedNbrs(net, from, to, type);
            for(int j=0;j<nstats;j++){
                this->stats[j] += espi==esps[j];
            }
        }

    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        int nstats = esps.size();
        int espi = sharedNbrs(net, from, to, type);
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        for(int j=0;j<nstats;j++){ //main edge change from-to
            this->stats[j] += change*(espi==esps[j]);
        }
        if(type == 1 || !net.isDirected()){
            NeighborIterator fit, fend, tit, tend;
            if(!net.isDirected()){
                fit = net.begin(from);
                fend = net.end(from);
                tit = net.begin(to);
                tend = net.end(to);
            }else{
                fit = net.inBegin(from);
                fend = net.inEnd(from);
                tit = net.outBegin(to);
                tend = net.outEnd(to);
            }
            while(fit != fend && tit != tend){
                if(*tit == *fit){ //it's a shared neighbor
                    int fnsn = sharedNbrs(net, *fit, from,  type);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += (fnsn+change)==esps[j];
                        this->stats[j] -= fnsn==esps[j];
                    }
                    int tnsn = sharedNbrs(net, to, *fit, type);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += (tnsn+change)==esps[j];
                        this->stats[j] -= tnsn==esps[j];
                    }
                    tit++;
                    fit++;
                }else if(*tit < *fit)
                    tit++;
                else
                    fit++;
            }
        }else{
            // A bit brute force. Will work with any definition of shared nbr

            BinaryNet<Engine>* pnet = const_cast< BinaryNet<Engine>* > (&net);
            std::set<int> nbrs;


            nbrs.insert(net.inBegin(to), net.inEnd(to));
            for(std::set<int>::iterator it = nbrs.begin(); it != nbrs.end(); it++){
                if(net.hasEdge(from, *it) || net.hasEdge(*it, from)){
                    int curShared = sharedNbrs(net, *it, to,  type);
                    pnet->toggle(from, to);
                    int newShared = sharedNbrs(net, *it, to,  type);
                    pnet->toggle(from, to);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += newShared==esps[j];
                        this->stats[j] -= curShared==esps[j];
                    }
                }
            }
            nbrs.clear();

            nbrs.insert(net.inBegin(from), net.inEnd(from));
            for(std::set<int>::iterator it = nbrs.begin(); it != nbrs.end(); it++){
                if(net.hasEdge(to, *it) || net.hasEdge(*it, to)){
                    int curShared = sharedNbrs(net, *it, from,  type);
                    pnet->toggle(from, to);
                    int newShared = sharedNbrs(net, *it, from,  type);
                    pnet->toggle(from, to);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += newShared==esps[j];
                        this->stats[j] -= curShared==esps[j];
                    }
                }
            }
            nbrs.clear();

            nbrs.insert(net.outBegin(to), net.outEnd(to));
            for(std::set<int>::iterator it = nbrs.begin(); it != nbrs.end(); it++){
                if(net.hasEdge(from, *it) || net.hasEdge(*it, from)){
                    int curShared = sharedNbrs(net, to, *it,  type);
                    pnet->toggle(from, to);
                    int newShared = sharedNbrs(net, to, *it, type);
                    pnet->toggle(from, to);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += newShared==esps[j];
                        this->stats[j] -= curShared==esps[j];
                    }
                }
            }
            nbrs.clear();

            nbrs.insert(net.outBegin(from), net.outEnd(from));
            for(std::set<int>::iterator it = nbrs.begin(); it != nbrs.end(); it++){
                if(net.hasEdge(to, *it) || net.hasEdge(*it, to)){
                    int curShared = sharedNbrs(net, from, *it,  type);
                    pnet->toggle(from, to);
                    int newShared = sharedNbrs(net, from, *it,  type);
                    pnet->toggle(from, to);
                    for(int j=0;j<nstats;j++){ // side edge change +/-1
                        this->stats[j] += newShared==esps[j];
                        this->stats[j] -= curShared==esps[j];
                    }
                }
            }
            nbrs.clear();

        }

    }

    bool isOrderIndependent(){
        return true;
    }

};

typedef Stat<Directed, Esp<Directed> > DirectedEsp;
typedef Stat<Undirected, Esp<Undirected> > UndirectedEsp;

/*!
 * Great circle distance between two long-lat points
 */
template<class Engine>
class GeoDist : public BaseStat< Engine > {
protected:
    std::string latVarName;
    int latIndex;
    std::string longVarName;
    int longIndex;
    std::vector<double> distCuts;
public:

    GeoDist() : latIndex(-1), longIndex(-1){}

    virtual ~GeoDist(){};

    GeoDist(List params) : latIndex(-1), longIndex(-1){
        ParamParser p(name(), params);
        longVarName = p.parseNext< std::string >("long");
        latVarName = p.parseNext< std::string >("lat");
        distCuts = p.parseNext("distCuts", std::vector<double>(1, 41000.0));
        p.end();
    }

    std::string name(){
        return "geoDist";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<distCuts.size();i++){
            int e = distCuts[i];
            std::string nm = "geoDist."+asString(e);
            statnames.push_back(nm);
        }
        return statnames;
    }

    double dist(double th1, double ph1, double th2, double ph2)
    {
        double dx, dy, dz;
        ph1 -= ph2;
        ph1 *= (3.1415926536 / 180.0), th1 *= (3.1415926536 / 180.0), th2 *= (3.1415926536 / 180.0);

        dz = sin(th1) - sin(th2);
        dx = cos(ph1) * cos(th1) - cos(th2);
        dy = sin(ph1) * cos(th1);
        return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * 6371.0;
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        std::vector<std::string> vars = net.continVarNames();
        //int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == longVarName){
                longIndex = i;
            }
            if(vars[i] == latVarName){
                latIndex = i;
            }
        }
        if(latIndex<0)
            ::Rf_error("latitude attribute not found in network");
        for(int i=0;i<net.size();i++){
            double deg = net.continVariableValue(latIndex,i);
            if(deg<-90 || deg>90)
                Rf_error("Latitude values out of range.");
        }

        if(longIndex<0)
            ::Rf_error("longitude attribute not found in network");
        for(int i=0;i<net.size();i++){
            double deg = net.continVariableValue(longIndex,i);
            if(deg<-180 || deg>180)
                Rf_error("Longitude values out of range.");
        }

        int nstats = distCuts.size();
        this->init(nstats);

        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            double distance = dist(
                    net.continVariableValue(latIndex,from),
                    net.continVariableValue(longIndex,from),
                    net.continVariableValue(latIndex,to),
                    net.continVariableValue(longIndex,to)
            );
            for(int j=0;j<distCuts.size();j++){
                this->stats[j] += std::min(distCuts[j], distance);
            }
        }
        //this->stats[0] = result / (double) net.nEdges();
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        double distance = dist(
                net.continVariableValue(latIndex,from),
                net.continVariableValue(longIndex,from),
                net.continVariableValue(latIndex,to),
                net.continVariableValue(longIndex,to)
        );
        for(int j=0;j<distCuts.size();j++){
            this->stats[j] += change * std::min(distCuts[j], distance);
        }
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};

typedef Stat<Directed, GeoDist<Directed> > DirectedGeoDist;
typedef Stat<Undirected, GeoDist<Undirected> > UndirectedGeoDist;


/*!
 * distance measure
 */
template<class Engine>
class Dist : public BaseStat< Engine > {
protected:
    EdgeDirection direction;
    std::vector< std::string > varNames;
    std::vector<int> indices;
public:

    Dist() : direction(UNDIRECTED){}

    virtual ~Dist(){};

    Dist(List params){
        ParamParser p(name(), params);
        varNames = p.parseNext< std::vector<std::string> >("varNames");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        p.end();
    }

    std::string name(){
        return "dist";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"dist");
        return statnames;
    }

    double dist(const BinaryNet<Engine>& net, int from, int to){
        double ssq = 0.0;
        for(int j=0;j<indices.size();j++){
            ssq += pow(net.continVariableValue(indices[j],from) -
                    net.continVariableValue(indices[j],to), 2.0);
        }
        return sqrt(ssq);
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        std::vector<std::string> vars = net.continVarNames();
        //int variableIndex = -1;
        indices = std::vector<int>(varNames.size(),-1);
        for(int i=0;i<vars.size();i++){
            for(int j=0;j<varNames.size();j++){
                if(vars[i] == varNames[j]){
                    indices[j] = i;
                }
            }
        }
        for(int i=0;i<varNames.size();i++)
            if(indices[i] < 0)
                ::Rf_error("dist: variable not found in network");

        int nstats = 1;
        this->init(nstats);


        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
        double result = 0.0;
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            result += dist(net, from,to);
        }
        this->stats[0] = result;
        //this->stats[0] = result / (double) net.nEdges();
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        this->stats[0] = this->stats[0] + change * dist(net, from,to);
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};

typedef Stat<Directed, Dist<Directed> > DirectedDist;
typedef Stat<Undirected, Dist<Undirected> > UndirectedDist;



/*!
 * distance measure
 */
template<class Engine>
class AbsDiff : public BaseStat< Engine > {
protected:
    std::vector< std::string > varNames;
    std::vector<int> indices;
    double power;
public:

    AbsDiff() : power(1.0){}

    virtual ~AbsDiff(){};

    AbsDiff(List params){
        ParamParser p(name(), params);
        varNames = p.parseNext< std::vector<std::string> >("varNames");
        power = p.parseNext("power", 1.0);
        p.end();
    }

    std::string name(){
        return "absDiff";
    }

    std::vector<std::string> statNames(){
        std::string nm = "absDiff";
        for(int i=0;i<varNames.size();i++)
            nm = nm + "." + varNames.at(i);
        std::vector<std::string> statnames(1, nm);
        return statnames;
    }

    double dist(const BinaryNet<Engine>& net, int from, int to){
        double ssq = 0.0;
        for(int j=0;j<indices.size();j++){
            ssq += pow(abs(net.continVariableValue(indices[j],from) -
                    net.continVariableValue(indices[j],to)), power);
        }
        return ssq;
    }

    virtual void calculate(const BinaryNet<Engine>& net){
        std::vector<std::string> vars = net.continVarNames();
        //int variableIndex = -1;
        indices = std::vector<int>(varNames.size(),-1);
        for(int i=0;i<vars.size();i++){
            for(int j=0;j<varNames.size();j++){
                if(vars[i] == varNames[j]){
                    indices[j] = i;
                }
            }
        }
        for(int i=0;i<varNames.size();i++)
            if(indices[i] < 0)
                ::Rf_error("dist: variable not found in network");

        int nstats = 1;
        this->init(nstats);


        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();
        double result = 0.0;
        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            result += dist(net, from,to);
        }
        this->stats[0] = result;
        //this->stats[0] = result / (double) net.nEdges();
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        this->stats[0] = this->stats[0] + change * dist(net, from, to);
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};

typedef Stat<Directed, AbsDiff<Directed> > DirectedAbsDiff;
typedef Stat<Undirected, AbsDiff<Undirected> > UndirectedAbsDiff;



/*!
 * barabasi-albert type change statistic. Order dependent
 */
template<class Engine>
class PreferentialAttachment : public BaseStat<Engine>{
    EdgeDirection direction;
    double k;
public:
    PreferentialAttachment(){
        k = 1.0;
        direction = IN;
    }

    /*!
     * constructor. params is unused
     */
    PreferentialAttachment(List params){
        ParamParser p(name(), params);
        k = p.parseNext("k", 1.0);
        direction = p.parseNextDirection("direction", IN);
        p.end();
    }


    std::string name(){
        return "preferentialAttachment";
    }

    double degree(const BinaryNet<Engine> & net, int i){
        if(net.isDirected()){
            if(direction == IN)
                return net.indegree(i);
            else if(direction == OUT)
                return net.outdegree(i);
            else
                return net.indegree(i) + net.outdegree(i);
        }else
            return net.degree(i);
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"preferentialAttachment");
        return statnames;
    }

    void calculate(const BinaryNet<Engine>& net){
        this->init(1);
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        bool hasEdge = net.hasEdge(from,to);
        double direction = hasEdge ? -1.0 : 1.0;
        double totDegree = (net.nEdges() - hasEdge) * 2.0;
        int alter = order[actorIndex] == from ? to : from;
        double deg = degree(net, alter) - hasEdge;
        double netSize = actorIndex + 1.0;

        this->stats[0] += direction * log( (k + deg) / (k*netSize + totDegree));

        //BaseOffset<Engine>::update(direction * log( (k + net.degree(alter) - hasEdge) / (k*(actorIndex + 1.0) + (net.nEdges() - hasEdge) * 2.0)), 0);
    }

};

typedef Stat<Directed, PreferentialAttachment<Directed> > DirectedPreferentialAttachment;
typedef Stat<Undirected, PreferentialAttachment<Undirected> > UndirectedPreferentialAttachment;


template<class Engine>
class SharedNbrs : public BaseStat< Engine > {
protected:
    typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
    double k;
public:


    SharedNbrs(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas = t;
        k=1.0;
    }
    SharedNbrs(List params){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats=v;
        this->thetas = t;

        ParamParser p(name(), params);
        k = p.parseNext("k", 1.0);
        p.end();
    }

    std::string name(){
        return "sharedNbrs";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"sharedNbrs");
        return statnames;
    }




    void calculate(const BinaryNet<Engine>& net){
        this->init(1);
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,
            const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double shared = sharedNbrs(net, from, to);
        bool hasEdge = net.hasEdge(from,to);
        int deg = net.degree(order[actorIndex]) - hasEdge;
        int alter = order[actorIndex] == from ? to : from;
        double altDeg = net.degree(alter) - hasEdge;
        double netSize = actorIndex + 1.0;
        double totDegree = (net.nEdges() - hasEdge) * 2.0;
        //if(deg < 0.5){
        //	shared = altDeg;
        //	deg = 1.0;
        //}
        //double value = log((k + shared) / (netSize*k + deg*totDegree));
        double sdegs = 0.0;
        NeighborIterator fit = net.begin(order[actorIndex]);
        NeighborIterator fend = net.end(order[actorIndex]);
        //while(fit != fend){
        //	sdegs += std::min(net.degree(*fit), deg + !hasEdge) - 1.0;
        //	fit++;
        //}

        double d = std::min((double)altDeg, (double)deg);
        if(d < .5) d++;

        double value = log(k + shared / d);

        //double maxShared = std::min(altDeg, (double) deg) + !hasEdge - 1.0;
        //double value = log((k + shared) / (k + maxShared));
        //double value = log((k + shared) / (k + 0.5 * sdegs));
        //if(sdegs < .5) sdegs++;
        //double value = (shared) / (0.5 * sdegs);
        if(hasEdge){
            BaseOffset<Engine>::update(-value, 0);
        }else{
            BaseOffset<Engine>::update(value, 0);
        }
    }

};

typedef Stat<Directed, SharedNbrs<Directed> > DirectedSharedNbrs;
typedef Stat<Undirected, SharedNbrs<Undirected> > UndirectedSharedNbrs;




/*!
 */
template<class Engine>
class NodeLogMaxCov : public BaseStat< Engine > {
protected:
    EdgeDirection direction;
    std::string variableName;
    int varIndex;
    bool isDiscrete;
public:

    NodeLogMaxCov(){
        varIndex =  0;
        direction = UNDIRECTED;
        isDiscrete = false;
    }

    NodeLogMaxCov(std::string name,EdgeDirection d){
        varIndex = 0;
        direction = d;
        variableName = name;
        isDiscrete = false;
    }

    NodeLogMaxCov(std::string name){
        varIndex = 0;
        direction = UNDIRECTED;
        variableName = name;
        isDiscrete = false;
    }


    NodeLogMaxCov(List params){
        varIndex = 0;
        isDiscrete=false;

        ParamParser p(name(), params);
        variableName = p.parseNext< std::string >("name");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        p.end();

    }

    std::string name(){
        return "nodeLogMaxCov";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        statnames.assign(1,"nodeLogMaxCov."+variableName);
        return statnames;
    }

    double getValue(const BinaryNet<Engine>& net, int ind){
        double val;
        if(isDiscrete)
            val = net.discreteVariableValue(varIndex,ind);
        else
            val = net.continVariableValue(varIndex,ind);
        return val;
    }

    void calculate(const BinaryNet<Engine>& net){
        isDiscrete = false;
        std::vector<std::string> vars = net.continVarNames();
        int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName){
                variableIndex = i;
            }
        }
        if(variableIndex == -1){
            isDiscrete = true;
            vars = net.discreteVarNames();
            for(int i=0;i<vars.size();i++){
                if(vars[i] == variableName){
                    variableIndex = i;
                }
            }
        }
        if(variableIndex<0)
            ::Rf_error("nodal attribute not found in network");
        varIndex = variableIndex;
        int nstats = 1;
        this->init(nstats);
        boost::shared_ptr< std::vector<std::pair<int,int> > > el = net.edgelist();

        for(int i=0;i<el->size();i++){
            int from = el->at(i).first;
            int to = el->at(i).second;
            double val1 = getValue(net,from);
            double val2 = getValue(net,to);
            double val = val1 > val2 ? val1 : val2;
            for(int j=0;j<nstats;j++){
                this->stats[j] += log(val);
            }
        }
    }

    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
        double val1 = getValue(net,from);
        double val2 = getValue(net,to);
        double val = val1 > val2 ? val1 : val2;
        this->stats[0] += change * log(val);
    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }

};

typedef Stat<Directed, NodeLogMaxCov<Directed> > DirectedNodeLogMaxCov;
typedef Stat<Undirected, NodeLogMaxCov<Undirected> > UndirectedNodeLogMaxCov;



/*!
 * Differential activity by group
 *
 */
template<class Engine>
class NodeFactor : public BaseStat< Engine > {
protected:
    EdgeDirection direction;
    std::string variableName;
    int varIndex;
    int nstats;
public:

    NodeFactor(){
        varIndex = nstats = 0;
        direction = UNDIRECTED;
    }

    NodeFactor(std::string name,EdgeDirection d){
        varIndex = nstats = 0;
        direction = d;
        variableName = name;
    }

    NodeFactor(std::string name){
        varIndex = nstats = 0;
        direction = UNDIRECTED;
        variableName = name;
    }


    NodeFactor(List params){
        varIndex = nstats = 0;

        ParamParser p(name(), params);
        variableName = p.parseNext< std::string >("name");
        direction = p.parseNextDirection("direction", UNDIRECTED);
        p.end();
    }

    std::string name(){
        return "nodeFactor";
    }

    std::vector<std::string> statNames(){
        std::vector<std::string> statnames;
        for(int i=0;i<nstats;i++){
            std::string nm = "nodeFactor."+variableName+"."+asString(i+1);
            if(direction == IN)
                nm = "in-" + nm;
            if(direction == OUT)
                nm = "out-" + nm;
            statnames.push_back(nm);
        }
        return statnames;
    }


    inline int degree(const BinaryNet<Engine>& net,int i){
        int result = 0;
        if(net.isDirected()){
            if(direction==OUT || direction==UNDIRECTED)
                result += net.outdegree(i);
            if(direction==IN || direction==UNDIRECTED)
                result += net.indegree(i);
        }else
            result = net.degree(i);
        return result;
    }

    void calculate(const BinaryNet<Engine>& net){
        std::vector<std::string> vars = net.discreteVarNames();
        int variableIndex = -1;
        for(int i=0;i<vars.size();i++){
            if(vars[i] == variableName){
                variableIndex = i;
            }
        }
        if(variableIndex<0)
            ::Rf_error("nodal attribute not found in network");
        varIndex = variableIndex;
        int nlevels = net.discreteVariableAttributes(variableIndex).labels().size();
        nstats = nlevels-1;
        this->init(nstats);
        double n = net.size();
        double deg = 0.0;
        for(int i=0;i<n;i++){
            deg = degree(net,i);
            int val = net.discreteVariableValue(varIndex,i) - 1;
            if(val<nstats)
                this->stats[val] += deg;
        }
    }


    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        int fromVal = net.discreteVariableValue(varIndex,from)-1;
        int toVal = net.discreteVariableValue(varIndex,to)-1;
        int change;
        change = !net.hasEdge(from,to)? 1 : -1;

        if( (direction==UNDIRECTED || direction==OUT) && fromVal<nstats)
            this->stats[fromVal] += change;
        if( (direction==UNDIRECTED || direction==IN) && toVal<nstats)
            this->stats[toVal] += change;

    }

    bool isOrderIndependent(){
        return true;
    }

    bool isDyadIndependent(){
        return true;
    }
};

typedef Stat<Directed, NodeFactor<Directed> > DirectedNodeFactor;
typedef Stat<Undirected, NodeFactor<Undirected> > UndirectedNodeFactor;

/**
* An example lolog statistic, defined as the sum of dcov over the values that have edges
*/
template<class Engine>
class EdgeCov : public BaseStat< Engine > {
protected:
    NumericMatrix dcov; //the dyadic matrix
    std::string termName; /*!< the name of the term variable */
public:
    
    //Constructor
    EdgeCov(){}
    
    //Parse parameters
    EdgeCov(List params){
      ParamParser p(name(), params);
      dcov = p.parseNext< NumericMatrix >("x");
      termName = p.parseNext< std::string >("name","");
      p.end();
    }
    
    //The name 
    std::string name(){return "edgeCov";}
    
    std::vector<std::string> statNames(){
      std::vector<std::string> statnames(1,"edgeCov."+termName);
      return statnames;
    }
    
    //Calculate the statistic
    virtual void calculate(const BinaryNet<Engine>& net){
        if((dcov.nrow() != net.size()) || (dcov.ncol() != net.size())){
          ::Rf_error("EdgeCov error: the dyadic covariate matrix should have the same dimensions as the adjacency matrix.");
        }
        std::vector<double> v(1,0);
        this->stats=v;
        this->lastStats = std::vector<double>(1,0.0);
        if(this->thetas.size()!=1)
            this->thetas = v;
        if(net.isDirected()){
         for(int i=0;i<net.size();i++){
          for(int j=0;j<net.size();j++){
            this->stats[0] += net.hasEdge(i,j)*dcov(i,j);
          }
         }
        }else{
         for(int i=1;i<net.size();i++){
          for(int j=0;j<i;j++){
            this->stats[0] += net.hasEdge(i,j)*dcov(i,j);
          }
         }
        }
    }
    
    //Update the statistic given a dyad toggle
    virtual void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
      BaseOffset<Engine>::resetLastStats();
      bool addingEdge = !net.hasEdge(from,to);
      double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
      this->stats[0] += change * dcov(from,to);
    }
    
    //Declare that this statistic is order independent
    bool isOrderIndependent(){
      return true;
    }
    
    //Declare that this statistic is dyad independent
    bool isDyadIndependent(){
      return true;
    }
    
};

typedef Stat<Undirected, EdgeCov<Undirected> > UndirectedEdgeCov;
typedef Stat<Directed, EdgeCov<Directed> > DirectedEdgeCov;


/*!
 * twoPath
 */
template<class Engine>
class TwoPath : public BaseStat< Engine > {
public:
    TwoPath(){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats = v;
        this->thetas = t;
    }
    
    /*!
    * constructor. params is unused
    */
    TwoPath(List params){
        std::vector<double> v(1,0.0);
        std::vector<double> t(1,0.0);
        this->stats = v;
        this->thetas = t;
    }
    
    std::string name(){
        return "twoPath";
    }
    std::vector<std::string> statNames(){
        std::vector<std::string> statnames(1,"twoPath");
        return statnames;
    }
    
    void calculate(const BinaryNet<Engine>& net){
        this->init(1);
        double rec = 0.0;
        int from, to;
        boost::shared_ptr< std::vector< std::pair<int,int> > > edges = net.edgelist();
        if(!net.isDirected()){
            for(int i=0; i<net.size();i++){
                double nEd = net.degree(i);
                rec += nchoosek(nEd, 2.0);
            }
        }else{
            for(int i=0;i<edges->size();i++){
                from = (*edges)[i].first;
                to = (*edges)[i].second;
                rec+=(net.outdegree(to)-net.hasEdge(to,from));
            }
        }
        std::vector<double> v(1,rec);
        this->stats=v;
    }
    
    void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
        BaseOffset<Engine>::resetLastStats();
        bool removingEdge = net.hasEdge(from,to);
        bool hasReverse = net.hasEdge(to,from);
        double change;
        if(!net.isDirected()){
            change = net.degree(from) + net.degree(to) - 2.0 * removingEdge;
        }else{
            change = net.indegree(from) + net.outdegree(to) - 2.0*hasReverse;
        }
        if(removingEdge) change = -change;
        BaseOffset<Engine>::update(change,0);
    }
    
    bool isOrderIndependent(){
        return true;
    }
  
};

typedef Stat<Directed, TwoPath<Directed> > DirectedTwoPath;
typedef Stat<Undirected, TwoPath<Undirected> > UndirectedTwoPath;



template<class Engine>
class EdgeCovSparse : public BaseStat< Engine > {
protected:
  boost::unordered::unordered_map<std::pair<int, int>, double> map;
  std::string termName; /*!< the name of the term variable */
    
  void convertMatrix(SEXP x){
    Rcpp::Language call("getNamespace","Matrix");
    //call.eval();
    //Environment Matrix("package:Matrix"); 
    Environment Matrix = call.eval();
    Function summary = Matrix["summary"];
    //Rcpp::Function summary("summary");
    //SEXP sumry = summary(x);
    //sumry.attr("class") = "data.frame"
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(summary(x));
    Rcpp::IntegerVector ivals = df["i"];
    Rcpp::IntegerVector jvals = df["j"];
    Rcpp::NumericVector values = df["x"];
    for(int i=0; i < df.nrow(); i++){
      map[std::make_pair(ivals[i] - 1,jvals[i] - 1)] = values[i];
    }
  }
  
  double dcov(int i, int j, bool directed){
    std::pair<int,int> p = std::make_pair(i,j);
    double v = 0.0;
    //Rcpp::Rcout << "\n\n\n" << map.count(std::make_pair(25,1));
    if(map.count(p))
      v = map.at(p);
    else if(!directed){
      p = std::make_pair(j,i);
      if(map.count(p))
        v = map.at(p);
    }
    return v;
  }
public:
  
  //Constructor
  EdgeCovSparse(){}
  
  //Parse parameters
  EdgeCovSparse(List params){
    ParamParser p(name(), params);
    Rcpp::RObject x = p.parseNext< Rcpp::RObject >("x");
    convertMatrix(x);
    termName = p.parseNext< std::string >("name","");
    p.end();
  }
  
  //The name 
  std::string name(){return "edgeCovSparse";}
  
  std::vector<std::string> statNames(){
    std::vector<std::string> statnames(1,"edgeCovSparse."+termName);
    return statnames;
  }
  
  //Calculate the statistic
  virtual void calculate(const BinaryNet<Engine>& net){
    //if(dcov.n_rows != net.size() | dcov.n_cols != net.size()){
    //  ::Rf_error("EdgeCov error: the dyadic covariate matrix should have the same dimensions as the adjacency matrix.");
    //}
    std::vector<double> v(1,0);
    this->stats=v;
    this->lastStats = std::vector<double>(1,0.0);
    if(this->thetas.size()!=1)
      this->thetas = v;
    if(net.isDirected()){
      for(int i=0;i<net.size();i++){
        for(int j=0;j<net.size();j++){
          this->stats[0] += net.hasEdge(i,j)*dcov(i,j, true);
        }
      }
    }else{
      for(int i=1;i<net.size();i++){
        for(int j=0;j<i;j++){
          this->stats[0] += net.hasEdge(i,j)*dcov(i,j, false);
        }
      }
    }
  }
  
  //Update the statistic given a dyad toggle
  virtual void dyadUpdate(const BinaryNet<Engine>& net,const int &from,const int &to,const std::vector<int> &order,const int &actorIndex){
    BaseOffset<Engine>::resetLastStats();
    bool addingEdge = !net.hasEdge(from,to);
    double change = 2.0 * (!net.hasEdge(from,to) - 0.5);
    this->stats[0] += change * dcov(from, to, net.isDirected());
  }
  
  //Declare that this statistic is order independent
  bool isOrderIndependent(){
    return true;
  }
  
  //Declare that this statistic is dyad independent
  bool isDyadIndependent(){
    return true;
  }
  
};

typedef Stat<Undirected, EdgeCovSparse<Undirected> > UndirectedEdgeCovSparse;
typedef Stat<Directed, EdgeCovSparse<Directed> > DirectedEdgeCovSparse;



#include <Rcpp.h>


}


#endif /* STATS_H_ */
