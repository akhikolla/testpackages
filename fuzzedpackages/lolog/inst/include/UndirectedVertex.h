#ifndef UNDIRECTEDVERTEXH_
#define UNDIRECTEDVERTEXH_
#include "Rcpp.h"
#include "Vertex.h"
#include <iostream>
#include <assert.h>
namespace lolog {

/*!
 *
 */
class UndirectedVertex: public Vertex {

protected:
    Set edgs;		        //!a set of edges

    Set miss;		        //!a set of missing dyads
    Set obs;
    bool useMissingSet;     //!should missing or observed be used to keep track of the
    //!missing dyads
    int nverts;	            //! the number of vertices in the network

    void refreshMissingRepresentation(){
        bool um = useMissingSet;
        if(um && miss.size() > 0.6*nverts){
            obs = Set();
            Set::iterator it = miss.begin();
            Set::iterator end = miss.end();
            Set::iterator lastInsertedLoc = obs.begin();
            for(int i=0;i<nverts;i++){
                if(i==this->idNum)
                    continue;
                if(it!=end && i==*it){
                    it++;
                    continue;
                }
                //std::cout << i<<" ";
                lastInsertedLoc = obs.insert(lastInsertedLoc,i);
            }
            useMissingSet = false;
            miss = Set();
            //oobserved.insert(this->idNum);
            //std::cout<<"to obs ";
        }else if(!um && obs.size() > 0.6*nverts){
            miss = Set();
            Set::iterator it = obs.begin();
            Set::iterator end = obs.end();
            Set::iterator lastInsertedLoc = miss.begin();
            for(int i=0;i<nverts;i++){
                if(it!=end && i==*it){
                    it++;
                    continue;
                }
                //std::cout << i<<" ";
                lastInsertedLoc = miss.insert(lastInsertedLoc,i);
            }
            useMissingSet = true;
            obs = Set();
            //std::cout<<"to miss ";
            miss.erase(this->idNum);
        }
    }

public:
    UndirectedVertex(int numVerts){
        nverts = numVerts;
        useMissingSet = true;
    }
    virtual ~UndirectedVertex(){

    }


    bool addEdge(int from){
        return edgs.insert(from).second;
    }

    bool hasEdge(int from){
        Set::iterator it = edgs.find(from);
        return it!=edgs.end();
    }

    int degree(){ return edgs.size();}

    bool removeEdge(int from){
        return edgs.erase(from)==1;
    }

    const Set& edges() const{
        return edgs;
    }

    void clearEdges(){
        edgs.clear();
    }

    int networkSize() const{return nverts;}
    void setNetworkSize(int netSize){nverts = netSize;}

    bool isMissing(int to){
        if(to == this->idNum)
            return false;
        if(useMissingSet){
            Set::iterator it = miss.find(to);
            return it!=miss.end();
        }else{
            Set::iterator it = obs.find(to);
            return it==obs.end();
        }
    }

    /*!
     * \returns true if was missing
     */
    bool setMissing(int to,bool value){
        bool ret;
        if(value){
            if(useMissingSet)
                ret = !miss.insert(to).second;
            else
                ret = obs.erase(to)==0;
        }else{
            if(!useMissingSet)
                ret = obs.insert(to).second;
            else
                ret = miss.erase(to)!=0;
        }
        refreshMissingRepresentation();
        return ret;
    }

    void setAllMissing(){
        useMissingSet=false;
        miss=Set();
        obs=Set();
    }

    void setAllObserved(){
        useMissingSet=true;
        miss=Set();
        obs=Set();
    }

    Set missing() const{
        if(useMissingSet){
            Set tmp = miss;
            tmp.erase(this->idNum);
            return tmp;
        }else{
            Set tmp = Set();
            Set::const_iterator it = obs.begin();
            for(int i=0;i<nverts;i++){
                if(it!=obs.end() && i==*it){
                    it++;
                    continue;
                }
                tmp.insert(tmp.end(),i);
            }
            return tmp;
        }
    }

    int nMissing(){
        if(useMissingSet)
            return miss.size();
        else
            return nverts - 1 - obs.size();
    }

    int randomMissingDyad(){
        assert(nMissing()>0);
        double percMissing = nMissing()/(nverts - 1.0);
        if(percMissing>0.05){
            for(int i=0;i<15;i++){
                int nbr = floor(Rf_runif(0,nverts - 1.0));
                if(nbr>=this->idNum)
                    nbr++;
                if(isMissing(nbr))
                    return nbr;
            }
        }
        int index = floor(Rf_runif(0,(double)nMissing()));
        if(useMissingSet){
            Set::iterator it = miss.begin();
            for(int i=0;i<index;i++)
                it++;
            return *it;
        }else{
            Set::iterator it = obs.begin();
            for(;it!=obs.end();it++){
                if(*it>index && index!=this->idNum)
                    return index;
                index++;
            }
            return index;
        }
        ::Rf_error("randomMissingDyad: logic error");
        return -1;
    }
};

} /* namespace lolog */
#endif /* UNDIRECTEDVERTEXH_ */
