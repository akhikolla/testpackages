#ifndef DVERTEXH_
#define DVERTEXH_
#include "Rcpp.h"
#include <vector>
#include <set>
#include "Vertex.h"
#include <iostream>
#include <assert.h>

//note randomNonEdge requires std::set
//#define Set std::set<int>
namespace lolog {


/*!
 * A directed vertex, for use in a DirectedNet.
 *
 * Has a sparse representation of edges (a binary tree) and a sparse representation
 * of missingness. If all or no out-dyads are missing, the missingness representation
 * takes up no additional space.
 */
class DirectedVertex : public Vertex {
protected:
    Set iedges;		//!a set of in edges
    Set oedges;		//!a set of out edges

    Set omissing;	//!a set of missing out dyads
    Set oobserved;	//!a set of observed dyads
    bool useMissingSet; 	//!should omissing or oobserved be used to keep track of the
    //!missing dyads
    int nverts;	//! the number of vertices in the network

    void refreshMissingRepresentation(){
        bool um = useMissingSet;
        if(um && omissing.size() > 0.6*nverts){
            oobserved = Set();
            Set::iterator it = omissing.begin();
            Set::iterator end = omissing.end();
            Set::iterator lastInsertedLoc = oobserved.begin();
            for(int i=0;i<nverts;i++){
                if(i==this->idNum)
                    continue;
                if(it!=end && i==*it){
                    it++;
                    continue;
                }
                //std::cout << i<<" ";
                lastInsertedLoc = oobserved.insert(lastInsertedLoc,i);
            }
            useMissingSet = false;
            omissing = Set();
            //oobserved.insert(this->idNum);
            //std::cout<<"to obs ";
        }else if(!um && oobserved.size() > 0.6*nverts){
            omissing = Set();
            Set::iterator it = oobserved.begin();
            Set::iterator end = oobserved.end();
            Set::iterator lastInsertedLoc = omissing.begin();
            for(int i=0;i<nverts;i++){
                if(it!=end && i==*it){
                    it++;
                    continue;
                }
                //std::cout << i<<" ";
                lastInsertedLoc = omissing.insert(lastInsertedLoc,i);
            }
            useMissingSet = true;
            oobserved = Set();
            //std::cout<<"to miss ";
            omissing.erase(this->idNum);
        }
    }

public:
    DirectedVertex(int netSize){
        nverts = netSize;
        useMissingSet = true;
    }
    virtual ~DirectedVertex(){}

    bool addInedge(int from){
        return iedges.insert(from).second;
    }
    bool addOutedge(int to){
        return oedges.insert(to).second;
    }

    bool hasInedge(int from){
        Set::iterator it = iedges.find(from);
        return it!=iedges.end();
    }

    bool hasOutedge(int to){
        Set::iterator it = oedges.find(to);
        return it!=oedges.end();
    }

    int indegree(){ return iedges.size();}
    int outdegree(){ return oedges.size();}

    bool removeInedge(int from){
        return iedges.erase(from)==1;
    }
    bool removeOutedge(int to){
        return oedges.erase(to)==1;
    }

    const Set& inedges() const{
        return iedges;
    }
    const Set& outedges() const{
        return oedges;
    }

    void clearInedges(){
        iedges.clear();
    }

    void clearOutedges(){
        oedges.clear();
    }

    int networkSize() const{return nverts;}
    void setNetworkSize(int netSize){nverts = netSize;}


    bool isOutmissing(int to){
        if(to==this->idNum)
            return false;
        if(useMissingSet){
            Set::iterator it = omissing.find(to);
            return it!=omissing.end();
        }else{
            Set::iterator it = oobserved.find(to);
            return it==oobserved.end();
        }
    }

    bool setOutmissing(int to,bool miss){
        if(to==this->idNum)
            return false;
        bool ret;
        if(miss){
            if(useMissingSet)
                ret = !omissing.insert(to).second;
            else
                ret = oobserved.erase(to)==0;
        }else{
            if(!useMissingSet)
                ret = oobserved.insert(to).second;
            else
                ret = omissing.erase(to)!=0;
        }
        refreshMissingRepresentation();
        return ret;
    }

    void setAllMissing(){
        useMissingSet=false;
        omissing=Set();
        oobserved=Set();
    }

    void setAllObserved(){
        useMissingSet=true;
        omissing=Set();
        oobserved=Set();
    }

    Set outmissing() const{

        if(useMissingSet){
            Set tmp = omissing;
            tmp.erase(this->idNum);
            return tmp;
        }else{
            Set tmp = Set();
            Set::const_iterator it = oobserved.begin();
            for(int i=0;i<nverts;i++){
                if(it!=oobserved.end() && i==*it){
                    //std::cout<<i<<" ";
                    it++;
                    continue;
                }
                tmp.insert(tmp.end(),i);
                //std::cout<<nverts<<" here";
                //return Set();
            }
            tmp.erase(this->idNum);
            return tmp;
        }
    }

    int nMissing(){
        if(useMissingSet)
            return omissing.size();
        else
            return nverts - 1 - oobserved.size();
    }

    int randomMissingDyad(){
        assert(nMissing()>0);
        double percMissing = nMissing()/(nverts-1.0);

        if(percMissing>0.05){
            for(int i=0;i<15;i++){
                int nbr = floor(Rf_runif(0,nverts-1.0));
                if(nbr>=this->idNum)
                    nbr++;
                if(isOutmissing(nbr))
                    return nbr;
            }
        }
        int index = floor(Rf_runif(0,(double)nMissing()));
        if(useMissingSet){
            Set::iterator it = omissing.begin();
            for(int i=0;i<index;i++)
                it++;
            return *it;
        }else{
            Set::iterator it = oobserved.begin();
            for(;it!=oobserved.end();it++){
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




}
#endif /* DVERTEXH_ */
