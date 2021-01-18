//
//  Tree.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/7/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Tree_h
#define Tree_h

#include <string>
#include <vector>
#include <iostream>
#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

class Node
{
    private:
        std::shared_ptr<Node>    ldes;
        std::shared_ptr<Node>    rdes;
        std::shared_ptr<Node>    anc;
        std::shared_ptr<Node>    sib;
        int     indx, Lindx;
        std::vector<unsigned int> hosts;
        int     flag;
        std::string name;
        bool    isRoot;
        bool    isTip;
        bool    isExtant, isExtinct;
        bool    isDuplication;
        double  birthTime, deathTime;
        double  branchLength;
        int     locusID;

    public:
                Node();
                ~Node();
        void    setAsRoot(bool t) {isRoot = t; }
        void    setBirthTime(double bt) {birthTime = bt; }
        void    setIsTip(bool t) {isTip = t; }
        void    setDeathTime(double dt) {deathTime = dt; }
        void    setIsExtant(bool t) {isExtant = t; }
        void    setIsExtinct(bool t) {isExtinct = t; }
        void    setLdes(std::shared_ptr<Node> l) {ldes = l; }
        void    setRdes(std::shared_ptr<Node> r) {rdes = r; }
        void    setAnc(std::shared_ptr<Node> a) {anc = a; }
        void    setSib(std::shared_ptr<Node> s) {sib = s; }
        void    setName(std::string f) { name = f; }
        void    setBranchLength(double bl) {branchLength = bl; }
        void    setFlag(int d) { flag = d; }
        void    setIndx(int i) {indx = i; }
        void    setLindx(int li ) {Lindx = li; }
        void    addHost(int hostIndx) { hosts.push_back(hostIndx); }
        void    setIsDuplication(bool t) { isDuplication = t; }
        void    setLocusID(int a) { locusID = a; }
        int     getFlag() {return flag; }
        std::shared_ptr<Node>   getLdes() {return ldes; }
        std::shared_ptr<Node>   getRdes() {return rdes; }
        std::shared_ptr<Node>   getAnc() {return anc; }
        std::shared_ptr<Node>   getSib() {return sib; }
        bool    getIsRoot() {return isRoot; }
        bool    getIsTip() {return isTip; }
        bool    getIsExtinct() {return isExtinct; }
        bool    getIsExtant() { return isExtant; }
        std::string getName() { return name; }
        double  getBranchLength() { return branchLength; }
        double  getDeathTime() {return deathTime; }
        double  getBirthTime() { return birthTime; }
        int     getIndex() {return indx; }
        int     getLindx() { return Lindx; }
        std::vector<unsigned int> getHosts() { return hosts; }
        void    setHosts(std::vector<unsigned int> hs) { hosts = hs; }
        bool    getIsDuplication() { return isDuplication; }
        int     getLocusID() { return locusID; }
};



class Tree
{
    protected:
        std::shared_ptr<Node> root;
        std::shared_ptr<Node> extantRoot;
        std::vector<std::shared_ptr<Node>> nodes;
        std::vector<std::shared_ptr<Node>> extantNodes;
        int numTaxa;
        int numNodes;
        int numTotalTips;
        int numExtant, numExtinct;
        double  currentTime;
        std::vector<double> branchLengths;

    public:
                    Tree(unsigned numExtant, double cTime);
                    Tree(unsigned numTaxa);
                    Tree(SEXP rtree);
        virtual      ~Tree();
        std::shared_ptr<Node>    getRoot() {return root; }
        std::shared_ptr<Node>    getExtantRoot() { return extantRoot; }
        void        setExtantRoot(std::shared_ptr<Node> r) { extantRoot = r; }
        void        setRoot(std::shared_ptr<Node> r) { root = r; }
        unsigned int         getNumExtant() {return numExtant; }
        int         getNumTips() { return extantNodes.size(); }
        int         getNumExtinct() {return numExtinct; }
        int         getNodesSize() { return (int) nodes.size(); }
        double      getTotalTreeLength();
        double      getTreeDepth();
        double      getCurrentTime() {return currentTime; }
        double      getEndTime();
        void        setNumExtant();
        void        setNumExtinct();
        void        clearNodes(std::shared_ptr<Node> r);
        void        zeroAllFlags();
        void        setWholeTreeFlags();
        void        setExtantTreeFlags();
        void        setSampleFromFlags();
        void        getRootFromFlags(bool isGeneTree = false);
        void        getExtantTree();

        std::vector<std::shared_ptr<Node>> getNodes() { return nodes; }
        std::vector<std::shared_ptr<Node>> getExtantNodes() { return extantNodes; }
        void        scaleTree(double scVal);
        void        scaleTreeDepthToValue(double scVal);

        void        reconstructTreeFromSim(std::shared_ptr<Node> oRoot);
        void        reconstructLineageFromSim(std::shared_ptr<Node> currN,
                                              std::shared_ptr<Node> prevN,
                                              unsigned &tipCounter,
                                              unsigned &intNodeCounter);


        void        reindexForR();
        std::vector<std::string>    getTipNames();
        std::vector<std::string>    getNodeLabels();
        NumericMatrix getEdges();
        std::vector<double> getEdgeLengths();
        int         getNnodes() { return nodes.size() - (numExtant + numExtinct);}
        void        setTipsFromRtree();
        double      findMaxNodeHeight();
        int         getIndexFromNodes(int indx) {return nodes[indx]->getIndex(); }
        void        switchIndicesFirstToSecond(std::map<int,int> mappy);
        virtual double  getTimeToNextEvent() { return 0.0; }
        virtual void    lineageBirthEvent(unsigned int indx) { return; }
        virtual void    lineageDeathEvent(unsigned int indx) { return; }
        virtual void    setTreeTipNames()  { return; }
        virtual void    ermEvent(double ct) { return; }
        virtual void    setBranchLengths() { return; }
        virtual int     calculatePatristicDistance(std::shared_ptr<Node>, std::shared_ptr<Node>);
        friend class Node;

};
#endif /* Tree_hpp */
