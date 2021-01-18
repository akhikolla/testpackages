#ifndef SpeciesTree_h
#define SpeciesTree_h

#include "Tree.h"
#include <sstream>
#include <map>
#include <set>
#include <RcppArmadillo.h>

using namespace Rcpp;

class SpeciesTree : public Tree
{
    private:

        double        speciationRate, extinctionRate;
        unsigned      extantStop;

    public:
                      SpeciesTree(unsigned numTaxa, double curTime, double specRate, double extRate);
                      SpeciesTree(unsigned numTaxa);
                      SpeciesTree(SEXP rtree);
                      SpeciesTree(const SpeciesTree& speciestree, unsigned numTaxa);
        virtual       ~SpeciesTree();

        std::shared_ptr<SpeciesTree>  clone() const { return std::shared_ptr<SpeciesTree>(new SpeciesTree(*this)); }
        void          setSpeciationRate(double sr) {speciationRate = sr; }
        void          setExtinctionRate(double er) {extinctionRate = er; }
        void          setCurrentTime(double et) { currentTime = et; }
        // tree-building functions
        double        getTimeToNextEvent() override;
        double                getTimeToNextEventMoran();
        void          lineageBirthEvent(unsigned indx) override;
        void          lineageDeathEvent(unsigned indx) override;
        void          ermEvent(double curTime) override;
        void          setNewLineageInfo(unsigned indx,
                                        std::shared_ptr<Node> r,
                                        std::shared_ptr<Node> l);

        // set node parameters across tree
        void          setBranchLengths() override;
        void          setPresentTime(double currentT);
        void          setTreeTipNames() override;
        void          recTipNamer(std::shared_ptr<Node> p, unsigned &extinctCount, unsigned &tipCount);

        // simulation functions
        void          setGSATipTreeFlags();
        void          reconstructTreeFromGSASim(std::shared_ptr<Node> oRoot);
        void          setTreeInfo();
        void          popNodes();
        void          recPopNodes(std::shared_ptr<Node> p);
        void          reconstructLineageFromGSASim(std::shared_ptr<Node> currN,
                                                   std::shared_ptr<Node> prevN,
                                                   unsigned &tipCounter,
                                                   unsigned &intNodeCounter);
      //  void          setSampleFromFlags();
        std::map<int,int>           makeIndxMap();
        std::map<int, std::string>  makeTipMap();
        std::map<int,double>        getBirthTimesFromNodes();
        std::map<int,double>        getDeathTimesFromNodes();
        double                      getCurrentTimeFromExtant() {return extantNodes[0]->getDeathTime();}
        double                      getCurrentTime();
        bool                        getIsExtantFromIndx(int indx) { return nodes[indx]->getIsExtant(); }
        bool                        macroEvent(int indx);

        std::pair<int, int>         preorderTraversalStep(int index);
        int                         postOrderTraversalStep(int index);
        int           findLastToGoExtinct(double eventTime);
        int           getNodesIndxFromExtantIndx(int extanIndx) {return extantNodes[extanIndx]->getIndex(); }
        friend class LocusTree;

};


#endif /* SpeciesTree_h */
