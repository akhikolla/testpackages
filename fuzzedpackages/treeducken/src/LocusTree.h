#ifndef LocusTree_h
#define LocusTree_h

#include "SpeciesTree.h"
#include <algorithm>
#include <set>

class LocusTree : public Tree
{
    private:
        double geneBirthRate, geneDeathRate, transferRate;
        double currentTime;
        double stopTime;
        unsigned numTaxa;
        unsigned numTransfers;
        unsigned numDuplications;
        std::vector<std::string> speciesNames;

    public:
        LocusTree(unsigned nt, double stop, double gbr, double gdr, double lgtr);
        LocusTree(const LocusTree& locustree, unsigned numTaxa);
        LocusTree(const SpeciesTree& speciestree, unsigned numTaxa, double gbr, double gdr, double ltr);
        virtual         ~LocusTree();
        double  getTimeToNextEvent() override;
        void    lineageBirthEvent(unsigned indx) override;
        void    lineageDeathEvent(unsigned indx) override;
        void    setNewLineageInfo(int indx, std::shared_ptr<Node> r, std::shared_ptr<Node> s);
        void    lineageTransferEvent(int indx, bool randTrans);
        void    ermEvent(double ct) override;

        int     speciationEvent(int indx, double time, std::pair<int,int> sibs);
        void    extinctionEvent(int indx, double time);
        void    setNewIndices(int indx, std::pair<int,int> sibs, int count);
        void    setSpeciesNames(std::vector<std::string> spNames) { speciesNames = spNames; }
        std::vector<std::string> getSpeciesNames() {return speciesNames;}
        std::string   printNewickTree();
        void    setTreeTipNames() override;
        void    recTipNamer(std::shared_ptr<Node> p, unsigned &copyNumber);
        void    recGetNewickTree(std::shared_ptr<Node> r, std::stringstream &ss);
        void    setBranchLengths() override;
        void    setPresentTime(double currentT);
        void    setStopTime(double st) {stopTime = st; currentTime = 0;}
        double  getCurrentTime() { return currentTime; }
        void    setCurrentTime(double ct) {currentTime = ct; }
        int     getNumberTransfers();
        unsigned int     getNumberDuplications() {return numDuplications;}
        int     chooseRecipientSpeciesID(std::shared_ptr<Node> s);
        std::map<int,double>     getBirthTimesFromNodes();
        std::set<int>            getExtLociIndx();
        std::set<int>            getCoalBounds();
        std::multimap<int,double>     getDeathTimesFromNodes();
        std::multimap<int,double>     getDeathTimesFromExtinctNodes();
        std::map<int,int>             getLocusToSpeciesMap();
        std::vector< std::vector<int> >     getExtantLoci(std::set<double, std::greater<double> > epochSet);
        std::vector< std::string >    printSubTrees();
        int     postOrderTraversalStep(int indx);
        void   setNamesBySpeciesID(std::map<int,std::string> tipMap);
        void   recursiveSetNamesBySpeciesID(std::shared_ptr<Node> n,
                                            int duplicationCount,
                                            std::map<int, std::string> tipMap);
        int    calculatePatristicDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2) override;

        bool   checkLocusTreeParams();
        friend class SpeciesTree;

};
#endif /* LocusTree_h*/
