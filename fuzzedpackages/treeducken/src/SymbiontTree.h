//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#ifndef SRC_SYMBIONTTREE_H
#define SRC_SYMBIONTTREE_H

#include "Tree.h"
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include "SpeciesTree.h"

class SymbiontTree : public Tree {

    private:
      double symbSpecRate, symbExtRate, hostExpanRate;
      double currentTime;
      double stopTime;
      unsigned numTaxa;
      unsigned numExpansions;
      unsigned hostLimit;
      std::map<unsigned int,std::vector<unsigned int>> symbHostMap; // keys are symb indices

    public:
      SymbiontTree(int nt,
                   double currSimTime,
                   double symbsr,
                   double symber,
                   double hostExpanRate,
                   int hostLimit);
      SymbiontTree(const SymbiontTree& symbionttree,
                   unsigned numTaxa);

      virtual         ~SymbiontTree();
      double  getTimeToNextJointEvent(double hostSpecRate,
                                         double hostExtRate,
                                         double cospeciaRate,
                                         arma::umat assocMat);
      void    lineageBirthEvent(unsigned indx) override;
      void    lineageDeathEvent(unsigned indx) override;
      virtual void    setNewLineageInfo(unsigned int indx, std::shared_ptr<Node> r, std::shared_ptr<Node> s);
      void            setNewLineageInfoExpan(unsigned int indx,
                                             std::shared_ptr<Node> r,
                                             std::shared_ptr<Node> s,
                                             unsigned int hostIndx);
      void            hostExpansionEvent(unsigned int indx, unsigned int hostIndx);
      arma::umat       ermJointEvent(double ct, arma::umat assocMat);

      void            setSymbTreeInfoSpeciation(unsigned int ancIndx, unsigned int desIndx);
      void            setSymbTreeInfoExtinction(unsigned int deadIndx);

      //std::string     printNewickTree();
      void            setTreeTipNames() override;
      void            recTipNamer(std::shared_ptr<Node> p, unsigned &extinctCount, unsigned &tipCount);

//      void            recGetNewickTree(Node *r, std::stringstream &ss);
      void            setBranchLengths() override;
      void            setPresentTime(double currentT);
      void            setStopTime(double st) { stopTime = st; currentTime = 0;}
      double          getCurrentTime() { return currentTime; }
      void            setCurrentTime(double ct) { currentTime = ct; }
      int             getNumberExpansions() {return numExpansions; }
      int             getNumHostSymbPairs() { return symbHostMap.size(); }

      std::vector<unsigned int>  getSymbsOnHost(unsigned int hostIndx);
      void            updateCurrentMap(unsigned int oldHostIndx, unsigned int newHostIndx);
      int             getExtantIndxFromNodes(unsigned int extantNodesIndx);
      void            cospeciationMapUpdate(unsigned int oldHostIndx,
                                            unsigned int numHosts,
                                            unsigned int symbIndx);
      void            updateHostsInNodes();
      unsigned int             getNodesIndxFromExtantIndx(int i) { return extantNodes[i]->getIndex(); }

};


#endif //SRC_SYMBIONTTREE_H
