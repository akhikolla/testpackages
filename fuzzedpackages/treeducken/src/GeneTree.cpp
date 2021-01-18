//
//  GeneTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "GeneTree.h"
#include <iostream>
#include <cmath>

#include "math.h"

#include "math.h"
#include <Rcpp.h>

GeneTree::GeneTree(unsigned nt, unsigned ipp, double ne, double genTime) : Tree(nt){
    numTaxa = nt;
    individualsPerPop = ipp;
    popSize = ne;
    generationTime = genTime;
}

GeneTree::~GeneTree(){

}


//TODO:  go back and speed this up by removing push_back calls
void GeneTree::initializeTree(std::vector< std::vector<int> > extantLociInd, double presentTime){
    int num_loci_in_prsent = 0;
    nodes.clear();
    extantNodes.clear();
    int k = 0;

    if(extantLociInd[0].empty()) {
        while(extantLociInd[k].empty()) {
            k++;            
        }
        num_loci_in_prsent = extantLociInd[k].size();
    }
    else{
        num_loci_in_prsent = extantLociInd[0].size();
    }
    for(int i = 0; i < num_loci_in_prsent; i++){
        for(int j = 0; j < individualsPerPop; j++){
            auto p = std::shared_ptr<Node>(new Node());
            p->setDeathTime(presentTime);
            p->setLindx(extantLociInd[k][i]);
            p->setLdes(NULL);
            p->setRdes(NULL);
            p->setAnc(NULL);
            p->setIsExtant(true);
            p->setIsTip(true);
            p->setIsExtinct(false);
            extantNodes.push_back(p);
            nodes.push_back(p);
            p->setIndx((int) nodes.size());
        }
    }

}

double GeneTree::getCoalTime(int n){
    double ct = NAN;
    double lambda = (double)(n * (n - 1)) / (popSize) ;
    ct = -log(unif_rand()) / (lambda);
    return ct;
}

bool GeneTree::censorCoalescentProcess(double startTime, double stopTime, int contempSpeciesIndx, int ancSpIndx, bool chck){
    int leftInd = 0;
    int rightInd = 0;
    int leftIndExtN = 0;
    int rightIndExtN = 0;
    int extIndx = 0;
    std::shared_ptr<Node> l = nullptr;
    std::shared_ptr<Node> r = nullptr;
    std::shared_ptr<Node> n = nullptr;
    double t = startTime;
    bool all_coalesced = false;
    // search extantNodes for members with Lindx = contempSpecisIndx
    std::vector<int> indInExtNodes;
    for(auto it = extantNodes.begin(); it != extantNodes.end(); ++it){
        if((*it)->getLindx() == contempSpeciesIndx){
            extIndx = std::distance(extantNodes.begin(), it);
            indInExtNodes.push_back(extIndx);
        }
    }
  // the coalescent part
    if(indInExtNodes.size() > 1){
        while(t > stopTime){
            t -= getCoalTime(indInExtNodes.size()); // dra a time
            // is the time older than the end point?
            if(t < stopTime){
                if(chck){
                    t = stopTime;
                    all_coalesced = true;
                    break;
                }
                else{
                    t = stopTime;
                    all_coalesced = false;
                    break;
                }
            }
            // randomly choose two nodes to coalesce in this locus
            rightInd = unif_rand() * (indInExtNodes.size() - 1);
            iter_swap(indInExtNodes.begin() + rightInd, indInExtNodes.begin());
            rightIndExtN = indInExtNodes[0];
            r = extantNodes[rightIndExtN];

            std::reverse(indInExtNodes.begin(), indInExtNodes.end());
            leftInd = unif_rand() * (indInExtNodes.size() - 2);
            iter_swap(indInExtNodes.begin() + leftInd, indInExtNodes.begin());
            leftIndExtN = indInExtNodes[0];
            l = extantNodes[leftIndExtN];
            std::reverse(indInExtNodes.begin(), indInExtNodes.end());
            // do the coalescing
            n = coalescentEvent(t, l, r);
            // book-keeping that almost certainly could be done better:
            // delete two nodes from erase
            if(leftIndExtN > rightIndExtN){
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
            }
            else{
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
            }
            // clear the indices in extant nodes vector
            indInExtNodes.clear();
            // put in the new node
            extantNodes.insert(extantNodes.begin(), n);

            // populate the indices in extant nodes vector again, this time to see if it is empty
            for(auto it = extantNodes.begin(); it != extantNodes.end(); ++it){
                if((*it)->getLindx() == contempSpeciesIndx){
                    extIndx = std::distance(extantNodes.begin(), it);
                    indInExtNodes.push_back(extIndx);
                }
            }
            // if only one is left get out of the loop
            if(indInExtNodes.size() == 1){
                all_coalesced = true;
                break;
            }
        }
    }
    // only one member so nothing to do besides progress time
    else if (indInExtNodes.size() == 1){
        t = stopTime;
        all_coalesced = true;
        extantNodes[indInExtNodes[0]]->setLindx(ancSpIndx);
    }
    else{
        // this is 0 to catch any stragglers and in Simulator::simulateCoalescentProcess those will be deleted from the contempSpecies listing
        all_coalesced = true;
    }
    // if everything coalesced loop through and change the locus indices in geneTree.nodes
    // TODO: refactor this to make clear when species indices are being used (they aren't) and locus indices are (they are)
    if(all_coalesced == true){
        for(int i = 0; i < indInExtNodes.size(); ++i){
            extantNodes[indInExtNodes[i]]->setLindx(ancSpIndx);
        }
    }
    // clear this again
    indInExtNodes.clear();
    return all_coalesced;
}

std::shared_ptr<Node> GeneTree::coalescentEvent(double t, 
                                std::shared_ptr<Node> p,
                                std::shared_ptr<Node> q){
    auto n = std::shared_ptr<Node>(new Node());
    n->setDeathTime(t);
    n->setLdes(p);
    n->setRdes(q);
    n->setIsExtant(false);
    n->setIsTip(false);
    n->setIsExtinct(false);
    n->setLindx(p->getLindx());
    nodes.push_back(n);
    n->setIndx((int) nodes.size());
    p->setBirthTime(t);
    p->setAnc(n);
   // p->setSib(q);

    q->setBirthTime(t);
    q->setAnc(n);
    // q->setSib(p);


    return n;
}

std::multimap<int, double> GeneTree::rescaleTimes(std::multimap<int, double> timeMap){
    std::multimap<int, double> rescaledTimeMap;
    std::pair<int, double> p;
    for(std::multimap<int, double>::iterator it = timeMap.begin(); it != timeMap.end(); ++it){
        p.first = (*it).first;
        p.second = ((*it).second);
        rescaledTimeMap.insert(p);
    }

    return rescaledTimeMap;

}


void GeneTree::rootCoalescentProcess(double startTime){
    double t = startTime;
    for(auto en : extantNodes){
        en->setLindx(0);
    }
    while(extantNodes.size() > 1){
        t -= getCoalTime(extantNodes.size());

        int rightInd = unif_rand() * (extantNodes.size() - 1);
        auto r = std::move(extantNodes[rightInd]);
        extantNodes.erase(extantNodes.begin() + rightInd);

        int leftInd = unif_rand() * (extantNodes.size() - 1);
        auto l = std::move(extantNodes[leftInd]);
        extantNodes.erase(extantNodes.begin() + leftInd);

        auto n = coalescentEvent(t, l, r);
        extantNodes.push_back(n);
    }
    extantNodes[0]->setAsRoot(true);
    extantNodes[0]->setBirthTime(t);
    setRoot(extantNodes[0]);
}

void GeneTree::recursiveRescaleTimes(std::shared_ptr<Node> r, double add){
    if(r != NULL){
        if( r->getRdes() == NULL){
            r->setBirthTime(r->getBirthTime() + add);
            r->setDeathTime(r->getDeathTime() + add);
        }
        else{

            r->getLdes()->setBirthTime(r->getLdes()->getBirthTime() + add);
            r->getLdes()->setDeathTime(r->getLdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getLdes(), add);

            r->getRdes()->setBirthTime(r->getRdes()->getBirthTime() + add);
            r->getRdes()->setDeathTime(r->getRdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getRdes(), add);

        }
    }
}

void GeneTree::setBranchLengths(){
    double brlen = NAN;
    numExtant = 0;
    numExtinct = 0;
    for(auto node : nodes){
        brlen = node->getDeathTime() - node->getBirthTime();
        node->setBranchLength(brlen);

        if(node->getIsTip()){
          if(node->getIsExtant())
            numExtant++;
          else
            numExtinct++;
          branchLengths.push_back(brlen);
        }
        else if(node->getIsRoot()){
          branchLengths.emplace(branchLengths.begin(), brlen);
        }
        else{
          branchLengths.push_back(brlen);
        }
    }
    this->setTreeTipNames();
}

void GeneTree::addExtinctSpecies(double bt, int indx){
    std::shared_ptr<Node> p = nullptr;
    for(int i = 0; i < individualsPerPop; i++){
        std::shared_ptr<Node> p = std::shared_ptr<Node>(new Node());
        p->setDeathTime(bt);
        p->setLindx(indx);
        p->setLdes(NULL);
        p->setRdes(NULL);
        p->setAnc(NULL);
        p->setIsExtant(false);
        p->setIsTip(true);
        p->setIsExtinct(true);
        extantNodes.push_back(p);
        nodes.push_back(p);
        p->setIndx((int) nodes.size() + 1);

    }
    //delete p;
}


void GeneTree::setIndicesBySpecies(std::map<int, int> spToLocusMap){
    numExtant = 0;
    numExtinct = 0;
    for(auto node : nodes){
        if(node->getIsTip()){
            // indx = (*it)->getLindx();
            // spIndx = spToLocusMap.find(indx)->second;
            //(*it)->setIndx((*it)->getLindx());

            if(node->getIsExtant())
              numExtant++;
            else
              numExtinct++;
        }
    }
    for(auto node : nodes){
      unsigned int notTipCount = numExtant + numExtinct + 1;
      if(!(node->getIsTip())){
        node->setIndx(notTipCount);
        notTipCount++;
      }
    }
}



void GeneTree::setTreeTipNames(){
    int indNumber = 0;
    std::stringstream tn;
    std::string name;
    int locusIndxCounter = 0;

    for(auto node : nodes){
        if(node->getIsTip()){
            tn << locusIndxCounter + 1;
            name = tn.str();
            tn.clear();
            tn.str(std::string());
            indNumber++;
            tn << indNumber;
            name += "_" + tn.str();
            tn.clear();
            tn.str(std::string());
            node->setName(name);
            if(indNumber == individualsPerPop){
              indNumber = 0;
              locusIndxCounter++;
            }
        }

    }

}

void GeneTree::reindexForR(){
  unsigned int intNodeCount = numExtant + numExtinct + 1;
  int tipCount = 1;
  for(int i = nodes.size() - 1; i > -1; i--){
    if(nodes[i]->getIsTip()){
      nodes[i]->setIndx(tipCount);
      tipCount++;
    }
    else{
      nodes[i]->setIndx(intNodeCount);
      intNodeCount++;
    }
  }
}


NumericMatrix GeneTree::getGeneEdges(){
  this->GeneTree::reindexForR();
  int numRows = (int) nodes.size() - 1;
  NumericMatrix edgeMat(numRows, 2);
  for(int i=0; i < nodes.size()-1; i++){
    if(!(nodes[i]->getIsRoot())){

      NumericMatrix::Row row = edgeMat(i, _);

      row[0] = nodes[i]->getAnc()->getIndex();
      row[1] = nodes[i]->getIndex();
    }
  }
  return edgeMat;
}

