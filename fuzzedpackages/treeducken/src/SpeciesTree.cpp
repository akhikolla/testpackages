#include "SpeciesTree.h"
#include <iostream>

#include "math.h"

using namespace Rcpp;


SpeciesTree::SpeciesTree(unsigned numTaxa, double ct, double br, double dr) : Tree(numTaxa, 0.0){
    extantStop = numTaxa;
    speciationRate = br;
    extinctionRate = dr;

}

SpeciesTree::SpeciesTree(unsigned numTaxa) : Tree(numTaxa){
    extantStop = numTaxa;
}

SpeciesTree::SpeciesTree(SEXP rtree) : Tree(rtree){
  speciationRate = 0.0;
  extinctionRate = 0.0;
}

SpeciesTree::SpeciesTree(const SpeciesTree& speciestree, unsigned numTaxa) : Tree(numTaxa) {
  extantStop = numTaxa;
  nodes = speciestree.nodes;
  extantNodes = speciestree.extantNodes;
  root = speciestree.root;
  speciationRate = speciestree.speciationRate;
  extinctionRate = speciestree.extinctionRate;
  extantStop = speciestree.extantStop;
  extantRoot = speciestree.extantRoot;
  currentTime = speciestree.currentTime;
  numNodes = speciestree.numNodes;
  numTotalTips = speciestree.numTotalTips;
  numExtant = speciestree.numExtant;
  numExtinct = speciestree.numExtinct;
}


SpeciesTree::~SpeciesTree(){

}

double SpeciesTree::getTimeToNextEvent(){
    double sumrt = speciationRate + extinctionRate;
    double returnTime = 0.0;
    NumericVector randNum = Rcpp::runif(1);

    returnTime = -log(randNum[0]) / (double(numExtant) * sumrt);
    return returnTime;
}

void SpeciesTree::lineageBirthEvent(unsigned indx){
    std::shared_ptr<Node> sis;
    std::shared_ptr<Node> right;
    right = std::shared_ptr<Node>(new Node());
    sis = std::shared_ptr<Node>(new Node());
    setNewLineageInfo(indx, right, sis);
}

void SpeciesTree::lineageDeathEvent(unsigned int indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void SpeciesTree::ermEvent(double cTime){
    currentTime = cTime;
    RNGScope scope;
    NumericVector randNum = Rcpp::runif(2);
    int nodeInd = randNum[0]*(numExtant - 1);
    double relBr = speciationRate / (speciationRate + extinctionRate);
    bool isBirth = (randNum[1] < relBr ? true : false);
    if(isBirth)
        lineageBirthEvent(nodeInd);
    else
        lineageDeathEvent(nodeInd);
}

void SpeciesTree::setNewLineageInfo(unsigned int indx,
                                    std::shared_ptr<Node> r,
                                    std::shared_ptr<Node> l){
    extantNodes[indx]->setLdes(l);
    extantNodes[indx]->setRdes(r);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsExtant(false);

    r->setLdes(NULL);
    r->setRdes(NULL);
    r->setSib(l);
    r->setAnc(extantNodes[indx]);
    r->setBirthTime(currentTime);
    r->setIsTip(true);
    r->setIsExtant(true);
    r->setIsExtinct(false);

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);

    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(r);
    extantNodes.push_back(l);
    nodes.push_back(r);
    nodes.push_back(l);
    numNodes = (int) nodes.size();
    numExtant = (int) extantNodes.size();
    r->setIndx(numNodes - 2);
    l->setIndx(numNodes - 1);

}

void SpeciesTree::setBranchLengths(){
    double bl = NAN;
    for(auto node : nodes){
      bl = node->getDeathTime() - node->getBirthTime();
      branchLengths.push_back(bl);
      node->setBranchLength(bl);
    }
}

void SpeciesTree::setPresentTime(double currentT){
    for(auto extantNode : extantNodes){
        extantNode->setDeathTime(currentT);
        extantNode->setIsExtant(true);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}

void SpeciesTree::setTreeInfo(){
  //  double trDepth = this->getTreeDepth();
    std::set<double> deathTimes;
    auto it = nodes.begin();
    (*it)->setBirthTime(0.0);
    (*it)->setDeathTime((*it)->getBranchLength() + (*it)->getBirthTime());
    (*it)->setIndx(0);
    ++it;
    for(; it != nodes.end(); ++it){
        (*it)->setBirthTime((*it)->getAnc()->getDeathTime());
        (*it)->setDeathTime((*it)->getBranchLength() + (*it)->getBirthTime());
        deathTimes.insert(deathTimes.begin(),(*it)->getBranchLength() + (*it)->getBirthTime());
        (*it)->setIndx((int)std::distance(nodes.begin(), it));
    }
    it = nodes.begin();
    std::set<double>::iterator set_iter = deathTimes.end();
    --set_iter;
    double currentTime = *(set_iter);
    for(; it != nodes.end(); ++it){
        if((*it)->getIsTip()){
          auto placeholder = 0.1;
          if (std::abs((*it)->getDeathTime() - currentTime) < placeholder) {
            (*it)->setIsExtant(true);
            (*it)->setIsExtinct(false);
            (*it)->setDeathTime(currentTime);
            numTaxa++;
            extantNodes.push_back(*it);
          } else {
            (*it)->setIsExtant(false);
            (*it)->setIsExtinct(true);
          }
        }
    }
    return;
}

void SpeciesTree::setTreeTipNames(){
  unsigned nodeIndx = numExtant + numExtinct;
  unsigned tipIt = 0;
  std::stringstream tn;

  for(unsigned int i=0; i < nodes.size(); i++){
    if(nodes[i]->getIsTip()){
      tipIt++;
      nodes[i]->setIndx(tipIt);
      if(nodes[i]->getIsExtant()){
        tn << nodes[i]->getIndex();
        std::string name = "H" + tn.str();
        nodes[i]->setName(name);

      }
      else{
        tn << nodes[i]->getIndex();
        std::string name = "X" + tn.str();
        nodes[i]->setName(name);
      }
    }
    else{
      nodeIndx++;
      nodes[i]->setIndx(nodeIndx);
    }
    tn.clear();
    tn.str(std::string());
  }
}


void SpeciesTree::recTipNamer(std::shared_ptr<Node> p,
                              unsigned &nodeIndx, 
                              unsigned &tipIndx){
  if(p != NULL){
    std::stringstream tn;
    if(p->getIsTip()){
      tipIndx++;
      p->setIndx(tipIndx);
      if(p->getIsExtinct()){
        tn << p->getIndex();
        std::string name = "X" + tn.str();
        p->setName(name);

      }
      else{
        tn << p->getIndex();
        std::string name = "H" + tn.str();
        p->setName(name);
      }
    }
    else{
      nodeIndx++;
      p->setIndx(nodeIndx);
      recTipNamer(p->getLdes(), nodeIndx, tipIndx);
      recTipNamer(p->getRdes(), nodeIndx, tipIndx);

    }
  }
}


void SpeciesTree::setGSATipTreeFlags(){
    zeroAllFlags();
    numTotalTips = 0;
    for(auto node : nodes){
        if(node->getIsTip()){
            numTotalTips++;
            node->setFlag(1);

        }
        else{
            node->setFlag(2);
        }
    }
    setSampleFromFlags();
}


void SpeciesTree::popNodes(){
    nodes.clear();
    extantNodes.clear();

    recPopNodes(this->getRoot());

}

void SpeciesTree::recPopNodes(std::shared_ptr<Node> p){
    if(p != nullptr){
        if(p->getIsTip()){
            if(p->getIsExtant()){
                extantNodes.push_back(p);
                nodes.push_back(p);
            }
            else{
                nodes.push_back(p);
            }
        }
        else{
            nodes.push_back(p);
            recPopNodes(p->getLdes());
            recPopNodes(p->getRdes());
        }
    }
}

void SpeciesTree::reconstructTreeFromGSASim(std::shared_ptr<Node> oRoot){
    std::shared_ptr<Node> n = std::shared_ptr<Node>(new Node());
    unsigned tipCounter = 0;
    unsigned intNodeCounter = extantStop;
    reconstructLineageFromGSASim(n, oRoot, tipCounter, intNodeCounter);
}

void SpeciesTree::reconstructLineageFromGSASim(std::shared_ptr<Node> currN, 
                                               std::shared_ptr<Node> prevN, 
                                               unsigned &tipCounter, 
                                               unsigned &intNodeCounter){
    std::shared_ptr<Node> p = nullptr;
    bool rootN = prevN->getIsRoot();
    double brlen = prevN->getBranchLength();
    int oFlag = prevN->getFlag();
    if(prevN->getIsTip() && oFlag == 1){
        // need to recalculate branchlength
        std::shared_ptr<Node> prevAnc = prevN->getAnc();
        int ancFlag = prevAnc->getFlag();
        if(ancFlag == 1){
            brlen += prevAnc->getBranchLength();
            while(!prevAnc->getIsRoot() && ancFlag < 2){
                prevAnc = prevAnc->getAnc();
                ancFlag = prevAnc->getFlag();
                if(ancFlag == 1)
                    brlen += prevAnc->getBranchLength();
            }
        }

        std::shared_ptr<Node> p = std::shared_ptr<Node>(new Node());
        tipCounter++;
        p->setIndx(tipCounter);
        p->setBranchLength(brlen);
        p->setIsTip(true);
        p->setBirthTime(prevN->getBirthTime());
        p->setDeathTime(prevN->getDeathTime());
        p->setIsExtant(prevN->getIsExtant());
        p->setIsExtinct(prevN->getIsExtinct());
        p->setAnc(currN);
        if(currN->getLdes() == NULL)
            currN->setLdes(p);
        else if(currN->getRdes() == NULL)
            currN->setRdes(p);
        else{
            stop("ERROR: Problem adding a tip to the tree!");
        }

    }
    else{
        if(oFlag > 1){
            std::shared_ptr<Node> s1 = std::shared_ptr<Node>(new Node());
            intNodeCounter++;
            s1->setIndx(intNodeCounter);
            if(prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromGSASim(s1, prevN->getLdes(), tipCounter, intNodeCounter);
            if(prevN->getRdes()->getFlag() > 0)
                reconstructLineageFromGSASim(s1, prevN->getRdes(), tipCounter, intNodeCounter);


            if(rootN == false){
                std::shared_ptr<Node> prevAnc = prevN->getAnc();
                int ancFlag = prevAnc->getFlag();
                if(ancFlag == 1){
                    brlen += prevAnc->getBranchLength();
                    while(!(prevAnc)->getIsRoot() && ancFlag < 2){
                        prevAnc = prevAnc->getAnc();
                        ancFlag = prevAnc->getFlag();
                        if(ancFlag == 1)
                            brlen += prevAnc->getBranchLength();
                    }
                }

                if(currN != NULL){
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                    s1->setAnc(currN);
                    if(currN->getLdes() == NULL)
                        currN->setLdes(s1);
                    else if(currN->getRdes() == NULL)
                        currN->setRdes(s1);
                    else{
                        stop("ERROR: Probem adding an internal node to the tree");
                    }
                }
                else{
                    s1->setAsRoot(true);
                    setRoot(s1);
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                }

            }
            else{
                s1->setAsRoot(true);
                setRoot(s1);
                s1->setBranchLength(0.0);
                s1->setBirthTime(prevN->getBirthTime());
                s1->setDeathTime(prevN->getDeathTime());
            }

        }
        else if(oFlag == 1){
            if(prevN->getRdes()->getFlag() == 0 && prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromGSASim(currN, prevN->getLdes(), tipCounter, intNodeCounter);
            else
                reconstructLineageFromGSASim(currN, prevN->getRdes(), tipCounter, intNodeCounter);
        }
    }
}


std::map<int,int> SpeciesTree::makeIndxMap(){
  std::map<int,int> indxMap;
  for(unsigned int i=0; i < nodes.size(); i++){
    int rIndx = nodes[i]->getIndex();
    int tdckenIndx = i;
    indxMap.insert(std::pair<int,int>(rIndx, tdckenIndx));
  }
  return indxMap;
}


std::map<int, std::string> SpeciesTree::makeTipMap(){
  std::map<int, std::string> tipMap;
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    if(nodes[i]->getIsTip())
    {
      int j = nodes[i]->getIndex();
      std::string tipName = nodes[i]->getName();
      tipMap.insert(std::pair<int, std::string>(j, tipName));
    }
  }
  return tipMap;
}

std::map<int,double> SpeciesTree::getBirthTimesFromNodes(){
    int indx = -1;
    double birthTime = NAN;
    std::map<int,double> birthTimeMap;
    for(auto node : nodes){
        indx = node->getIndex();
        birthTime = node->getBirthTime();
        birthTimeMap.insert(std::pair<int,double>(indx, birthTime));
    }
    return birthTimeMap;
}

std::map<int,double> SpeciesTree::getDeathTimesFromNodes(){
    int indx = -1;
    double deathTime = NAN;
    std::map<int,double> deathTimeMap;
    for(auto node : nodes){
        if(!(node->getIsExtant())){
            indx = node->getIndex();
            deathTime = node->getDeathTime();

            deathTimeMap.insert(std::pair<int,double>(indx, deathTime));
        }
    }
    return deathTimeMap;
}

std::pair<int,int> SpeciesTree::preorderTraversalStep(int indx){
    std::pair<int,int> sibs;
    sibs.first = nodes[indx]->getLdes()->getIndex();
    sibs.second = nodes[indx]->getRdes()->getIndex();
    return sibs;
}

int SpeciesTree::postOrderTraversalStep(int index){
    int d = -1;
    d = nodes[index]->getAnc()->getIndex();
    return d;
}

bool SpeciesTree::macroEvent(int indx){
    bool isSpec = 0;
    std::shared_ptr<Node> n = nodes[indx];

    if(n->getIsTip())
        isSpec = false;
    else
        isSpec = true;
    return isSpec;
}


int SpeciesTree::findLastToGoExtinct(double EventTime){
  int indxExtinct = -1;
  double epsi = std::numeric_limits<double>::epsilon();
  bool is_near = false;
  for(unsigned int i=0; i < nodes.size(); i++){
    if(nodes[i]->getIsTip() && nodes[i]->getIsExtinct()){
      double scale = std::max(abs(EventTime), abs(nodes[i]->getDeathTime()));
      is_near = abs(nodes[i]->getDeathTime() - EventTime) <= scale *(2*epsi);
      if(is_near){
        indxExtinct = i;
        break;
      }
    }
  }

  return indxExtinct;
}

double SpeciesTree::getCurrentTime() {
    std::vector<double> tempDeathTimes;
    tempDeathTimes.resize(nodes.size());
    for(unsigned int i = 0; i < nodes.size(); i++) {
        tempDeathTimes[i] = nodes[i]->getDeathTime();
    }
    sort(tempDeathTimes.begin(), tempDeathTimes.end(), std::greater<double>());
    return(tempDeathTimes[0]);
}
