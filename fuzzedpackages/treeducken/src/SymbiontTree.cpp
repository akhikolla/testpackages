//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#include "SymbiontTree.h"

#include "math.h"
SymbiontTree::SymbiontTree(int nt,
                           double ct,
                           double br,
                           double dr,
                           double her,
                           int K) : Tree(nt, 0.0){
    numTaxa = nt;
    symbSpecRate = br;
    symbExtRate = dr;
    hostExpanRate = her;
    numExpansions = 0;
    hostLimit = K;
    root->addHost(0);
    std::vector<unsigned> initialHosts;
    initialHosts.resize(1);
    initialHosts[0] = 0;
    symbHostMap.insert(std::pair<unsigned,std::vector<unsigned>> (0,initialHosts));
}

SymbiontTree::SymbiontTree(const SymbiontTree& symbionttree, unsigned numTaxa) : Tree(numTaxa) {
    nodes = symbionttree.nodes;
    extantNodes = symbionttree.extantNodes;
    root = symbionttree.root;
    symbSpecRate = symbionttree.symbSpecRate;
    symbExtRate = symbionttree.symbExtRate;
    extantRoot = symbionttree.extantRoot;
    currentTime = symbionttree.currentTime;
    numNodes = symbionttree.numNodes;
    numTotalTips = symbionttree.numTotalTips;
    numExtant = symbionttree.numExtant;
    numExtinct = symbionttree.numExtinct;
}

SymbiontTree::~SymbiontTree(){

}


double SymbiontTree::getTimeToNextJointEvent(double hostSpecRate,
                                        double hostExtRate,
                                        double cospeciaRate,
                                        arma::umat assocMat){
    int numHosts = assocMat.n_cols;
// arma::uvec numPairs = find(assocMat);
    double sumrt_host =  (hostSpecRate + hostExtRate) * numHosts;
    double sumrt_symb = (symbSpecRate + symbExtRate + hostExpanRate) * numExtant;
    double sumrt_both = cospeciaRate * numHosts;
    double returnTime = -log(unif_rand()) / (sumrt_host + sumrt_symb + sumrt_both);
    return returnTime;
}

void SymbiontTree::setSymbTreeInfoSpeciation(unsigned int indxToFind, unsigned int indxToReplace){
    for(auto s = extantNodes.begin(); s != extantNodes.end(); ++s){
        std::vector<unsigned int> hostsOfS = (*s)->getHosts();
        for(auto hosts : hostsOfS){
            if(hosts == indxToFind)
                hosts = indxToReplace;
        }
    }
}

void SymbiontTree::setSymbTreeInfoExtinction(unsigned int deadIndx){
    std::vector<unsigned int> toBeExtincted;
    for(unsigned int i = 0; i < extantNodes.size(); i++){
        std::vector<unsigned int> hostsOf = extantNodes[i]->getHosts();
        for(unsigned int j = 0; j < hostsOf.size(); j++){
            if(hostsOf[j] == deadIndx){
                std::swap(hostsOf.back(), j);
                hostsOf.pop_back();
            }
             
        }
        if(hostsOf.empty())
            toBeExtincted.push_back(i);
    }
    if(!(toBeExtincted.empty())){
        for(unsigned int i = 0; i < toBeExtincted.size(); i++){
            this->lineageDeathEvent(toBeExtincted[i]);
        }
    }
}

std::vector<unsigned int> SymbiontTree::getSymbsOnHost(unsigned int hostIndx){
    auto symbs = symbHostMap[hostIndx];
    return symbs;
}

void SymbiontTree::lineageBirthEvent(unsigned indx){
    std::shared_ptr<Node> right = std::shared_ptr<Node>(new Node());
    std::shared_ptr<Node> sis = std::shared_ptr<Node>(new Node());
    setNewLineageInfo(indx, right, sis);
}

void SymbiontTree::lineageDeathEvent(unsigned indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void SymbiontTree::setNewLineageInfo(unsigned int indx, std::shared_ptr<Node> r, std::shared_ptr<Node> l){
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
    //r->setHosts(extantNodes[indx]->getHosts());
    //std::vector<int> hostCheckVec = r->getHosts();

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);
    //l->setHosts(extantNodes[indx]->getHosts());

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

arma::umat SymbiontTree::ermJointEvent(double ct, arma::umat assocMat){
    currentTime = ct;
    this->setCurrentTime(ct);

    // pick a row at random
    int nodeInd = unif_rand()*(numExtant);

    arma::urowvec rvec = assocMat.row(nodeInd);
    assocMat.shed_row(nodeInd);

    // which event
    double relBr = symbSpecRate / (symbExtRate + symbSpecRate + hostExpanRate);
    double relDr = relBr + (symbExtRate / (symbExtRate + symbSpecRate + hostExpanRate));
    double dec = unif_rand();
    if(dec < relBr){
        // its a birth
        this->lineageBirthEvent(nodeInd);
        assocMat.resize(numExtant, assocMat.n_cols);
        assocMat(numExtant - 2, arma::span::all) = rvec;
        assocMat(numExtant - 1, arma::span::all) = rvec;
    }
    else if(dec < relDr)
        this->lineageDeathEvent(nodeInd);

    else{
        int hostInd = unif_rand() * assocMat.n_cols;
        this->hostExpansionEvent(nodeInd, hostInd);
        assocMat.resize(numExtant, assocMat.n_cols);
        assocMat(numExtant - 2, arma::span::all) = rvec;
        rvec(hostInd) = 1;
        assocMat(numExtant - 1, arma::span::all) = rvec;
    }
    return assocMat;
}

void SymbiontTree::hostExpansionEvent(unsigned int indx, unsigned int hostIndx){
    std::shared_ptr<Node> right = std::shared_ptr<Node>(new Node());
    std::shared_ptr<Node> sis = std::shared_ptr<Node>(new Node());
    this->setNewLineageInfoExpan(indx, right, sis, hostIndx);
}

void SymbiontTree::setNewLineageInfoExpan(unsigned int indx,
                                         std::shared_ptr<Node> r, 
                                         std::shared_ptr<Node> l, 
                                         unsigned int hostIndx){
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
    numExtant = (int) extantNodes.size();
    r->setIndx(numExtant - 2);
    l->setIndx(numExtant - 1);
}

void SymbiontTree::setTreeTipNames(){
    unsigned nodeIndx = numExtant + numExtinct;
    unsigned tipIt = 0;
    std::stringstream tn;

    for(unsigned int i=0; i < nodes.size(); i++){
        if(nodes[i]->getIsTip()){
            tipIt++;
            nodes[i]->setIndx(tipIt);
            if(nodes[i]->getIsExtant()){
                tn << nodes[i]->getIndex();
                std::string name = "S" + tn.str();
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

void SymbiontTree::recTipNamer(std::shared_ptr<Node> p, unsigned &nodeIndx, unsigned &tipIndx){
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
                std::string name = "S" + tn.str();
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

void SymbiontTree::setBranchLengths(){
    double bl = NAN;
    for(auto node : nodes){
        bl = node->getDeathTime() - node->getBirthTime();
        branchLengths.push_back(std::move(bl));
        node->setBranchLength(bl);
    }
}

void SymbiontTree::setPresentTime(double currentT){
    for(auto extantNode : extantNodes){
        extantNode->setDeathTime(currentT);
        extantNode->setIsExtant(true);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}


void SymbiontTree::updateCurrentMap(unsigned int oldHostIndx, unsigned int newHostIndx){
    auto found = symbHostMap.find(oldHostIndx);
    if (found != symbHostMap.end()) {
        // Swap value from oldKey to newKey, note that a default constructed value
        // is created by operator[] if 'm' does not contain newKey.
        std::swap(symbHostMap[newHostIndx], found->second);
        // Erase old key-value from map
        symbHostMap.erase(found);
    }

}



void SymbiontTree::cospeciationMapUpdate(unsigned int oldHostIndx,
                                         unsigned int numNodesHost,
                                         unsigned int oldSymbIndx){

    std::vector<unsigned int> symbsOnHost = symbHostMap[oldHostIndx];
    std::vector<unsigned int> leftHostSymbiontsValues;
    std::vector<unsigned int> rightHostSymbiontsValues;
    std::vector<std::shared_ptr<Node>> nodesForUpdating = this->getNodes();
    for(unsigned i = 0; i < symbsOnHost.size(); i++){
        std::vector<unsigned int> hostsInSymb = nodesForUpdating[symbsOnHost[i]]->getHosts();
        if(oldSymbIndx == symbsOnHost[i]){
            leftHostSymbiontsValues.push_back(this->getNodesSize() - 1);
            rightHostSymbiontsValues.push_back(this->getNodesSize() - 2);
        }
        else{
            double which = unif_rand();
            if(which < 0.5){
                leftHostSymbiontsValues.push_back(symbsOnHost[i]);
                for(unsigned int i=0; i < hostsInSymb.size(); i++){
                    if(hostsInSymb[i] == oldHostIndx)
                        hostsInSymb[i] = numNodesHost - 1;
                }
            }
            else{
                rightHostSymbiontsValues.push_back(symbsOnHost[i]);
                for(unsigned int i=0; i < hostsInSymb.size(); i++){
                    if(hostsInSymb[i] == oldHostIndx)
                        hostsInSymb[i] = numNodesHost - 2;
                }
            }

        }
        nodesForUpdating[symbsOnHost[i]]->setHosts(hostsInSymb);
    }
    symbHostMap[numNodesHost - 2] = rightHostSymbiontsValues;
    symbHostMap[numNodesHost - 1] = leftHostSymbiontsValues;
    // need to update nodes[i].hosts
    auto it = symbHostMap.find(oldHostIndx);
    symbHostMap.erase(it);


}

void SymbiontTree::updateHostsInNodes(){

}

int SymbiontTree::getExtantIndxFromNodes(unsigned int nodesIndx){
    int count = 0;
    for(auto extantNode : extantNodes){
        if((unsigned) extantNode->getIndex() == nodesIndx)
            break;
        count++;
    }
    return count;
}
