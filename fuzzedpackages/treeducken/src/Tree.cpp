#include "Tree.h"
#include <vector>
#include <string>
#include <cmath>

using namespace Rcpp;

Node::Node()
{
    ldes = nullptr;
    rdes = nullptr;
    anc = nullptr;
    sib = nullptr;
    indx = -1;
    Lindx = -1;
    flag = -1;
    isRoot = false;
    isTip = false;
    isExtant = false;
    isDuplication = false;
    isExtinct = false;
    branchLength = 0.0;
    birthTime = 0.0;
    deathTime = 0.0;

}

Node::~Node(){

}

Tree::Tree(unsigned numExta, double curTime){
    numNodes = 0;
    // intialize tree with root
    root = std::shared_ptr<Node>(new Node());
    root->setAsRoot(true);
    root->setBirthTime(0.0);
    root->setIndx(0);
    root->setIsExtant(true);
    
    nodes.push_back(root);
    extantNodes.push_back(root);
    numExtant = 1;
    numTaxa = numExta;
    numExtinct = 0;
    currentTime = curTime;
    numTotalTips = 0;
}

Tree::Tree(unsigned numTax){
    numTaxa = numTax; 
    numNodes = 2 * numTax - 1;
    // intialize tree with root
    // root = new Node();
    // root->setAsRoot(true);
    // root->setBirthTime(0.0);
    // root->setIndx(0);
    // root->setIsExtant(true);
    // nodes.push_back(root);
    // extantNodes.push_back(root);
    root = nullptr;
    //numExtant = 1;
    //currentTime = 0.0;
}
// Implicit converter from R tree object (ala APE) into C++ tree class
Tree::Tree(SEXP rtree){
    Rcpp::List tr(rtree);
    Rcpp::NumericMatrix edge_mat = tr["edge"];
    std::vector<double> edge_lengths = tr["edge.length"];
    std::vector<std::string> tip_names = tr["tip.label"];
    double root_edge = tr["root.edge"];
    numNodes = tr["Nnode"];

    std::map<int,int> indMap;
    numTaxa = (int) tip_names.size();
    nodes.resize(numNodes + numTaxa);
    extantNodes.resize(numTaxa);
    branchLengths.resize(numNodes + numTaxa);

    root = std::shared_ptr<Node>(new Node());
    root->setAsRoot(true);
    root->setBirthTime(0.0);
    root->setBranchLength(root_edge);
    root->setDeathTime(root_edge - 0.0);
    root->setIsExtant(false);
    root->setIndx(numTaxa + 1);
    nodes[0] = root;
    indMap[numTaxa] = 0;
    branchLengths[0] = root_edge;

    int i = 0;
    std::vector<int> nodeIndices;
    while(i < numTaxa + numNodes - 1){
        NumericMatrix::Row edgeMatRow = edge_mat(i,_);
        int indx1 = edgeMatRow[0] - 1;
        int indx2 = edgeMatRow[1] - 1;
        std::shared_ptr<Node> p = std::shared_ptr<Node>(new Node());
        p->setBranchLength(edge_lengths[i]);
        branchLengths[i + 1] = edge_lengths[i];

        p->setIndx(indx2 + 1);
        indMap[indx2] = i + 1;
        p->setAnc(nodes[indMap[indx1]]);
        p->setBirthTime(nodes[indMap[indx1]]->getDeathTime());
        p->setDeathTime(edge_lengths[i] + nodes[indMap[indx1]]->getDeathTime());
        p->setBranchLength(p->getDeathTime() - p->getBirthTime());
        nodes[i + 1] = p;
        if(indx2 < numTaxa){
            nodes[indMap[indx2]]->setName(tip_names[indx2]);
            nodes[indMap[indx2]]->setIsTip(true);
            std::string tipName = nodes[indMap[indx2]]->getName();
            if(tipName.find("X") == 0)
            {
                nodes[indMap[indx2]]->setIsExtinct(true);
            }
            else
                nodes[indMap[indx2]]->setIsExtant(true);
            if(nodes[indMap[indx1]]->getLdes()){
                nodes[indMap[indx1]]->setRdes(nodes[indMap[indx2]]);
            }
            else{
                nodes[indMap[indx1]]->setLdes(nodes[indMap[indx2]]);
            }
        }
        else{
            if(nodes[indMap[indx1]]->getLdes()){
                nodes[indMap[indx1]]->setRdes(nodes[indMap[indx2]]);
            }
            else{
                nodes[indMap[indx1]]->setLdes(nodes[indMap[indx2]]);
            }
        }
        i++;
    }
    this->setTipsFromRtree();
}




Tree::~Tree(){

    clearNodes(root);
    extantNodes.clear();
    nodes.clear();
}


void Tree::setTipsFromRtree(){
    int extTipCount = 0;
    numExtant = 0;
    numExtinct = 0;
    for(auto node : nodes){
        if(node->getIsTip()){
            if(node->getIsExtinct()){
                numExtinct++;
            }
            else{
                extantNodes[extTipCount] = node;
                numExtant++;
                extTipCount++;
            }
        }
        if(extTipCount == numTaxa)
            break;
    }

}


double Tree::findMaxNodeHeight(){
    std::shared_ptr<Node> p = root;
    double brlen = p->getBranchLength();
    while(p->getLdes()){
        p = p->getLdes();
        brlen += p->getBranchLength();
    }
    return brlen;
}

void Tree::clearNodes(std::shared_ptr<Node> currNode){
    if(currNode == nullptr){
        return;
    }

    clearNodes(currNode->getRdes());
    clearNodes(currNode->getLdes());
    currNode = nullptr;

}

void Tree::zeroAllFlags(){
    for(auto node : nodes) {
        node->setFlag(0);
    }
}

void Tree::setWholeTreeFlags(){
    this->zeroAllFlags();
    for(auto node : nodes){
        if(node->getIsTip()){
            node->setFlag(1);
        }
    }
    setSampleFromFlags();
}


void Tree::setExtantTreeFlags(){
    this->zeroAllFlags();
    for(auto node : nodes){
        if(node->getIsExtant())
            node->setFlag(1);
    }

    this->setSampleFromFlags();
}


void Tree::setSampleFromFlags(){
    int flag = -1;
    std::shared_ptr<Node> q = nullptr;
    for(auto node : nodes) {
        if(node->getIsTip()) {
            flag = node->getFlag();
            q = node;
            if(flag == 1){
                do{
                    q = q->getAnc();
                    flag = q->getFlag();
                    flag++;
                    q->setFlag(flag);
                }while (q->getIsRoot() == false && flag < 2);
            }
        }
    }
}


double Tree::getTotalTreeLength(){
    double sum = 0.0;
    for(auto node : nodes){
        std::shared_ptr<Node> n = node;
        sum += n->getBranchLength();
    }
    return sum;
}

double Tree::getTreeDepth(){
    double td = 0.0;
    std::shared_ptr<Node> r = this->getRoot();
    while(r->getIsTip() == false){
        if(!(r->getLdes()->getIsExtinct()))
            r = r->getLdes();
        else
            r = r->getRdes();
    }
    while(r->getIsRoot() == false){
        td += r->getBranchLength();
        r = r->getAnc();
    }
    td += r->getBranchLength();
    return td;
}

void Tree::reconstructTreeFromSim(std::shared_ptr<Node> oRoot){
    std::shared_ptr<Node> n = std::shared_ptr<Node>(new Node());
    unsigned tipCounter = numExtant;
    unsigned intNodeCounter = 0;
    reconstructLineageFromSim(n, oRoot, tipCounter, intNodeCounter);
}

void Tree::reconstructLineageFromSim(std::shared_ptr<Node> currN,
                                    std::shared_ptr<Node> prevN, 
                                    unsigned &tipCounter, 
                                    unsigned &intNodeCounter) {
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

        p = std::shared_ptr<Node>(new Node());
        tipCounter++;
        p->setIndx(tipCounter);
        p->setBranchLength(brlen);
        p->setIsTip(true);
        p->setName(prevN->getName());
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
                reconstructLineageFromSim(s1, prevN->getLdes(), tipCounter, intNodeCounter);
            if(prevN->getRdes()->getFlag() > 0)
                reconstructLineageFromSim(s1, prevN->getRdes(), tipCounter, intNodeCounter);


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
                        stop("ERROR: Problem adding a tip to the tree!");
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
                reconstructLineageFromSim(currN, prevN->getLdes(), tipCounter, intNodeCounter);
            else
                reconstructLineageFromSim(currN, prevN->getRdes(), tipCounter, intNodeCounter);
        }
    }
}
// Gene tree version only
void Tree::getRootFromFlags(bool isGeneTree){
    std::shared_ptr<Node> p = nullptr;

    this->setExtantTreeFlags();
    int numNodes = nodes.size() - 1;
    if(isGeneTree){
        for(int i=numNodes; i > 0; i--){
            p = nodes[i];
            if(p->getFlag() >= 2){
                extantRoot = p;
                p->setAsRoot(true);
                break;
            }

        }
    }
    else{

        for(int i=0; i < numNodes; i++){
            p = nodes[i];
            if(p->getFlag() >= 2){
                extantRoot = p;
                p->setAsRoot(true);
                break;
            }
        }
    }
}



double Tree::getEndTime(){
    double tipDtime = 0.0;
    for(auto node : nodes){
        if(node->getIsTip() && node->getIsExtant()){
            tipDtime = node->getDeathTime();
            break;
        }
    }

    return tipDtime;
}

void Tree::scaleTree(double scVal){
	double scaler = scVal;
	for(auto node : nodes){
		double bt = node->getBirthTime();
        double dt = node->getDeathTime();

		node->setBirthTime(bt * scaler);
        node->setDeathTime(dt * scaler);
        node->setBranchLength(node->getDeathTime() - node->getBirthTime());
	}
    currentTime = getTreeDepth();
}

void Tree::scaleTreeDepthToValue(double scVal){
	
	double depth = getTreeDepth();
	double scaler = scVal / depth;

	for(auto node : nodes){
		double bt = node->getBirthTime();
        double dt = node->getDeathTime();

		node->setBirthTime(bt * scaler);
        node->setDeathTime(dt * scaler);
        node->setBranchLength(node->getDeathTime() - node->getBirthTime());
	}
}

int Tree::calculatePatristicDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
    int count = 0;
    if(n1 != n2){
        while(n1->getIndex() != n2->getIndex()){
            count++;
            n1 = n1->getAnc();
            n2 = n2->getAnc();
        }
    }
    return count;
}


std::vector<std::string> Tree::getTipNames(){
    std::vector<std::string> tipNames;

    for(auto node : nodes){
        if(node->getIsTip())
            tipNames.push_back(node->getName());
    }
    return tipNames;
}

std::vector<std::string> Tree::getNodeLabels()
{
    std::vector<std::string> nodeLabels;
    for(auto node : nodes)
    {
        if(!(node->getIsTip())){
            if(node->getIsDuplication())
                nodeLabels.push_back(node->getName());
            else
                nodeLabels.push_back("");
        }
    }
    return nodeLabels;
}

void Tree::setNumExtant(){
    numExtant = 0;
    for(auto node : nodes){
        if(node->getIsTip() && node->getIsExtant())
            numExtant++;
    }
}

void Tree::setNumExtinct(){
    numExtinct = 0;
    for(auto node : nodes){
        if(node->getIsTip() && node->getIsExtinct())
            numExtinct++;
    }
}

// TODO: write a function to convert bdsa sims to format to read out into R

void Tree::reindexForR(){
    unsigned int intNodeCount = numExtant + numExtinct + 1;
    unsigned int tipCount = 1;
    for(unsigned int i = 0; i < nodes.size(); i++){
        if(nodes[i]->getIsTip()){
            nodes[i]->setIndx(tipCount);
            tipCount++;
        }
        else if(nodes[i]->getIsRoot()){
            // do nothing
        }
        else{
            nodes[i]->setIndx(intNodeCount);
            intNodeCount++;
        }
    }
}
// remember if something breaks you edited the notoriously
//  sketchy reconstruct lineages WTD

NumericMatrix Tree::getEdges(){
    int numRows = (int) nodes.size() - 1;
    NumericMatrix edgeMat(numRows, 2);
    for(unsigned int i=1; i < nodes.size(); i++){
        if(!(nodes[i]->getIsRoot())){

            NumericMatrix::Row row = edgeMat(i - 1, _);

            row[0] = nodes[i]->getAnc()->getIndex();
            row[1] = nodes[i]->getIndex();
        }
    }
    return edgeMat;
}



std::vector<double> Tree::getEdgeLengths(){
    std::vector<double> edgeLengths;
    edgeLengths = branchLengths;
    edgeLengths.erase(edgeLengths.begin());
    return edgeLengths;
}

void Tree::switchIndicesFirstToSecond(std::map<int,int> mappy){
    for(unsigned int i = 0; i < nodes.size(); i++){
        int newIndx = mappy[nodes[i]->getIndex()];
        nodes[i]->setIndx(newIndx);
    }
}

