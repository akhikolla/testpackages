#include "Simulator.h"
#include <iostream>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
Simulator::Simulator(unsigned nt, double lambda, double mu, double rho)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 1;
    currentSimTime = 0.0;
    numTaxaToSim = nt;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;

    treeScale = -1;
    numLoci = 0;
    numGenes = 0;
    geneBirthRate = 0.0;
    geneDeathRate = 0.0;
    transferRate = 0.0;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0.0;

}


Simulator::Simulator(unsigned ntax,
                     double lambda,
                     double mu,
                     double rho,
                     unsigned numLociToSim,
                     double gbr,
                     double gdr,
                     double lgtr,
                     std::string transfType)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 2;
    currentSimTime = 0.0;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;

    numLoci = numLociToSim;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0;
    transferType = transfType;

}

Simulator::Simulator(unsigned ntax,
                     double lambda,
                     double mu,
                     double rho,
                     unsigned numLociToSim,
                     double gbr,
                     double gdr,
                     double lgtr,
                     unsigned ipp,
                     double Ne,
                     double genTime,
                     int ng,
                     double og,
                     double ts,
                     bool sout)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 3;
    currentSimTime = 0.0;
    numTaxaToSim = ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    numLoci = numLociToSim;
    numGenes = ng;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = ipp;
    popSize = Ne;
    printSOUT = sout;
    generationTime = genTime;
    geneTrees.resize(ng);
    treeScale = ts;
}

Simulator::Simulator(double stopTime,
          double hostSpeciationRate,
          double hostExtinctionRate,
          double symbSpeciationRate,
          double symbExtinctionRate,
          double switchingRate,
          double csr,
          double rho,
          int hl){

    speciationRate = hostSpeciationRate;
    extinctionRate = hostExtinctionRate;
    samplingRate = rho;
    cospeciationRate = csr;
    geneBirthRate = symbSpeciationRate;
    geneDeathRate = symbExtinctionRate;
    transferRate = switchingRate;
    timeToSim = stopTime;

    hostLimit = hl;

    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    symbiontTree = nullptr;
    
}



Simulator::~Simulator(){
    gsaTrees.clear();
    for(auto locusTree : locusTrees){
        geneTrees.clear();
    }
    locusTrees.clear();
    assocMat.clear();
}

void Simulator::initializeSim(){
    spTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate));
}


void Simulator::initializeEventVector(){
  inOrderVecOfHostIndx.push_back(0);
  inOrderVecOfSymbIndx.push_back(0);
  inOrderVecOfEvent.push_back("I");
  inOrderVecOfEventTimes.push_back(0.0);
}
/*
 Main general sampling algorithm (GSA) for simulating a species tree to an expected number
 of tips under the constant-rate birth-death process.

 Below is the machinery to use GSA sampling (Hartmann 2010) to simulate a species tree.
 Much of this code is modified from FossilGen (written by Tracy Heath)
 */
bool Simulator::gsaBDSim(){
    double timeIntv = NAN;
    double sampTime = NAN;
    bool treeComplete = false;
    // make a species tree object with the number of taxa to sim to, currsimtime (0.0)
    // and speciation and extinction rate
    spTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate));
    double eventTime = NAN;
    // runs until the number of extant tips reaches gsaStop (set by users, default is 10*number to sim to)
    while(spTree->getNumExtant() < gsaStop){
        // find the time to the next event with rate speciation rate + extinction rate * number of current tips
        eventTime = spTree->getTimeToNextEvent();
        // add this to the sim time tracker
        currentSimTime += eventTime;
        // speciation or extinction occurs
        spTree->ermEvent(currentSimTime);
        if(spTree->getNumExtant() < 1){
            // if the tree goes to 0 tips prematurely end
            // return false for a non-tree
            treeComplete = false;
            return treeComplete;
        } // otherwise if the number of current tips is the number to sim to
        else if(spTree->getNumExtant() == numTaxaToSim){
            // get time to another event and add this to the end
            timeIntv = spTree->getTimeToNextEvent();
            // randomly choose some proportion of the above 'timeIntv' and add
            // to the sim time tracker
            sampTime = (unif_rand() * timeIntv) + currentSimTime;
            // set this as the present time for this sub-tree
            spTree->setPresentTime(sampTime);
            // reconstruct this tree from the root of the whole tree to the
            // number of extant tips (i.e. numTaxaToSim) and add to a vector
            // gsaTrees
            processGSASim();
        }

    }
    // randomly pick one of the gsaTrees
    unsigned gsaRandomTreeID = unif_rand() * (gsaTrees.size() - 1);
    spTree = gsaTrees[gsaRandomTreeID];
    // process this one
    processSpTreeSim();
    // set branch length variable in each Node of spTree.nodes
    spTree->setBranchLengths();
    // set tip names
    spTree->setTreeTipNames();
    // set currentSimTime
    currentSimTime = spTree->getCurrentTimeFromExtant();
    treeComplete = true;

    return treeComplete;
}

// prune sub trees from larger GSA tree to randomly choose from so
// that we are correctly sampling from the birth-death distribution
void Simulator::processGSASim(){
  // make a new species tree with numTaxaToSim + however many extinct tips there are
    auto tt = std::shared_ptr<SpeciesTree>(new SpeciesTree(numTaxaToSim + spTree->getNumExtinct()));
  // prep this by setting flags for nodes so that we know what is extant tip, extinct tip,
  // internode, and root
    this->prepGSATreeForReconstruction();
  // get our root from the original tree
   auto simRoot = spTree->getRoot();
  // set that node to by the root of our new tree
    tt->setRoot(simRoot);
  // reconstruct the tree recursively from the root
    tt->reconstructTreeFromGSASim(simRoot);
  // add to vector of SpeciesTree
    gsaTrees.push_back(tt);
}

// post simulation processing function since we overwrite the tree held in
// spTree
void Simulator::processSpTreeSim(){
  // set the speciation and extinction rate
    spTree->setSpeciationRate(speciationRate);
    spTree->setExtinctionRate(extinctionRate);
  // populate nodes vector
    spTree->popNodes();
  // set the number of extant and extinct tips respecitively
    spTree->setNumExtant();
    spTree->setNumExtinct();

}

// Wrapper for SpeciesTree::setGSATipTreeFlags
void Simulator::prepGSATreeForReconstruction(){
    spTree->setGSATipTreeFlags();
}

// Wrapper to make sure that the tree output by gsaBDSim is a proper tree with
// the correct number of tips and not just a tree with no tips
bool Simulator::simSpeciesTree(){
    bool good = false;
    while(!good){
        good = gsaBDSim();
    }
    return good;
}

// Wrapper to make sure that the tree output by bdSimpleSim is a proper tree with
//  tips and not just a tree with no tips
bool Simulator::simSpeciesTreeTime(){
  bool good = false;
  while(!good){
    good = bdSimpleSim();
  }
  return good;
}
// Simulation function for birth-death tree simulation to a set time (note that
// this can produce trees with no tips)
bool Simulator::bdSimpleSim(){
  bool treeComplete = false;
  // set time to 0.0
  currentSimTime = 0.0;
  // get the time to simulate to
  double stopTime = this->getTimeToSim();
  double eventTime = NAN;
  // make a variable SpeciesTree to hold our tree. numTaxaToSim is set to 1 here.
  // this is arbitrary and I should've planned out my classes and constructors better
  spTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate));
  while(currentSimTime < stopTime){
    // get the time to the next event as a function of speciation and extinction rates and number of currently alive
    // tips
    eventTime = spTree->getTimeToNextEvent();
    // add this to the currentSimTime
    currentSimTime += eventTime;
    // check if that went over stop time if it did set it to that time
    if(currentSimTime >= stopTime){
      currentSimTime = stopTime;
    }
    else{ // otherwise choose an event at random
      spTree->ermEvent(currentSimTime);
    }
    // if tree goes to 0 living tips end prematurely
    if(spTree->getNumExtant() < 1){
      treeComplete = false;
      return treeComplete;
    }
  }

  if(spTree->getNumExtant() <= 1){
    treeComplete = false;
    return treeComplete;
  }

  treeComplete = true;
  // set currentSimTime to stopTime in case it hasn't been
  currentSimTime = stopTime;

  // set the tree to the end time
  spTree->setPresentTime(currentSimTime);
  return treeComplete;
}

// wrapper for the paired birth-death process
bool Simulator::simHostSymbSpeciesTreePair(){
  bool good = false;
  while(!good){
    good = pairedBDPSim();

  }
  return good;
}

// Function simulate host and symbiont tree at the same time to a set time
bool Simulator::pairedBDPSim(){
  bool treePairGood = false;

  currentSimTime = 0.0;
  // set stopTime
  double stopTime = this->getTimeToSim();
  // make a SpeciesTree (this is the host tree)
  spTree = std::shared_ptr<SpeciesTree>(new SpeciesTree(1, currentSimTime, speciationRate, extinctionRate));

  // and a SymbiontTree (this is the symbiont tree)
  symbiontTree = std::shared_ptr<SymbiontTree>( new SymbiontTree(1,
                                                                currentSimTime,
                                                                geneBirthRate,
                                                                geneDeathRate,
                                                                transferRate,
                                                                hostLimit));

  double eventTime = NAN;
  // initialize the four vectors that are output in R as the event dataframe
  this->initializeEventVector();
  // set the association matrix to start with the host and symbiont being associated
  // a 1x1 matrix of 1
  assocMat = arma::ones<arma::umat>(1,1);
  while(currentSimTime < stopTime){
    // get time to the the next joint event based on the speciation rate, extinction rate,
    // and cospeciation rate
    eventTime = symbiontTree->getTimeToNextJointEvent(speciationRate,
                                                 extinctionRate,
                                                 cospeciationRate,
                                                 assocMat);
    currentSimTime += eventTime;
    // if we exceed the sim time set to stopTime so as not to go over
    if(currentSimTime >= stopTime){
      currentSimTime = stopTime;
    }
    else{
      // otherwise a cophylogenetic event occurs, this can be three things:
      // host event (host speciation or extinction)
      // symbiont event (symbiont speciation or extinction)
      // or a joint event (a.k.a. a cospeciation)
      // this returns the association matrix
      assocMat = this->cophyloEvent(currentSimTime, assocMat);
    }
    // if either tree goes to 0 or the association matrix becomes malformed
    // prematurely end the simulation, clearing the event dataframe vectors
    if(spTree->getNumExtant() < 1 ||
       symbiontTree->getNumExtant() < 1 ||
       assocMat.n_rows < 1 ||
       assocMat.n_cols < 1){
      treePairGood = false;
      this->clearEventDFVecs();
      return treePairGood;
    }
  }
  // TODO: not sure if this is needed
  if(spTree->getNumExtant() <= 1 || symbiontTree->getNumExtant() <= 1){
    treePairGood = false;
    this->clearEventDFVecs();
    return treePairGood;
  }
  treePairGood = true;
  currentSimTime = stopTime;
  // set the present time in both host and symbiont tree
  symbiontTree->setPresentTime(currentSimTime);
  spTree->setPresentTime(currentSimTime);

  return treePairGood;
}

// cophyloEvent - chooses which event occurs based on the rates of the 6 different events
// note that this can likely be simplified mathematically
arma::umat Simulator::cophyloEvent(double eventTime, arma::umat assocMat){
  double hostEvent = speciationRate + extinctionRate;
  double symbEvent = geneBirthRate + geneDeathRate + transferRate;
  double cospecEvent = cospeciationRate;
  // which tree then ermEvent whichever or cospeciation event
  // probability of a host event
  double hostEventProb = hostEvent / (hostEvent + symbEvent + cospecEvent);
  // probability of symb event
  double symbEventProb = symbEvent / (hostEvent + symbEvent + cospecEvent);
  symbEventProb += hostEventProb;
  double whichEvent = unif_rand();
  // randomly chooose host, symb, or cospeciation event
  if(whichEvent < hostEventProb){
    assocMat = this->cophyloERMEvent(eventTime, assocMat);
  }
  else if(whichEvent < symbEventProb){
    assocMat = this->symbiontTreeEvent(eventTime, assocMat);
  }
  else{

    assocMat = this->cospeciationEvent(eventTime, assocMat);
  }
  return assocMat;
}

// Function that creates the dataframe out of the vectors that record events
Rcpp::DataFrame Simulator::createEventDF(){
  this->updateEventIndices();
  DataFrame df = DataFrame::create(Named("Symbiont Index") = inOrderVecOfSymbIndx,
                                   Named("Host_Index") = inOrderVecOfHostIndx,
                                   Named("Event_Type") = inOrderVecOfEvent,
                                   Named("Event_Time") = inOrderVecOfEventTimes);
  return df;
}

// Function that clears the vectors that record events
void Simulator::clearEventDFVecs(){
  inOrderVecOfHostIndx.erase(inOrderVecOfHostIndx.begin(),
                             inOrderVecOfHostIndx.end());
  inOrderVecOfSymbIndx.erase(inOrderVecOfSymbIndx.begin(),
                             inOrderVecOfSymbIndx.end());
  inOrderVecOfEvent.erase(inOrderVecOfEvent.begin(),
                          inOrderVecOfEvent.end());
  inOrderVecOfEventTimes.erase(inOrderVecOfEventTimes.begin(),
                               inOrderVecOfEventTimes.end());
}


// update the indices of the event vector from the C++ indexing where the root is
// index 0, to the APE package indexing where the root is numTips+1
void Simulator::updateEventIndices(){

  for(int i = 0; i < inOrderVecOfHostIndx.size(); i++){
    int oldHostIndx = inOrderVecOfHostIndx(i);
    int oldSymbIndx = inOrderVecOfSymbIndx(i);
    int newHostIndx = spTree->getIndexFromNodes(oldHostIndx);

    int newSymbIndx = symbiontTree->getIndexFromNodes(oldSymbIndx);
    inOrderVecOfHostIndx(i) = newHostIndx;
    inOrderVecOfSymbIndx(i) = newSymbIndx;
  }
}

// add an element to each event vector for an event
void Simulator::updateEventVector(int h, int s, int e, double time){
  inOrderVecOfHostIndx.push_back(h);
  inOrderVecOfSymbIndx.push_back(s);
  switch(e) {
    case 0:
      // symbiont loss
      inOrderVecOfEvent.push_back(std::move("SL"));
      break;
    case 1:
      // host loss
      inOrderVecOfEvent.push_back(std::move("HL"));
      break;
    case 2:
      // symbiont gain
      inOrderVecOfEvent.push_back(std::move("SG"));
      break;
    case 3:
      // host gain
      inOrderVecOfEvent.push_back("HG");
      break;
    case 4:
      // association gain
      inOrderVecOfEvent.push_back("AG");
      break;
    case 5:
      // association loss
      inOrderVecOfEvent.push_back("AL");
      break;
    case 6:
      // cospeciation
      inOrderVecOfEvent.push_back("C");
      break;
    default:
      Rcout << "not sure what happened there folks." << std::endl;
  }
  inOrderVecOfEventTimes.push_back(time);
}



// Event occurring on the symbiont tree at eventTiem with matrix assocMat
arma::umat Simulator::symbiontTreeEvent(double eventTime, arma::umat assocMat){
  // get the number of tips on the symbiont tree
  unsigned int numExtantSymbs = symbiontTree->getNumExtant();
  // randomly choose one of these to have an event on
  // arma::uword nodeInd = unif_rand()*(numExtantSymbs - 1);
  arma::uword nodeInd = 0;
  if(numExtantSymbs > 1)
    nodeInd = arma::randi<arma::uword>(arma::distr_param(0, numExtantSymbs - 1));

  // relative birth rate, uses geneBirthRate, geneDeathRate and transferRate
  // for the only purpose so that I did not need to add more members to the class
  // geneBirthRate = symbiont speciation rate
  // geneDeathRate = symbiont extinction rate
  // transferRate = symbiont host expansion rate
  double relBr = geneBirthRate / (geneBirthRate
                                    + geneDeathRate
                                    + transferRate);
  // relative death rate
  double relDr = relBr + (geneDeathRate / (geneBirthRate
                                            + geneDeathRate
                                            + transferRate));
  double decid = unif_rand();
  // make sure our times are correctly set
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);

  unsigned int numExtantHosts = spTree->getNumExtant();
  
  arma::urowvec rvec = assocMat.row(nodeInd);

  // delete the appropriate row of the association matrix
  assocMat.shed_row(nodeInd);
  // randomly decide between birth, death, and transfer
  if(decid < relBr){
    // update the event vectors
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      2,
                      eventTime);
    // birth event on the symbiont tree
    symbiontTree->lineageBirthEvent(nodeInd);
    numExtantSymbs = symbiontTree->getNumExtant();


    assocMat.resize(numExtantSymbs, numExtantHosts);
    // add two new rows to the association matrix that are identical to the deleted row
    assocMat(numExtantSymbs-2, arma::span::all) = rvec;
    assocMat(numExtantSymbs-1, arma::span::all) = rvec;

    // sort symbs on new hosts
    for(arma::uword i = 0; i < rvec.n_cols; i++){
      if(rvec(i) == 1){
        // update the event vectors for the sorting events
        updateEventVector(spTree->getNodesIndxFromExtantIndx((int) i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 2),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx((int) i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                          4,
                          eventTime);
      }
    }
  }
  else if(decid < relDr){
    // update the event vectors for the main event
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      0,
                      eventTime);
    numExtantSymbs = symbiontTree->getNumExtant();
    for(arma::uword i = 0; i< rvec.n_cols; i++){
      if(rvec(i) == 1){
        // update the event vectors for the sorting events
        updateEventVector(spTree->getNodesIndxFromExtantIndx((int) i),
                              symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                              4,
                              eventTime);
      }
    }

    // death event
    symbiontTree->lineageDeathEvent(nodeInd);
  }
  else{
    // expansion event (a.k.a. birth event with the addition of one host in a descendent symbiont lineage)
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      2,
                      eventTime);
    // check if the sum of the row equals the number of columns
    // in other words is the symbiont at rvec occupying all the hosts already?
    // if no, if yes go to else
    if(sum(rvec) != numExtantHosts){
      //std::vector<arma::uword> hostIndices;
      arma::urowvec hostIndices(rvec.n_cols, arma::fill::zeros);

      // make a list of unoccupied hosts
      for(arma::uword i = 0; i < rvec.n_cols; i++){
        if(rvec(i) < 1)
          hostIndices(i) = 1;
      }
      // randomly choose from one of those unoccupied hosts
      arma::uvec unoccupiedHosts = arma::find(hostIndices > 0);
      //  arma::uvec hostsWithSymbs = arma::find(hostIndices > 0);
      arma::uword hostEndpoint = unoccupiedHosts.n_elem - 1;
      arma::uword hostInd = 0;
      if(unoccupiedHosts.n_elem > 1)
        hostInd = arma::randi<arma::uword>(arma::distr_param(0, hostEndpoint)); //col of assocMat

      // birth event
      symbiontTree->lineageBirthEvent(nodeInd);
      numExtantSymbs = symbiontTree->getNumExtant();
      // add two rows
      assocMat.resize(numExtantSymbs, numExtantHosts);
      // make one of these rows the same as the deleted row
      assocMat(numExtantSymbs-2, arma::span::all) = rvec;
      // change that row to have an extra one where the randomly picked unoccupied host was
      rvec(hostIndices(hostInd)) = 1;
      // make that a new row
      assocMat(numExtantSymbs-1, arma::span::all) = rvec;

      // sort symbs on new hosts
      for(arma::uword i = 0; i < rvec.n_cols; i++){
        if(rvec(i) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx((int) i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                            4,
                            eventTime);
            updateEventVector(spTree->getNodesIndxFromExtantIndx((int) i),
                              symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                              4,
                              eventTime);
        }
      }
    }
    else{ // if all hosts are occupied this is just a regular birth event
      symbiontTree->lineageBirthEvent(nodeInd);
      numExtantSymbs = symbiontTree->getNumExtant();


      assocMat.resize(numExtantSymbs, numExtantHosts);

      assocMat(numExtantSymbs-2, arma::span::all) = rvec;
      assocMat(numExtantSymbs-1, arma::span::all) = rvec;

      // sort symbs on new hosts
      for(arma::uword i = 0; i < rvec.n_cols; i++){
        if(rvec(i) == 1){
          // shuffle might be better here?
          updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                            4,
                            eventTime);
        }
      }

    }
  }
  return assocMat;
}


// cophylogenetic erm event is actually the function for the host event
// apologies for the misleading name
arma::umat Simulator::cophyloERMEvent(double eventTime, arma::umat assocMat){
  unsigned numExtantHosts = spTree->getNumExtant();
  // randomly pick a host
  arma::uword nodeInd = 0;
  if(numExtantHosts > 1)
    nodeInd = arma::randi<arma::uword>(arma::distr_param(0, numExtantHosts - 1));
  // choose event based on the relative birth rate
  double relBr = speciationRate / (speciationRate + extinctionRate);
  bool isBirth = (unif_rand() < relBr ? true : false);
  // set the times to keep up
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  unsigned numExtantSymbs = symbiontTree->getNumExtant();
  // delete a column 
  arma::ucolvec cvec = assocMat.col(nodeInd);

  assocMat.shed_col(nodeInd);

  if(isBirth){
    // add the birth event to event vectors
    updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      3,
                      eventTime);
    // birth event occur
    spTree->lineageBirthEvent(nodeInd);
    // recalculate num extant hosts
    numExtantHosts = spTree->getNumExtant();

    // add two rows
    assocMat.resize(numExtantSymbs, numExtantHosts);
    // make two new rows of data frame to be clear about which host speciated into what
    updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      4,
                      eventTime);
    updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      4,
                      eventTime);

   // sort symbs on new hosts
   for(arma::uword i = 0; i < cvec.n_rows; ++i) {
      if(cvec(i) == 1){
        arma::umat rr = arma::randi<arma::umat>(1,2, arma::distr_param(0,1));
        if(rr(0,0) == 0 && rr(0,1) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            5,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

        }
        else if(rr(0,0) == 1 && rr(0,1) == 0){

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

            updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            5,
                            eventTime);

        }
        else if(rr(0,0) == 1 && rr(0,1) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);
        }
        else{
          rr.replace(0,1);
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

        }

        assocMat(i, arma::span(numExtantHosts-2, numExtantHosts-1)) = rr;

      }
      else{
        arma::umat rr(1,2,arma::fill::zeros);// = arma::zeros<arma::umat>(1,2);

        assocMat(i, arma::span(numExtantHosts-2, numExtantHosts-1)) = rr;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);

      }
    }
  }
  else{ // otherwise death occurs
    // update event vectors
    updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      1,
                      eventTime);
    numExtantSymbs = symbiontTree->getNumExtant();
    // check rows for 0's, rows with 0s get deleted in symbiont tree
    arma::uword numRows = assocMat.n_rows;
    arma::uvec hostless(numExtantSymbs, arma::fill::zeros);// = arma::zeros<arma::uvec>(numExtantSymbs);
    for(arma::uword i = numRows; i != 0; i--){
      if(!(any(assocMat.row(i-1)))){
        updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                          symbiontTree->getNodesIndxFromExtantIndx(i-1),
                          0,
                          eventTime);
        symbiontTree->lineageDeathEvent(i-1);

        hostless(i-1) = 1;
      }
    }
    // host tree death event
    spTree->lineageDeathEvent(nodeInd);
    // delete the rows from the association matrix
    arma::uvec toBeDeleted = arma::find(hostless);
    assocMat.shed_rows(toBeDeleted);
  }
  return assocMat;
}

// Cospeciation event occur
arma::umat Simulator::cospeciationEvent(double eventTime, arma::umat assocMat){
  // draw index of host
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  unsigned numExtantHosts = spTree->getNumExtant();
  arma::urowvec hostIndices(assocMat.n_cols, arma::fill::zeros);
  // pick a host with symbionts at random
  for(arma::uword i = 0; i < assocMat.n_cols; i++){
    if(sum(assocMat.col(i)) > 0)
      hostIndices(i) = 1;
  }

  arma::uvec hostsWithSymbs = arma::find(hostIndices > 0);
  arma::uword hostEndpoint = hostsWithSymbs.n_elem - 1;
  arma::uword indxOfHost = 0;
  if(hostsWithSymbs.n_elem > 1)
    indxOfHost = arma::randi<arma::uword>(arma::distr_param(0, hostEndpoint)); //col of assocMat

  
  arma::ucolvec cvec = assocMat.col(hostsWithSymbs(indxOfHost));

  arma::uvec symbIndices = arma::find(cvec);
  arma::uword symbEndpoint = symbIndices.n_elem - 1;

  arma::uword indxOfSymb = 0;
  if(symbIndices.n_elem > 1)
    indxOfSymb = arma::randi<arma::uword>(arma::distr_param(0, symbEndpoint));
  arma::urowvec rvec = assocMat.row(symbIndices(indxOfSymb));
  // add a C to the event vectors
  updateEventVector(spTree->getNodesIndxFromExtantIndx(hostsWithSymbs(indxOfHost)),
                    symbiontTree->getNodesIndxFromExtantIndx(symbIndices(indxOfSymb)),
                    6,
                    eventTime);
  // birth in both trees at the same time
  spTree->lineageBirthEvent(hostsWithSymbs(indxOfHost));
  symbiontTree->lineageBirthEvent(symbIndices(indxOfSymb));

  numExtantHosts = spTree->getNumExtant();
  unsigned numExtantSymbs = symbiontTree->getNumExtant();
  // delete a row and column of assocaition matrix

  assocMat.shed_row(symbIndices(indxOfSymb));
  if(!(assocMat.is_empty()))
    assocMat.shed_col(hostsWithSymbs(indxOfHost));


  cvec.shed_row(symbIndices(indxOfSymb));

  rvec.shed_col(hostsWithSymbs(indxOfHost));
  // add 2 new rows and 2 new cols
  assocMat.resize(numExtantSymbs, numExtantHosts);
  // throw the identity matrix into the bottom right corner representing
  // that each new host is associated with one new symbiont

  assocMat(arma::span(numExtantSymbs - 2, numExtantSymbs - 1), 
           arma::span(numExtantHosts - 2, numExtantHosts -1)) = arma::eye<arma::umat>(2,2);

  // record in event vector
  updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                    symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                    4,
                    eventTime);

  updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                    symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                    4,
                    eventTime);
  // loop through cvec to sort the old hosts of the ancestor symbiont on new symbionts
  for(arma::uword i = 0; i < cvec.n_rows; i++){
    if(cvec(i) == 1){

      arma::umat rr(1,2,arma::fill::ones);// = arma::ones<arma::umat>(1,2);
      int randOne = unif_rand() * 2;
      if(randOne == 0){
        rr(0, 0) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          4,
                          eventTime);
      }
      else{
        rr(0, 1) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);
      }
      assocMat(i, arma::span(numExtantHosts - 2, numExtantHosts - 1)) = rr;
      // assocMat.submat(i, numExtantHosts-2, i, numExtantHosts-1) = rr;
    }
    else{
      arma::umat rr(1,2,arma::fill::zeros);// = arma::zeros<arma::umat>(1,2);
      assocMat(i, arma::span(numExtantHosts - 2, numExtantHosts - 1)) = rr;
      updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                        symbiontTree->getNodesIndxFromExtantIndx(i),
                        5,
                        eventTime);
      updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                        symbiontTree->getNodesIndxFromExtantIndx(i),
                        5,
                        eventTime);
    }
  }



  // loop through rvec to sort the old symbs of the ancestor host on new hosts
  for(arma::uword i = 0; i < rvec.n_elem; i++){
    if(rvec(i) == 1){
      arma::umat cr(2,1,arma::fill::ones); //= arma::ones<arma::umat>(1,2);

      int randOne = unif_rand() * 2;
      if(randOne == 0){
        cr(0,0) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                          5,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                          4,
                          eventTime);
      }
      else{
        cr(1, 0) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                          5,
                          eventTime);
      }

      assocMat(arma::span(numExtantSymbs - 2, numExtantSymbs - 1), i) = cr;
      // assocMat.submat(numExtantSymbs-2, i, numExtantSymbs-1,i) = rr;
    }
    else{
      arma::umat cr(2,1,arma::fill::zeros);// = arma::zeros<arma::umat>(1,2);
      assocMat(arma::span(numExtantSymbs - 2, numExtantSymbs - 1), i) = cr;
    }
  }
  return assocMat;
}

// locus tree simulation function, probably should be renamed
bool Simulator::bdsaBDSim(){
    bool treesComplete = false;
    // get the stop time from the spTree
    double stopTime = spTree->getCurrentTime();
    double eventTime = NAN;
    bool isSpeciation = 0;
    // start a new locus tree

    lociTree = std::shared_ptr<LocusTree>( new LocusTree(numTaxaToSim,
                                                        currentSimTime,
                                                        geneBirthRate,
                                                        geneDeathRate,
                                                        transferRate));


    // species tree is read in from r so convert index from R to C++ indexing
    std::map<int, int> rToTdckenIndxMap = spTree->makeIndxMap();

    spTree->switchIndicesFirstToSecond(rToTdckenIndxMap);
    
    // get the root
    auto spRoot = spTree->getRoot();

    std::map<int, std::string> tipMap = spTree->makeTipMap();
    // set the locus tree index to line up with the species tree
    lociTree->getRoot()->setLindx(spRoot->getIndex());
    // get a map of <index of node in species tree, death time of that node>
    std::map<int,double> speciesDeathTimes = spTree->getDeathTimesFromNodes();
    // make a set of species that are currently alive to keep track of which lineages
    // in speciesDeathTimes are around
    std::vector<int> contempSpecies;
    // placeholder for descendents of a member of contempSpecies
    std::pair<int, int> sibs;
    // set the stop time
    lociTree->setStopTime(stopTime);
    // if contempSpecies is not empty for some reason clear it
    if(!(contempSpecies.empty()))
        contempSpecies.clear();
    // insert the root of spTree index as the first species index to simulate within
    contempSpecies.push_back(spRoot->getIndex());

    while(currentSimTime < stopTime){
      // get time to next event based on rate parameters and number extant tips
      eventTime = lociTree->getTimeToNextEvent();
      currentSimTime += eventTime;
      // loop through currently living species to see if a species death time has been reached
      // if so speciation occurs and all locus tree lineages associated with that species index
      // speciate also
      // if extinction occurs those lineages go extinct
      for(std::vector<int>::iterator it = contempSpecies.begin(); it != contempSpecies.end();){
        if(currentSimTime >= speciesDeathTimes[(*it)]){
          isSpeciation = spTree->macroEvent((*it));
          if(isSpeciation){
          //  currentSimTime = speciesDeathTimes[(*it)];
            // tipward traverse to get the descendants of (*it) put these in sibs
            sibs = spTree->preorderTraversalStep(*it);
            // speciate the members of lociTree with indx of (*it) give those descendants
            // either sibs[0] or sib[1] as index
            lociTree->speciationEvent((*it), speciesDeathTimes[(*it)], sibs);
            //erase (*it) from contempSpecies and insert the descendant 2 where it was
            //and descendant 1 to one space after
            it = contempSpecies.erase(it);
            it = contempSpecies.insert(it, sibs.second);
            it =  contempSpecies.insert(it, sibs.first);



          }
          else{

            // species extinction or no event
            // check if the member of spTree is extinct or extant in final tree
            if(!(spTree->getIsExtantFromIndx(*it))){
             // currentSimTime = speciesDeathTimes[(*it)];
              lociTree->extinctionEvent(*it, speciesDeathTimes[(*it)]);
              it = contempSpecies.erase(it);
            }
            else{
              // otherwise move along
              it = contempSpecies.erase(it);
            }
          }
        }
        else{
          // otherwise move along
          ++it;
        }
      }
      // currentSimTime = placeholder;
      // if we pass stopTime change to be at stopTime and break
      if(currentSimTime >= stopTime){
        currentSimTime = stopTime;
        break;
      }
      // we get to 0 living nodes end sm
      if(lociTree->getNumTips() < 1){ // TODO: this should be refactored because the function name makes no sense
        treesComplete = false;
        return treesComplete;
      }
      // locus tree event (assuming parameters are NON-ZERO)
      // if parameters are all 0 no locus tree events occur and we get
      // the species tree back (relabeld)
      if(lociTree->checkLocusTreeParams()){
        lociTree->ermEvent(currentSimTime);
      }

    }

    currentSimTime = stopTime;
    // set the present time
    lociTree->setPresentTime(currentSimTime);

    // set the names based on their species ID, so tips are named
    // "T<SPECIES_INDX>_<LOCUS {A,B,C,...}>"
    lociTree->setNamesBySpeciesID(tipMap);

    treesComplete = true;

    return treesComplete;
}
// wrapper around locus tree sim to make sure we get a proper tree
bool Simulator::simLocusTree(){
  bool good = false;

  while(!good){
    good = bdsaBDSim();
  }
  return good;

}

// function to get the epochs of the locus tree (i.e. our coalescent breakpoints)
std::set<double, std::greater<double> > Simulator::getEpochs(){
    std::set<double, std::greater<double> > epochs;
    auto lociTreeNodes = lociTree->getNodes();
    for(auto it = lociTreeNodes.begin(); it != lociTreeNodes.end(); ++it){
        if(!((*it)->getIsExtinct())){
            if((*it)->getIsTip())
                epochs.insert((*it)->getDeathTime());
            epochs.insert((*it)->getBirthTime());
        }
        else
            epochs.insert((*it)->getDeathTime());
    }
    return epochs;
}

// multispecies coalescent simulator
bool Simulator::coalescentSim(){
    bool treeGood = false;
    geneTree = std::shared_ptr<GeneTree>(new GeneTree(numTaxaToSim, indPerPop, popSize, generationTime));
    std::map<int,int> spToLo;

    int ancIndx = -1;
    int epochCount = 0;

    double stopTime = NAN;
    double stopTimeEpoch = NAN;
    double stopTimeLoci = NAN;
    bool allCoalesced = false, deathCheck = false;
    bool is_ext = 0;
  // get the coalescent breakpoints from lociTree
    std::set<double, std::greater<double> > epochs = getEpochs();
    // how many epochs
    int numEpochs = (int) epochs.size();
  // get the indices of extinct loci
    std::set<int> extinctFolks = lociTree->getExtLociIndx();
    // this function is for adding in multilocus coalescent later
    std::set<int> coalescentBounds = lociTree->getCoalBounds();
  // get ContempLoci - the ones alive at the end of the locus tree sim (tips at present)
    std::vector< std::vector<int> > contempLoci = lociTree->getExtantLoci(epochs);
    // get the stop times of loci as a map with indices as keys
    std::map<int, double> stopTimes = lociTree->getBirthTimesFromNodes();
    // intialize the tree with individuals sampled from contempLoci and starting with the first coalescent bound
    geneTree->initializeTree(contempLoci, *(epochs.begin()));
    std::set<int>::iterator extFolksIt;
    //loop through the epochs in order
    for(std::set<double, std::greater<double> >::iterator epIter = epochs.begin(); epIter != epochs.end(); ++epIter){
        // set the time as the current epoch
        currentSimTime = *epIter;
      // if we aren't in the last epoch
        if(epochCount != numEpochs - 1){
          // iterate
            epIter = std::next(epIter, 1);
            stopTimeEpoch = *epIter;
            // loop through the contempLoci that are around during this epoch
            for(unsigned int j = 0; j < contempLoci[epochCount].size(); ++j){
                // find the extinct folks in this epoch
                extFolksIt = extinctFolks.find(contempLoci[epochCount][j]);
                // check if there were any found extinct
                is_ext = (extFolksIt != extinctFolks.end());
                // if so
                if(is_ext){
                  // add tips for the extinct species
                    geneTree->addExtinctSpecies(currentSimTime, contempLoci[epochCount][j]);
                  // erase from the extinct folks to mark that it was added
                    extinctFolks.erase(extFolksIt);
                }
                // get the stopTime
                stopTimeLoci = stopTimes[contempLoci[epochCount][j]];
                // if the current stop time is greater than the stop time of the epoch
                // it will not go extinct during this epoch so deathCheck keeps track of that
                if(stopTimeLoci > stopTimeEpoch){
                    stopTime = stopTimeLoci;
                    deathCheck = true;
                }
                else{
                    stopTime = stopTimeEpoch;
                    deathCheck = false;
                }
                // get the index of the ancestor of contempLoci[epochCount][j]
                ancIndx = lociTree->postOrderTraversalStep(contempLoci[epochCount][j]);
                // run the censored coalescent on memebers of geneTree with Lindx of contempLoci[epochCount][j]
                allCoalesced = geneTree->censorCoalescentProcess(currentSimTime,
                                                                 stopTime,
                                                                 contempLoci[epochCount][j],
                                                                 ancIndx,
                                                                 deathCheck);


                // if all coalesced remove that loci from the matrix of loci
                if(allCoalesced){
                    int check = contempLoci[epochCount][j];
                    for(unsigned int k = epochCount + 1; k <  (unsigned) numEpochs; k++){
                        for(unsigned m = 0; m < contempLoci[k].size(); ++m){
                            if(contempLoci[k][m] == check){
                                contempLoci[k].erase(contempLoci[k].begin() + m);
                                break;
                            }
                        }
                    }
                }
                // reset these
                allCoalesced = false;
                is_ext = false;
            }
            // go back one in epIter
            epIter = std::prev(epIter, 1);
        }
        else{
          // if we are in the last epoch do a coalescent until we have one lineage
            // finish coalescing
            geneTree->rootCoalescentProcess(currentSimTime);
            treeGood = true;
            geneTree->setBranchLengths();
        }
        // post incremenet
        epochCount++;

    }

    return treeGood;
}

// wrapper around simGeneTree takes the index of geneTrees as an argument and
// places the result of Simulator::coalescentSim into geneTrees[j]
// this assumes that most are simulating >1 geneTrees
bool Simulator::simGeneTree(int j){
  bool gGood = false;
  RNGScope scope;

  while(!gGood){
    gGood = coalescentSim();
  }
  geneTrees[j] = geneTree;
  return gGood;
}


// ##################################
//
// below here is I think dead code?
double Simulator::calcSpeciesTreeDepth(){
    return spTree->getTreeDepth();
}

double Simulator::calcExtantSpeciesTreeDepth(){
    std::shared_ptr<SpeciesTree> tt = std::shared_ptr<SpeciesTree> (new SpeciesTree(numTaxaToSim));
    spTree->getRootFromFlags(false);
    tt->setRoot(spTree->getExtantRoot());
    tt->setExtantRoot(tt->getRoot());
    tt->reconstructTreeFromSim(tt->getExtantRoot());
    double extTreeDepth = tt->getTreeDepth();
    tt = nullptr;
    return extTreeDepth;
}

double Simulator::calcLocusTreeDepth(int i){
    return locusTrees[i]->getTreeDepth();
}

int Simulator::findNumberTransfers(){
    int numTrans = 0;
    for(unsigned int i = 0; i < locusTrees.size(); i++){
        numTrans += locusTrees[i]->getNumberTransfers();
    }
    return numTrans;
}

// wrappers for SpeciesTree, symbionTree, lociTree, geneTree root edge calculation
double Simulator::getSpeciesTreeRootEdge(){
  return spTree->getRoot()->getDeathTime() - spTree->getRoot()->getBirthTime();
}

double Simulator::getLocusTreeRootEdge(){
  return lociTree->getRoot()->getDeathTime() - lociTree->getRoot()->getBirthTime();
}

double Simulator::getSymbiontTreeRootEdge(){
  return symbiontTree->getRoot()->getDeathTime() - symbiontTree->getRoot()->getBirthTime();
}


double Simulator::getGeneTreeRootEdge(int j){
  return geneTrees[j]->getRoot()->getBranchLength();
}
