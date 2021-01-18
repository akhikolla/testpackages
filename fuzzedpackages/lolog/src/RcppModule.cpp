#ifndef INSIDE

#include <Rcpp.h>
#include <BinaryNet.h>
#include <tests.h>
#include <LatentOrderLikelihood.h>

/*
 * Handles all functions and methods exported to R.
 */

RCPP_MODULE(lolog) {
    using namespace Rcpp;
    using namespace lolog;

    class_<DirectedNet >("DirectedNet")
    .constructor<Rcpp::IntegerMatrix,int>()
    .constructor<SEXP>()
    .method("clone",&DirectedNet::cloneR)
    .method("size",&DirectedNet::size)
    .method("isDirected",&DirectedNet::isDirected)
    .method("setDyads",&DirectedNet::setDyadsR)
    .method("getDyads",&DirectedNet::getDyadsR)
    .method("emptyGraph",&DirectedNet::emptyGraph)
    .method("edges",&DirectedNet::edgelistR1)
    .method("edges",&DirectedNet::edgelistR2)
    .method("[",&DirectedNet::getDyadMatrixR)
    .method("[<-",&DirectedNet::setDyadMatrixR)
    .method("variableNames",&DirectedNet::getVariableNamesR1)
    .method("variableNames",&DirectedNet::getVariableNamesR2)
    .method("[[",&DirectedNet::getVariableR)
    .method("getVariable",&DirectedNet::getVariableR)
    .method("getVariable",&DirectedNet::getVariableR1)
    .method("[[<-",&DirectedNet::setVariableR)
    .method("nMissing",&DirectedNet::nMissingR)
    .method("nEdges",&DirectedNet::nEdgesR1)
    .method("nEdges",&DirectedNet::nEdgesR2)
    .method("inDegree",&DirectedNet::indegreeR)
    .method("outDegree",&DirectedNet::outdegreeR)
    .method("outNeighbors",&DirectedNet::outneighborsR)
    .method("inNeighbors",&DirectedNet::inneighborsR)
    .method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR1)
    .method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR2)
    .method("setAllDyadsMissing",&DirectedNet::setAllDyadsMissingR3)
    ;

    class_<UndirectedNet >("UndirectedNet")
    .constructor<Rcpp::IntegerMatrix,int>()
    .constructor<SEXP>()
    .method("clone",&UndirectedNet::cloneR)
    .method("size",&UndirectedNet::size)
    .method("isDirected",&UndirectedNet::isDirected)
    .method("setDyads",&UndirectedNet::setDyadsR)
    .method("getDyads",&UndirectedNet::getDyadsR)
    .method("emptyGraph",&UndirectedNet::emptyGraph)
    .method("edges",&UndirectedNet::edgelistR1)
    .method("edges",&UndirectedNet::edgelistR2)
    .method("[",&UndirectedNet::getDyadMatrixR)
    .method("[<-",&UndirectedNet::setDyadMatrixR)
    .method("variableNames",&UndirectedNet::getVariableNamesR1)
    .method("variableNames",&UndirectedNet::getVariableNamesR2)
    .method("[[",&UndirectedNet::getVariableR)
    .method("getVariable",&UndirectedNet::getVariableR)
    .method("getVariable",&UndirectedNet::getVariableR1)
    .method("[[<-",&UndirectedNet::setVariableR)
    .method("nMissing",&UndirectedNet::nMissingR)
    .method("nEdges",&UndirectedNet::nEdgesR1)
    .method("nEdges",&UndirectedNet::nEdgesR2)
    .method("degree",&UndirectedNet::degreeR)
    .method("neighbors",&UndirectedNet::neighborsR)
    .method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR1)
    .method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR2)
    .method("setAllDyadsMissing",&UndirectedNet::setAllDyadsMissingR3)
    ;

    class_<Model<Undirected> >("UndirectedModel")
    .constructor()
    .constructor< Model<Undirected> >()
    .method("setNetwork",&Model<Undirected>::setNetworkR)
    .method("getNetwork",&Model<Undirected>::getNetworkR)
    .method("addStatistic",&Model<Undirected>::addStatistic)
    .method("addOffset",&Model<Undirected>::addOffset)
    .method("calculate",&Model<Undirected>::calculate)
    .method("statistics",&Model<Undirected>::statisticsR)
    .method("names",&Model<Undirected>::names)
    .method("offset",&Model<Undirected>::offset)
    .method("thetas",&Model<Undirected>::thetasR)
    .method("setThetas",&Model<Undirected>::setThetas)
    .method("setVertexOrder",&Model<Undirected>::setVertexOrderVector)
    .method("getVertexOrder",&Model<Undirected>::getVertexOrderVector)
    .method("isIndependent",&Model<Undirected>::isIndependent)
    ;
    class_<Model<Directed> >("DirectedModel")
    .constructor()
    .constructor< Model<Directed> >()
    .method("setNetwork",&Model<Directed>::setNetworkR)
    .method("getNetwork",&Model<Directed>::getNetworkR)
    .method("addStatistic",&Model<Directed>::addStatistic)
    .method("addOffset",&Model<Directed>::addOffset)
    .method("calculate",&Model<Directed>::calculate)
    .method("statistics",&Model<Directed>::statisticsR)
    .method("names",&Model<Directed>::names)
    .method("offset",&Model<Directed>::offset)
    .method("thetas",&Model<Directed>::thetasR)
    .method("setThetas",&Model<Directed>::setThetas)
    .method("setVertexOrder",&Model<Directed>::setVertexOrderVector)
    .method("getVertexOrder",&Model<Directed>::getVertexOrderVector)
    .method("isIndependent",&Model<Directed>::isIndependent)
    ;

    class_<LatentOrderLikelihood<Undirected> >("UndirectedLatentOrderLikelihood")
    .constructor< Model<Undirected> >()
    .method("setModel",&LatentOrderLikelihood<Undirected>::setModel)
    .method("getModel",&LatentOrderLikelihood<Undirected>::getModelR)
    .method("setThetas",&LatentOrderLikelihood<Undirected>::setThetas)
    .method("variationalModelFrame",&LatentOrderLikelihood<Undirected>::variationalModelFrame)
    .method("variationalModelFrameWithFunc",&LatentOrderLikelihood<Undirected>::variationalModelFrameWithFunc)
    .method("generateNetwork",&LatentOrderLikelihood<Undirected>::generateNetwork)
    ;

    class_<LatentOrderLikelihood<Directed> >("DirectedLatentOrderLikelihood")
    .constructor< Model<Directed> >()
    .method("setModel",&LatentOrderLikelihood<Directed>::setModel)
    .method("getModel",&LatentOrderLikelihood<Directed>::getModelR)
    .method("setThetas",&LatentOrderLikelihood<Directed>::setThetas)
    .method("variationalModelFrame",&LatentOrderLikelihood<Directed>::variationalModelFrame)
    .method("variationalModelFrameWithFunc",&LatentOrderLikelihood<Directed>::variationalModelFrameWithFunc)
    .method("generateNetwork",&LatentOrderLikelihood<Directed>::generateNetwork)
    ;

    function("initLologStatistics",&initStats);

    function("registerDirectedStatistic",&registerDirectedStatistic);
    function("registerUndirectedStatistic",&registerUndirectedStatistic);
    function("registerDirectedStatistic",&registerDirectedOffset);
    function("registerUndirectedStatistic",&registerUndirectedOffset);

    function("runLologCppTests",&tests::runLologTests);
}

#endif

