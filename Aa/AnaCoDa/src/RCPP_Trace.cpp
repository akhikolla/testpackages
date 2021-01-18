
#ifndef STANDALONE
#include "include/base/Trace.h"
#include <Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(Trace_mod)
{
  class_<Trace>( "Trace" )
  
    //Getter Functions:
    .method("getStdDevSynthesisRateAcceptanceRateTrace", &Trace::getStdDevSynthesisRateAcceptanceRateTrace)
    .method("getSynthesisRateTrace", &Trace::getSynthesisRateTrace)
    .method("getSynthesisRateAcceptanceRateTrace", &Trace::getSynthesisRateAcceptanceRateTrace)
    .method("getCodonSpecificAcceptanceRateTraceForAA", &Trace::getCodonSpecificAcceptanceRateTraceForAA)
    .method("getCodonSpecificAcceptanceRateTraceForCodon", &Trace::getCodonSpecificAcceptanceRateTraceForCodon)
    .method("getMixtureAssignmentTrace", &Trace::getMixtureAssignmentTrace)
    .method("getCodonSpecificAcceptanceRateTrace", &Trace::getCodonSpecificAcceptanceRateTrace)
    .method("getNseRateSpecificAcceptanceRateTrace", &Trace::getNseRateSpecificAcceptanceRateTrace)
    .method("getMixtureProbabilitiesTrace", &Trace::getMixtureProbabilitiesTrace)
    .method("getExpectedSynthesisRateTrace", &Trace::getExpectedSynthesisRateTrace)
    .method("getSynthesisOffsetAcceptanceRateTrace", &Trace::getSynthesisOffsetAcceptanceRateTrace)
    .method("getSynthesisOffsetAcceptanceRateTraceForIndex", &Trace::getSynthesisOffsetAcceptanceRateTraceForIndex)
    .method("getCodonSpecificParameterTrace", &Trace::getCodonSpecificParameterTraceByParamType)
    .method("getSynthesisRateAcceptanceRateTraceByMixtureElementForGene",
            &Trace::getSynthesisRateAcceptanceRateTraceByMixtureElementForGeneR)
    .method("getSynthesisRateTraceForGene", &Trace::getSynthesisRateTraceForGeneR)
    .method("getSynthesisRateTraceByMixtureElementForGene", &Trace::getSynthesisRateTraceByMixtureElementForGeneR)
    .method("getMixtureAssignmentTraceForGene", &Trace::getMixtureAssignmentTraceForGeneR)
    .method("getMixtureProbabilitiesTraceForMixture", &Trace::getMixtureProbabilitiesTraceForMixtureR)
    .method("getStdDevSynthesisRateTraces", &Trace::getStdDevSynthesisRateTraces)
    .method("getNumberOfMixtures", &Trace::getNumberOfMixtures)
    

    //Setter Functions:
    .method("setStdDevSynthesisRateTraces", &Trace::setStdDevSynthesisRateTraces)
    .method("setStdDevSynthesisRateAcceptanceRateTrace", &Trace::setStdDevSynthesisRateAcceptanceRateTrace)
    .method("setSynthesisRateTrace", &Trace::setSynthesisRateTrace)
    .method("setSynthesisRateAcceptanceRateTrace", &Trace::setSynthesisRateAcceptanceRateTrace)
    .method("setMixtureAssignmentTrace", &Trace::setMixtureAssignmentTrace)
    .method("setMixtureProbabilitiesTrace", &Trace::setMixtureProbabilitiesTrace)
    .method("setCodonSpecificAcceptanceRateTrace", &Trace::setCodonSpecificAcceptanceRateTrace)
    .method("setNseRateSpecificAcceptanceRateTrace", &Trace::setNseRateSpecificAcceptanceRateTrace)


    //ROC Specific:
    .method("getCodonSpecificParameterTraceByMixtureElementForCodon",
            &Trace::getCodonSpecificParameterTraceByMixtureElementForCodonR)
    .method("getSynthesisOffsetTrace", &Trace::getSynthesisOffsetTraceR)
    .method("getObservedSynthesisNoiseTrace", &Trace::getObservedSynthesisNoiseTraceR)
    .method("setSynthesisOffsetTrace", &Trace::setSynthesisOffsetTrace)
    .method("setSynthesisOffsetAcceptanceRateTrace", &Trace::setSynthesisOffsetAcceptanceRateTrace)
    .method("setObservedSynthesisNoiseTrace", &Trace::setObservedSynthesisNoiseTrace)
    .method("setCodonSpecificParameterTrace", &Trace::setCodonSpecificParameterTrace)

    //PANSE Specific
    .method("resizeNumberCodonSpecificParameterTrace", &Trace::resizeNumberCodonSpecificParameterTrace)
    .method("getPartitionFunctionTraces",&Trace::getPartitionFunctionTraces)
    .method("setPartitionFunctionTraces",&Trace::setPartitionFunctionTraces)


    //FONSE Specific
    .method("getInitiationCostTrace",&Trace::getInitiationCostTrace)
    .method("setInitiationCostTrace", &Trace::setInitiationCostTrace)
    ;
}
#endif

