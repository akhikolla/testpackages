#include <Rcpp.h>

#include "DynamicPrograms.h"

#include "Data.h"
#include "DataGauss.h"
#include "DataHsmuce.h"
#include "DataMDependentPS.h"
#include "DataJsmurf.h"
#include "DataJsmurfPS.h"
#include "DataJsmurfLR.h"
#include "DataHjsmurf.h"
#include "DataHjsmurfSPS.h"
#include "DataHjsmurfLR.h"
#include "DataLR.h"
#include "Data2Param.h"

#include "IntervalSystem.h"
#include "IntervalSystemAll.h"
#include "IntervalSystemAllLengths.h"
#include "IntervalSystemDyaLen.h"
#include "IntervalSystemDyaLenLengths.h"
#include "IntervalSystemDyaPar.h"
#include "IntervalSystemDyaParLengths.h"

using namespace Rcpp;

// [[Rcpp::export(name = ".callRoutines")]]
RObject callRoutines(RObject observations,
                     int routineType, List argumentsListRoutine,
                     int dataType, List argumentsListData,
                     int intervalSystemType, List argumentsListIntervalSystem) {
  Data* data = NULL;
  
  switch (dataType) {
  case 0:
    DataGauss::setData(observations, argumentsListData);
    data = new DataGauss();
    break;
  case 10:
    DataMDependentPS::setData(observations, argumentsListData);
    data = new DataMDependentPS();
    break;
  case 11:
    DataJsmurf::setData(observations, argumentsListData);
    data = new DataJsmurf();
    break;
  case 12:
    DataJsmurfPS::setData(observations, argumentsListData);
    data = new DataJsmurfPS();
    break;
  case 13:
    DataJsmurfLR::setData(observations, argumentsListData);
    data = new DataJsmurfLR();
    break;
  case 20:
    DataHsmuce::setData(observations);
    data = new DataHsmuce();
    break;
  case 21:
    DataHjsmurf::setData(observations, argumentsListData);
    data = new DataHjsmurf();
    break;
  case 22:
    DataHjsmurfSPS::setData(observations, argumentsListData);
    data = new DataHjsmurfSPS();
    break;
  case 23:
    DataHjsmurfLR::setData(observations, argumentsListData);
    data = new DataHjsmurfLR();
    break;
  case 100:
    Data2Param::setData(observations, argumentsListData);
    data = new Data2Param();
    break;
  case 102:
    DataLR::setData(observations, argumentsListData);
    data = new DataLR();
    break;
  default:
    stop("dataType %d is not defined", dataType);
  }
  
  if (dataType >= 100) {
    RObject ret;
    
    switch (routineType) {
    case 1:
      ret = computeStatistic(data, argumentsListRoutine);
      break;
    case 10:
      ret = findSmallScales(data, argumentsListRoutine);
      break;
    default:
      data -> cleanUpStaticVariables();
      delete data;
      stop("only computeStat can be called for this parametric family");
    }
    data -> cleanUpStaticVariables();
    delete data;
    
    return ret;
  }

  IntervalSystem* intervalSystem = NULL;
  switch(intervalSystemType) {
  case 0:
    intervalSystem = new IntervalSystemAll(data -> getN());
    break;
  case 1:
    intervalSystem = new IntervalSystemAllLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  case 10:
    intervalSystem = new IntervalSystemDyaLen(data -> getN());
    break;
  case 11:
    intervalSystem = new IntervalSystemDyaLenLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  case 20:
    intervalSystem = new IntervalSystemDyaPar(data -> getN());
    break;
  case 21:
    intervalSystem = new IntervalSystemDyaParLengths(data -> getN(), argumentsListIntervalSystem);
    break;
  default:
    data -> cleanUpStaticVariables();
    delete data;
    stop("intervalSystemType %d is not defined", intervalSystemType);
  }
  
  RObject ret;
  switch (routineType) {
  case 0:
    ret = intervalSystem -> computeMultiscaleStatisticNull(data);
    break;
  case 1:
    ret = intervalSystem -> computeMultiscaleStatistic(data, argumentsListRoutine);
    break;        
  case 2:
    Data::setCriticalValues(argumentsListRoutine);
    ret = intervalSystem -> computeBounds(data); 
    break;
  case 3:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitSimpleDynamicProgram(data, intervalSystem);
    break;
  case 4:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitIntervalDynamicProgram(data, intervalSystem);  
    break;
  case 5:
    Data::setCriticalValues(argumentsListRoutine);
    ret = fitBandDynamicProgram(data, intervalSystem);
    break;
  default:
    delete intervalSystem;
    data -> cleanUpStaticVariables();
    delete data;
    stop("routineType %d is not defined", routineType);
  }
  
  delete intervalSystem;
  data -> cleanUpStaticVariables();
  delete data;
  
  return ret;
}
