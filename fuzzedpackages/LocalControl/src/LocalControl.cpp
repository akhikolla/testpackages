/**
 * This file contains the c++ code for the R package localControl.
 * The contents of the file can be broken into 3 or 4 sections which are marked below.
 * Eventually, each of the sections should be broken up and linked together with an appropriate header file.
 *
 * @author Nick Lauve, Christophe Lambert, Robert Obenchain, Stan Young
 */

// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <unistd.h>

using namespace Rcpp;

// START: Helper functions
// Takes a vector of elements, returns ordered set of unique elements.
template<class T>
std::vector<double> getUniqueElements(T const &elements){

  std::set<double> uniques(elements.begin(), elements.end());
  return std::vector<double>(uniques.begin(), uniques.end());

}

// Takes a 2 dimensional vector and returns a dataframe for R output.
// Function is specialized to write "rad_" headers currently.
//  For the sake of modularity, we should change this to have the sting header passed as a parameter.
template <class T>
DataFrame v2d2df(const std::vector< std::vector <T> > &v2d){
  size_t nc = v2d.size();
  List dfl(nc);
  CharacterVector colnames(nc);
  std::stringstream nameBuilder;

  for(size_t i = 0; i < nc; i++){
    dfl[i] = v2d[i];

    nameBuilder.str("");
    nameBuilder << "rad_" << i+1;
    colnames[i] = nameBuilder.str();
  }

  DataFrame resFrame(dfl);
  resFrame.attr("names") = colnames;

  return resFrame;
}

// Used to store the distance between two patients.
struct DistanceElement{
  double dist;
  int idx;

  DistanceElement(){}
  DistanceElement(double dist1, int idx1): dist(dist1), idx(idx1){}
};

// Used to calculate the L2 norm between cluster variable vectors for two patients.
double distanceL2(std::vector<double> const &v1, std::vector<double> const &v2){
  std::size_t i;
  double sum = 0.0;

  for(i=0; i < v1.size(); i++){
    sum += (v1[i]-v2[i])*(v1[i]-v2[i]);
  }
  return(sqrt(sum));
}

// Used to sort distance elements.
bool compareDistanceElement(const DistanceElement &lhs, const DistanceElement &rhs){
  return (lhs.dist < rhs.dist);
}
// END: Helper functions

// START: Param & result objects
// These are data structures for localControl parameters and returned values.
// For this section, thinking about breaking this code into smaller files,
// I think we should put all survival oriented objects/functions into one file, then do the same for cross-sectional.
class Patient{
public:
  Patient(const std::vector<double> &cvars, double out, bool treat):
  clusterVars(cvars), outcome(out), treatment(treat){}

  std::vector<double> clusterVars;
  double outcome;
  bool treatment;

};

class SurvivalPatient{
public:
  SurvivalPatient( const std::vector<double> &cvars, double ot, int ti, int ri, bool cen, bool trt):
  clusterVars(cvars), outcomeTime(ot), timeIndex(ti), riskIndex(ri), censored(cen), treatment(trt){}

  std::vector<double> clusterVars;
  double outcomeTime;
  int timeIndex;
  int riskIndex;
  bool censored;
  bool treatment;
};

struct LCResult{
  std::vector<std::vector<double> > t1Outcomes;
  std::vector<std::vector<double> > t0Outcomes;
  std::vector<std::vector<int> > t1Counts;
  std::vector<std::vector<int> > t0Counts;
  //std::vector< std::vector< std::vector<bool> > > assignments
};

struct SurvivalResult{
  std::vector< std::vector< std::vector<double> > > t0HappenVects;
  std::vector< std::vector< std::vector<double> > > t1HappenVects;
  std::vector< std::vector<double> > t0CensorVects;
  std::vector< std::vector<double> > t1CensorVects;
  std::vector<int> informative;
  std::vector<int> t0Informative;
  std::vector<int> t1Informative;
};
// END: Param & result objects

// START: localControl execution objects
// This abstract class, 'LocalController' contains all of the parameters which we send from R, as well as the
//  mutex and counters which were originally contained in global scope. We create the derived clases:
//  'CSController', and 'SurvivalController' to specify the behavior of the functions in the innermost
//  section of the algorithm. (basically, how and where we're saving the results of the near-neighborhoods)
class LocalController{

public:

  size_t getPatientCount(){
    return numPatients;
  }

  unsigned int getThreadCount(){
    return numThreads;
  }

  void setLimits(std::vector<double> const &lims){
    limits = lims;
    numLimits = lims.size();
  }

  int getLoopCount(){
    pthread_mutex_lock(&lMutex);
    int lc = loopCounter;
    loopCounter = loopCounter + 1;
    pthread_mutex_unlock(&lMutex);
    return lc;
  }

  void initLC(){
    loopCounter = 0;
    threadCounter = 0;
    pthread_mutex_init(&lMutex, NULL);
    pthread_mutex_init(&tMutex, NULL);
    lockThreads();
  }

  void addThread(){
    threadCounter++;
  }

  void removeThread(){
    pthread_mutex_lock(&tMutex);
    threadCounter--;
    pthread_mutex_unlock(&tMutex);
  }

  void runLC(){
    unlockThreads();

    while(threadCounter > 0){
      usleep(1000000);
    }
  }

  virtual void setPatients(const DataFrame& patients) = 0;

  virtual void initResults(unsigned int nt) = 0;

  virtual void addCluster(int patIdx, int threadIdx) = 0;

  virtual List getResults() = 0;

protected:
  static const int treatmentColIdx = 0;
  static const int outcomeColIdx = 1;

  size_t numPatients, numLimits, clusterVarCount;
  std::vector<double> limits;
  unsigned int numThreads;

  std::vector<DistanceElement> getDistances(int centerIdx){
    std::vector<DistanceElement> distances(numPatients);

    for(size_t otherIdx = 0; otherIdx < numPatients; otherIdx++){
      distances[otherIdx].dist = distanceL2(getPatientVars(centerIdx), getPatientVars(otherIdx));
      distances[otherIdx].idx = otherIdx;
    }

    std::sort(distances.begin(),distances.end(), compareDistanceElement);
    return distances;
  }

  const virtual std::vector<double> &getPatientVars(int i) const = 0;

private:
  pthread_mutex_t lMutex;
  pthread_mutex_t tMutex;
  volatile int loopCounter;
  volatile int threadCounter;

  void lockThreads(){
    pthread_mutex_lock(&lMutex);
  }

  void unlockThreads(){
    pthread_mutex_unlock(&lMutex);
  }

};

class CSController: public LocalController{

public:

  virtual void setPatients(const DataFrame& patients){

    unsigned int i, j;

    numPatients = patients.nrows();
    patientVect.reserve(numPatients);

    NumericVector treatments = patients[treatmentColIdx];
    NumericVector outcomes = patients[outcomeColIdx];

    //Patients has 2 extra columns for treatment and outcome
    //After we add the vectors to clustercols, we subtract 2 from this value
    clusterVarCount = patients.size();

    std::vector<NumericVector> clustercols;
    for(j = 2; j < clusterVarCount; j++){
      clustercols.push_back(patients[j]);
    }

    clusterVarCount -=2;

    std::vector<double> cvars;
    cvars.resize(clusterVarCount);

    for(i = 0; i < numPatients; i++){
      for(j = 0; j < clusterVarCount; j++){
        cvars[j] = clustercols[j][i];
      }

      patientVect.push_back(Patient(cvars,outcomes[i],treatments[i]));
    }
  }

  virtual void initResults(unsigned int nt){

    numThreads = nt;

    results.t1Counts.resize(numLimits);
    results.t0Counts.resize(numLimits);
    results.t1Outcomes.resize(numLimits);
    results.t0Outcomes.resize(numLimits);
    for(size_t j = 0; j < numLimits; j ++){
      results.t1Counts[j].resize(numPatients);
      results.t0Counts[j].resize(numPatients);
      results.t1Outcomes[j].resize(numPatients);
      results.t0Outcomes[j].resize(numPatients);
    }
  }

  virtual void addCluster(int patIdx, int threadIdx){

    std::vector<DistanceElement> dists = getDistances(patIdx);

    size_t radIdx, nnIdx;

    for(radIdx = 0; radIdx < numLimits; radIdx++){

      // Add 'fully informative' notion here.
      // when n1 & n2 > 0, we can use cluster for variance calculations.

      //each iteration looks at a different rad
      double currentRad = limits[radIdx];

      double t1outcome = 0.0, t0outcome = 0.0;
      int n1Count = 0, n0Count = 0;
      nnIdx = 0;

      //check each patient in range, assign cluster values to our counters.
      while(nnIdx < dists.size() && dists[nnIdx].dist <= currentRad){
        if(patientVect[dists[nnIdx].idx].treatment){
          n1Count++;
          t1outcome += patientVect[dists[nnIdx].idx].outcome;
        }else {
          n0Count++;
          t0outcome += patientVect[dists[nnIdx].idx].outcome;
        }

        nnIdx++;
      }

      results.t1Counts[radIdx][patIdx] = n1Count;
      results.t0Counts[radIdx][patIdx] = n0Count;
      if(n1Count > 0 && n0Count > 0){
        results.t1Outcomes[radIdx][patIdx] += t1outcome / n1Count;
        results.t0Outcomes[radIdx][patIdx] += t0outcome / n0Count;
      } else{
        results.t0Outcomes[radIdx][patIdx] = std::numeric_limits<double>::quiet_NaN();
        results.t1Outcomes[radIdx][patIdx] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  virtual List getResults(){

    List countDF = List::create(Named("T1") = v2d2df(results.t1Counts),
                                            Named("T0") = v2d2df(results.t0Counts));

    List outcoDF = List::create(Named("T1") = v2d2df(results.t1Outcomes),
                                            Named("T0") = v2d2df(results.t0Outcomes));

    return List::create(Named("outcomes") = outcoDF,
                              Named("counts") = countDF);

  }

protected:
  const virtual std::vector<double> &getPatientVars(int patIdx) const {
    return patientVect[patIdx].clusterVars;
  }

private:
  std::vector<Patient> patientVect;
  LCResult results;

};

class SurvivalController: public LocalController{

public:
  void setCensorCode(int cc){
    censorCode =  cc;
  }

  void setFailTimes(const std::vector<double> &fts){
    failTimes = fts;
    uniqueFailTimes = getUniqueElements(fts);
    numTimes = uniqueFailTimes.size();
  }

  virtual void setPatients(const DataFrame& patients){

    t1_count = 0, t0_count = 0;

    std::size_t i,j;

    NumericVector treatments = patients[treatmentColIdx];
    NumericVector outcomes = patients[outcomeColIdx];

    clusterVarCount = patients.size();
    numPatients = patients.nrows();

    setRisks(outcomes);

    patientVect.reserve(numPatients);

    std::vector<double> cvars;
    cvars.resize(clusterVarCount);

    std::vector<NumericVector> clustercols;
    for(j = 2; j < clusterVarCount; j++){
      clustercols.push_back(patients[j]);
    }

    clusterVarCount -=2;

    for(i = 0; i < numPatients; i++){

      if(treatments[i] == 1){
        t1_count++;
      }else{
        t0_count++;
      }

      for(j = 0; j < clusterVarCount; j++){
        cvars[j] = clustercols[j][i];
      }

      int tmpRisk = outcomes[i];
      bool tmpCensor = (tmpRisk == censorCode);

      if(!tmpCensor){
        tmpRisk = getRiskIdx(outcomes[i]);
      }

      double tmpTime = failTimes[i];
      int uniqueTimeIndex = getUTI(tmpTime);

      patientVect.push_back(SurvivalPatient(cvars,
                                            tmpTime,
                                            uniqueTimeIndex,
                                            tmpRisk,
                                            tmpCensor,
                                            treatments[i]));

    }

  }

  virtual void initResults(unsigned int nt){
    numThreads = nt;
    results.resize(numThreads);

    for(size_t iThread = 0; iThread < numThreads; iThread++){
      results[iThread].t0HappenVects.resize(numLimits, std::vector<std::vector<double> >(numRisks, std::vector<double>(numTimes)));
      results[iThread].t1HappenVects.resize(numLimits, std::vector<std::vector<double> >(numRisks, std::vector<double>(numTimes)));
      results[iThread].t1CensorVects.resize(numLimits, std::vector<double> (numTimes));
      results[iThread].t0CensorVects.resize(numLimits, std::vector<double> (numTimes));
      results[iThread].informative.resize(numLimits,0);
      results[iThread].t0Informative.resize(numLimits,0);
      results[iThread].t1Informative.resize(numLimits,0);
    }
  }

  virtual void addCluster(int patIdx, int threadIdx){

    std::vector<DistanceElement> const &distances = getDistances(patIdx);

    size_t radIdx, nnIdx;

    for(radIdx = 0; radIdx < numLimits; radIdx++){
      double currentRad = limits[radIdx];

      std::vector<std::vector<double> > t1Happens(numRisks, std::vector<double>(numTimes));
      std::vector<std::vector<double> > t0Happens(numRisks, std::vector<double>(numTimes));
      std::vector<double> t1Censors(numTimes);
      std::vector<double> t0Censors(numTimes);

      int n1 = 0;
      int n0 = 0;
      nnIdx = 0;

      while(nnIdx < distances.size() && distances[nnIdx].dist <= currentRad){
        int uniqueIdx = patientVect[distances[nnIdx].idx].timeIndex;
        if(patientVect[distances[nnIdx].idx].treatment){
          n1++;
          if(patientVect[distances[nnIdx].idx].censored) {
            t1Censors[uniqueIdx]++;
          } else {
            t1Happens[patientVect[distances[nnIdx].idx].riskIndex][uniqueIdx]++;
          }
        } else {
          n0++;
          if(patientVect[distances[nnIdx].idx].censored) {
            t0Censors[uniqueIdx]++;
          } else {
            t0Happens[patientVect[distances[nnIdx].idx].riskIndex][uniqueIdx]++;
          }
        }
        nnIdx++;
      }

      if(n1 > 0 && n0 > 0){
        results[threadIdx].informative[radIdx]++;
        if(patientVect[patIdx].treatment){
          results[threadIdx].t1Informative[radIdx]++;
        } else {
          results[threadIdx].t0Informative[radIdx]++;
        }
        for(size_t j = 0;j < numTimes; j++){
          results[threadIdx].t1CensorVects[radIdx][j] += t1Censors[j]/(double) n1;
          results[threadIdx].t0CensorVects[radIdx][j] += t0Censors[j]/(double) n0;
          for(size_t l = 0; l < numRisks; l++){
            results[threadIdx].t0HappenVects[radIdx][l][j] += t0Happens[l][j]/(double)n0;
            results[threadIdx].t1HappenVects[radIdx][l][j] += t1Happens[l][j]/(double)n1;
          }
        }
      }
    }
  }

  virtual List getResults(){

    size_t i,j,k,l;

    //merge vectors from all threads
    for(i = 1; i < numThreads; i++){
      for(j = 0; j < numLimits; j++){
        results[0].informative[j] += results[i].informative[j];
        results[0].t0Informative[j] += results[i].t0Informative[j];
        results[0].t1Informative[j] += results[i].t1Informative[j];
        for(k = 0; k < numTimes; k++){
          for(l = 0; l < numRisks; l++){
            results[0].t0HappenVects[j][l][k] += results[i].t0HappenVects[j][l][k];
            results[0].t1HappenVects[j][l][k] += results[i].t1HappenVects[j][l][k];
          }
          results[0].t0CensorVects[j][k] += results[i].t0CensorVects[j][k];
          results[0].t1CensorVects[j][k] += results[i].t1CensorVects[j][k];
        }
      }
    }

    SurvivalResult oneResult = results[0];

    List events(2);
    CharacterVector tNames(2);
    tNames[0] = "T0";
    tNames[1] = "T1";

    List t1(numLimits);
    List t0(numLimits);
    CharacterVector limitNames(numLimits);
    CharacterVector eventNames(numRisks + 1);
    std::stringstream nameBuilder;

    for(size_t lIdx = 0; lIdx < numLimits; lIdx++){

      nameBuilder.str("");
      nameBuilder << "rad_" << (lIdx + 1);
      limitNames[lIdx] = nameBuilder.str();

      List radEvents1(numRisks + 1);
      List radEvents0(numRisks + 1);

      for(size_t rIdx = 0; rIdx < numRisks; rIdx++){
        radEvents1[rIdx] = oneResult.t1HappenVects[lIdx][rIdx];
        radEvents0[rIdx] = oneResult.t0HappenVects[lIdx][rIdx];

        nameBuilder.str("");
        nameBuilder << "Failcode_" << riskTypes[rIdx];
        eventNames[rIdx] = nameBuilder.str();
      }
      eventNames[numRisks] = "Censored";

      radEvents1[numRisks] = oneResult.t1CensorVects[lIdx];
      radEvents0[numRisks] = oneResult.t0CensorVects[lIdx];

      radEvents1.attr("names") = eventNames;
      radEvents0.attr("names") = eventNames;

      t1[lIdx] = radEvents1;
      t0[lIdx] = radEvents0;

    }

    t1.attr("names") = limitNames;
    t0.attr("names") = limitNames;

    events[0] = t0;
    events[1] = t1;
    events.attr("names") = tNames;

    List outList(6);
    CharacterVector outListNames(6);

    outListNames[0] = "Failtimes";
    outListNames[1] = "Events";
    outListNames[2] = "CIF";
    outListNames[3] = "SDR";
    outListNames[4] = "KM";
    outListNames[5] = "NumInf";

    outList[0] = uniqueFailTimes;
    outList[1] = events;
    outList.attr("names") = outListNames;

    //Vectors to store cumulative incidence curves for all risks x rads
    std::vector<std::vector<std::vector<double> > > cIF0s(numRisks, std::vector<std::vector<double> >(numLimits, std::vector<double>(numTimes)));
    std::vector<std::vector<std::vector<double> > > cIF1s(numRisks, std::vector<std::vector<double> >(numLimits, std::vector<double>(numTimes)));

    //Vectors to store standard error estimates for each of the cumulative incidence functions
    std::vector<std::vector<std::vector<double> > > sdr0s(numRisks, std::vector<std::vector<double> >(numLimits, std::vector<double>(numTimes)));
    std::vector<std::vector<std::vector<double> > > sdr1s(numRisks, std::vector<std::vector<double> >(numLimits, std::vector<double>(numTimes)));

    //Vectors storing the kaplan meier curves at each rad
    std::vector<std::vector<double> > t0KMs(numLimits, std::vector<double>(numTimes));
    std::vector<std::vector<double> > t1KMs(numLimits, std::vector<double>(numTimes));


    for(i = 0; i < numLimits; i++){
      
      SurvivalResult sResult = oneResult;

      for(k = 0; k < numTimes; k++){
        for(l = 0; l < numRisks; l++){
          sResult.t0HappenVects[i][l][k] *= (double) sResult.t0Informative[i] / (double) sResult.informative[i];
          sResult.t1HappenVects[i][l][k] *= (double) sResult.t1Informative[i] / (double) sResult.informative[i];
        }
        sResult.t0CensorVects[i][k] *= (double) sResult.t0Informative[i] / (double) sResult.informative[i];
        sResult.t1CensorVects[i][k] *= (double) sResult.t1Informative[i] / (double) sResult.informative[i];
      }
      
      std::vector<std::vector<double> > sp0s(numRisks, std::vector<double>(numTimes));
      std::vector<std::vector<double> > sp1s(numRisks, std::vector<double>(numTimes));
      
      std::vector<std::vector<double> > ws1s(numRisks, std::vector<double>(numTimes));
      std::vector<std::vector<double> > ws0s(numRisks, std::vector<double>(numTimes));
      
      std::vector<double> survSum1(numTimes);
      std::vector<double> survSum0(numTimes);

      std::vector<double> sv0s(numRisks);
      std::vector<double> sv1s(numRisks);

      std::vector<double> sq1s(numRisks);
      std::vector<double> sq0s(numRisks);

      ///Time 0, everyone is alive.
      t0KMs[i][0] = 1.0;
      t1KMs[i][0] = 1.0;

      double atrisk1 = sResult.t1Informative[i];
      double atrisk0 = sResult.t0Informative[i];

      int max0 = 1, max1 = 1;

      for(j=1; j < numTimes; j++){

        double events1 = 0.0, events0 = 0.0;

        for(l = 0; l < numRisks; l++){
          events1 += sResult.t1HappenVects[i][l][j];
          events0 += sResult.t0HappenVects[i][l][j];
        }

        if(events0 > 0 || sResult.t0CensorVects[i][j] > 0){
          max0 = j;
        }
        if(events1 > 0 || sResult.t1CensorVects[i][j] > 0){
          max1 = j;
        }
        
        t0KMs[i][j] = t0KMs[i][j-1] * ((atrisk0)-events0)/(atrisk0);
        t1KMs[i][j] = t1KMs[i][j-1] * ((atrisk1)-events1)/(atrisk1);
        
        survSum0[j] = survSum0[j-1] + (events0 / ((atrisk0)*((atrisk0)-events0)));
        survSum1[j] = survSum1[j-1] + (events1 / ((atrisk1)*((atrisk1)-events1)));

        for(l = 0; l < numRisks; l++){

          sp0s[l][j] = (sResult.t0HappenVects[i][l][j] / atrisk0) * t0KMs[i][j-1];
          sp1s[l][j] = (sResult.t1HappenVects[i][l][j] / atrisk1) * t1KMs[i][j-1];
         
          cIF0s[l][i][j] = cIF0s[l][i][j-1] + sp0s[l][j];
          cIF1s[l][i][j] = cIF1s[l][i][j-1] + sp1s[l][j];
   
          double wp0 = ((atrisk0 - sResult.t0HappenVects[i][l][j]) / (sResult.t0HappenVects[i][l][j]*atrisk0)) + survSum0[j-1];
          double wp1 = ((atrisk1 - sResult.t1HappenVects[i][l][j]) / (sResult.t1HappenVects[i][l][j]*atrisk1)) + survSum1[j-1];
          
          if(std::isinf(wp0)) wp0 = 0;
          if(std::isinf(wp1)) wp1 = 0;
          
          ws0s[l][j] = (-1/atrisk0) + survSum0[j-1];
          ws1s[l][j] = (-1/atrisk1) + survSum1[j-1];
                              
          if(std::isinf(ws0s[l][j])) ws0s[l][j] = 0;
          if(std::isinf(ws1s[l][j])) ws1s[l][j] = 0;

          double nc1 = 0, nc0 = 0;
          
          for(size_t n = 0; n < j; n++){ 
            nc1 += sp1s[l][n] * sp1s[l][j] * ws1s[l][n]; 
            nc0 += sp0s[l][n] * sp0s[l][j] * ws0s[l][n]; 
          }

          sq1s[l] += nc1 * 2; 
          sq0s[l] += nc0 * 2;
          
          sv1s[l] += sp1s[l][j] * sp1s[l][j] * wp1;
          sv0s[l] += sp0s[l][j] * sp0s[l][j] * wp0;
          
          sdr1s[l][i][j] = sqrt(sq1s[l] + sv1s[l]);
          sdr0s[l][i][j] = sqrt(sq0s[l] + sv0s[l]);
        
        }
        atrisk0 -= (events0 + sResult.t0CensorVects[i][j] );
        atrisk1 -= (events1 + sResult.t1CensorVects[i][j] );
      }

      for(j=max0+1 ; j < numTimes; j++){
        t0KMs[i][j] = std::numeric_limits<double>::quiet_NaN();
        for(l = 0; l < numRisks; l++){
          cIF0s[l][i][j] = std::numeric_limits<double>::quiet_NaN();
          sdr0s[l][i][j] = std::numeric_limits<double>::quiet_NaN();
        }
      }

      for(j=max1+1 ; j < numTimes; j++){
        t1KMs[i][j] = std::numeric_limits<double>::quiet_NaN();
        for(l = 0; l < numRisks; l++){
          cIF1s[l][i][j] = std::numeric_limits<double>::quiet_NaN();
          sdr1s[l][i][j] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
    
    //Rcerr << "\n";
    List CIFList(numRisks);
    List ERRList(numRisks);
    CharacterVector riskNames(numRisks);

    for(size_t i = 0; i < numRisks; i++){

      nameBuilder.str("");
      nameBuilder << "Failcode_" << riskTypes[i];
      riskNames[i] = nameBuilder.str();

      CIFList[i] = List::create(Named("T1") = v2d2df(cIF1s[i]),
                                Named("T0") = v2d2df(cIF0s[i]));

      ERRList[i] = List::create(Named("T1") = v2d2df(sdr1s[i]),
                                Named("T0") = v2d2df(sdr0s[i]));
    }

    CIFList.attr("names") = riskNames;
    ERRList.attr("names") = riskNames;

    outList[2] = CIFList;
    outList[3] = ERRList;
    outList[4] = List::create(Named("T1") = v2d2df(t1KMs),
                              Named("T0") = v2d2df(t0KMs));

    outList[5] = oneResult.informative;
    return outList;
}

protected:
  const virtual std::vector<double> &getPatientVars(int patIdx) const {
    return patientVect[patIdx].clusterVars;
  }

private:
  std::vector<SurvivalPatient> patientVect;
  std::vector<SurvivalResult> results;
  std::vector<double> uniqueFailTimes;
  std::vector<double> failTimes;
  std::vector<double> riskTypes;
  int censorCode;
  size_t numRisks, numTimes, t1_count, t0_count;

  void setRisks(const NumericVector &risks){
    std::vector<double> uniqueRisks = getUniqueElements(risks);
    for(std::size_t i = 0; i < uniqueRisks.size(); i++){
      if(uniqueRisks[i] != censorCode){
        riskTypes.push_back(uniqueRisks[i]);
      }
    }
    numRisks = riskTypes.size();
  }

  //unique times are sorted, so we can use lower_bound
  int getUTI(double ft){
    std::vector<double>::iterator lb = lower_bound(uniqueFailTimes.begin(), uniqueFailTimes.end(), ft);
    return (lb - uniqueFailTimes.begin());
  }

  //find to get the correct risk index because they might not be sorted
  int getRiskIdx(double risk){
    std::vector<double>::iterator fin = find(riskTypes.begin(),riskTypes.end(),risk);
    return (fin - riskTypes.begin());
  }

};
// END: localControl execution objects

// START: localControl Main (R-Bridge)
// The following functions are used to get things running.
// I think 'lcMainloop', and 'runLocalControl' can/should eventually be moved inside the
// localController objects. (TODO)
void lcMainLoop(LocalController &lci, unsigned int tIdx){

  size_t loopIdx;
  size_t loopMax = lci.getPatientCount();

  while(1){

    loopIdx = lci.getLoopCount();

    if(loopIdx >= loopMax){
      break;
    }
    else{
      if(loopIdx % 1000 == 0){
        Rcerr<<"Processing " << loopIdx << "\n";
      }

      lci.addCluster(loopIdx, tIdx);
    }
  }

  lci.removeThread();

}

struct ThreadParam{
  int threadNum;
  LocalController* lc;
};

void *lcThreadOp(void *params){

  ThreadParam *ps = (ThreadParam *) params;
  LocalController &lci = *(ps->lc);
  int myThread = ps->threadNum;

  lcMainLoop(lci, myThread);

  pthread_exit(NULL);

  return NULL;
}

//This function spins off the threads to do the work for local control.
int runLocalControl(LocalController &lci){

  int nt = lci.getThreadCount();

  // std::vector<pthread_t> threads = std::vector<pthread_t>(nt);
  // std::vector<ThreadParam> tp = std::vector<ThreadParam>(nt);
  pthread_t *threads = new pthread_t[nt];
  ThreadParam *tp = new ThreadParam[nt];
  lci.initLC();

  for(int i=0; i < nt; i++){

    tp[i].lc = &lci;
    tp[i].threadNum = i;

    int rc = pthread_create(&threads[i], NULL, lcThreadOp, (void *) &tp[i]);

    if(rc){
      Rcerr << "Error: unable to create thread, " << rc << "\n";
      return -1;
    }

    lci.addThread();

  }

  lci.runLC();

  for(int i=0; i < nt; i++){
    pthread_join(threads[i], NULL);
  }

  delete [] tp;
  delete [] threads;
  return 0;

}

// The following two functions are visible in R-package-scope.
// These functions are not exported to the user's global scope, it is assumed that
// they will only ever be called from our package functins in localControl.R.
// Attempting to call these functions manually will probably result in R crashing.

// [[Rcpp::export]]
List newLC(DataFrame& patients, std::vector<double> limits, unsigned int numThreads)
{
  CSController params;
  params.setPatients(patients);
  params.setLimits(limits);
  params.initResults(numThreads);

  if(runLocalControl(params) == -1){
     return List();
  }
  else return params.getResults();
}

// [[Rcpp::export]]
List newCRLC(DataFrame& patients, std::vector<double> limits, std::vector<double> fTimes, int cenCode, unsigned int numThreads)
{
  SurvivalController params;
  params.setCensorCode(cenCode);
  params.setFailTimes(fTimes);
  params.setPatients(patients);
  params.setLimits(limits);
  params.initResults(numThreads);

  if(runLocalControl(params) == -1){
     return List();
  }
  else return params.getResults();
}
// END: localControl Main (R-Bridge)
