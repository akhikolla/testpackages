#include <Rcpp.h> 
#include <math.h>
#include <time.h> 
#include <map>
#include <vector> 
#include <limits>

#include "constants.h"          //Global Constants
#include "functions.h"          //Functions
#include "structs.h"
#include "mtrand.h"
#include "convexHull_Monotone.h"

//Namespaces
using namespace Rcpp;
using namespace std;

//Forward Delcarations
NumericVector           checkAndUniqueLevels(NumericVector& input);
IntegerMatrix           checkAndSortDM(IntegerMatrix& dm, NumericMatrix& xyz);
NumericMatrix           pertubate(const IntegerMatrix& dm, NumericMatrix& xyz, const NumericVector& levels, double upToPercent);


Vec3 origin(0,0,0);
bool compare_Vec3_sort(const Vec3& a, const Vec3& b){
  return atan2(a.y-origin.y,a.x-origin.x) < atan2(b.y-origin.y,b.x-origin.x);
}


// [[Rcpp::export]]
NumericMatrix contourWalker(IntegerMatrix& dm, NumericMatrix& xyz, NumericVector& levels, 
                            double criticalRatio = 10.0, 
                            double maximumPertubation=1e-5){
  //Create Local Copies & Variables
  IntegerMatrix dmNew     = IntegerMatrix(dm);  
  NumericMatrix xyzNew    = NumericMatrix(xyz); 
  NumericVector levelsNew = NumericVector(levels); 
  
  //Run Checks
  dmNew         = checkAndSortDM(dmNew,xyzNew);
  levelsNew     = checkAndUniqueLevels(levelsNew);
  xyzNew        = pertubate(dmNew,xyzNew,levelsNew,maximumPertubation);
  
  //Variables
  vector<Node3>         nodes;
  vector<Node3>         chull;
  vector<int>           chullIndexes;
  vector<Vec3>          contour;
  vector<ContourData>   contourData;
  vector<Del>           dels;
  
  //Build Nodes
  for(int i = 0; i < xyzNew.nrow(); i++)
    nodes.push_back(Node3(i,Vec3(xyzNew(i,0), xyzNew(i,1), xyzNew(i,2))));
    
  //Determine the FULL Convex Hull, ie inclusive of points lying on the hull
  bool  includeColinear   = true,   //Include All Points on the Hull
        isZeroBased       = true;   //C++
  chullIndexes = convexHullAM_IndexesVector(xyzNew.column(0),xyzNew.column(1),includeColinear,isZeroBased);
  for(unsigned int i = 0; i < chullIndexes.size(); i++)
    chull.push_back( nodes[ chullIndexes[i] ] );
  
  //Build Deleaunay Triangles
  for(int i = 0; i < dmNew.nrow(); i++)
    dels.push_back(Del(nodes[dmNew(i,0)],nodes[dmNew(i,1)],nodes[dmNew(i,2)]));
  
  //Establish the Network
  //Network Sisters (Dels that share an edge), 
  //If they are not a sister, try and network them as a cousin (Dels that share a vertex)
  for(size_t i = 0, ds = dels.size(); i < ds - 1; i++){
    Del* di = &dels[i];
    for(size_t j = i+1; !(di->isFull()) && j < ds; j++){
      Del* dj = &dels[j];
      di->makeSisters(dj);
    }
  }
  
  //Draw the Contours and pass the result back to the contour container by reference
  size_t levelID, groupID, delID, pathID;
  for(levelID = 0; levelID < (size_t)levelsNew.size(); levelID++){
      for(delID = 0, groupID = 0; delID < dels.size(); delID++){
        double level = levelsNew[levelID];
        (&(dels[delID]))->drawContour(level,chull,contour,NULL,dels,criticalRatio);
        if(contour.size() > 2){
          for(pathID = 0; pathID < contour.size(); pathID++)
            contourData.push_back(ContourData(levelID,groupID,pathID,contour[pathID].x,contour[pathID].y,contour[pathID].z));
          groupID++;
        }
        contour.clear();
      }
      for(delID = 0; delID < dels.size(); delID++){
        dels[delID].reset();
      }
  }
  
  //Pre-allocate the result array
  NumericMatrix result(contourData.size(),6);
  
  //Put the contourData vector into numeric matrix and return the result 
  for(pathID = 0; pathID < contourData.size(); pathID++)
    result.row(pathID) = (&contourData[pathID])->toNumericVector();
  
  //Done
  return result;
}


NumericVector checkAndUniqueLevels(NumericVector& input){
  if(input.size() == 0)
    throw std::out_of_range("The Levels Vector is Empty, please specify at least one level to contour.");
  vector<double> tmp = vector<double>();
  for(int i = 0; i < input.size(); i++)
    tmp.push_back(input(i));
  uniqueOnly<double>(tmp);
  return NumericVector(tmp.begin(),tmp.end());
}

IntegerMatrix checkAndSortDM(IntegerMatrix& dm, NumericMatrix& xyz){
  if(dm.ncol() != 3 || xyz.ncol() != 3)
    throw std::out_of_range("Expecting 3 columns in both 'dm' and 'xyz' matrixes.");
  
  //Check 1-based indexing of the dm matrix (coming from R)
  if(min(dm) == 1)
    dm += -1;
  
  //Check values in dm, will not raise out of bounds exception in xyz
  if(max(dm) >= xyz.nrow())
    throw std::out_of_range("Values in 'dm' would result in out of bounds errors in xyz.");
  
  IntegerMatrix dmSorted(dm.nrow(),dm.ncol());
  Centroid* centroids = new Centroid[dm.nrow()]();
  Vec3* vecs = new Vec3[3]();
  
  //Assemble the centroids array, from the vectors
  for(int i = 0; i < dm.nrow(); i++){
     getVecsByRefZeroBased(dm,xyz,i,vecs);
     centroids[i] = Centroid(i); //To Sort
      for(int j = 0; j < 3; j++)
        centroids[i] += vecs[j]/3.0;
  }
  
  //Execute the sort using the comparitor which is part of the centroids struct
  std::sort(centroids,centroids + dm.nrow());
  
  //Restore the data
  for(int i = 0; i < dm.nrow(); i++)
    dmSorted.row(i) = dm.row(centroids[i].index);
  dm = dmSorted;
  
  //Cleanup the dynamic arrays
  delete [] centroids;
  delete [] vecs;
  
  //Done, Return
  return dm;
}


inline double nudge(double input, double upToPcnt){
  MTRand_closed drand;
  double  r1 = -drand(), 
          r2 = -drand();
  return input*(1.0 + r1*abs(upToPcnt)/100.0) + r2*D_TOL;
}

NumericMatrix pertubate(const IntegerMatrix& dm, NumericMatrix& xyz, const NumericVector& levels, double upToPercent){
  int pb, i,j, limit = dm.nrow(), cnt=0;
  double pcnt = abs(upToPercent);
  while(cnt <= limit && abs(upToPercent - D_TOL) > 0){ 
    for(i = 0, pb = 0; i < dm.nrow(); i++){
      int    ix[3] = {dm(i,0),dm(i,1),dm(i,2)};
      double  z[3] = {xyz(ix[0],2),xyz(ix[1],2),xyz(ix[2],2)};
      
      //Check Points
      if(isEqual(z[0],z[1])){
        xyz(ix[1],2) = z[1] = nudge(z[1],pcnt); pb++;
      }
      
      //Check Points
      if(isEqual(z[0],z[2]) || isEqual(z[1],z[2])){
        xyz(ix[2],2) = z[2] = nudge(z[2],pcnt); pb++;
      }
       
      //Check levels
      for(j = 0; j < 3; j++){
        z[j] = xyz(ix[j],2);
        if(std::find(levels.begin(), levels.end(), z[j]) != levels.end()){
          xyz(ix[j],2) = nudge(z[j],pcnt); pb++; 
        }
      }

    }; 
    cnt++;
    if(pb == 0) break;
  }
  return xyz;
}




