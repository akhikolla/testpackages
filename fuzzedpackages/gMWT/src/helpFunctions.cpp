#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath> 
#include <R.h>
#include <Rmath.h>
//#include <stl>

/* This is quite ugly code, but it still speeds up the calculation
 */

using namespace std;

extern "C" {

 
    double getP(double *x, double *y, int *nx, int *ny){
    int NX = nx[0];
    int NY = ny[0];
    int i;
    int j;
    double resAbs = 0;
    double resRel;
    
    for(i=0;i<NX;i++){
     for(j=0;j<NY;j++){
       if(x[i] < y[j]) resAbs = resAbs + 1; 
       if(x[i] == y[j]) resAbs = resAbs + 0.5;
     }
    }
    resRel = (double) resAbs/(NX*NY);
    return(resRel);
  }

  
  void getPR(double *x, double *y, int *nx, int *ny,double *result){
    int NX = nx[0];
    int NY = ny[0];
    int i;
    int j;
    double resAbs = 0;
    double resRel;
    
    for(i=0;i<NX;i++){
     for(j=0;j<NY;j++){
       if(x[i] < y[j]) resAbs = resAbs + 1; 
       if(x[i] == y[j]) resAbs = resAbs + 0.5;
     }
    }
    resRel = (double) resAbs/(NX*NY);
    result[0] = resRel;
  }

   void getPTripR(double *x, double *y, double *z, int *nx, int *ny, int *nz, double *result){
    int NX = nx[0];
    int NY = ny[0];
    int NZ = nz[0];
    int i;
    int j;
    int k;
    double resAbs = 0;
    double resRel;
    
    for(i=0;i<NX;i++){  
     for(j=0;j<NY;j++){
      for(k=0;k<NZ;k++){
       if( (x[i]<y[j]) && (y[j]<z[k])) {
	  resAbs = resAbs + 1; }
	  else if( (x[i]==y[j]) && (y[j]<z[k])){
	    resAbs = resAbs + 0.5;}
	    else if( (x[i]<y[j]) && (y[j]==z[k])){
	      resAbs = resAbs + 0.5;}
	      else if( (x[i]==y[j]) && (y[j]==z[k])){ 
		resAbs = resAbs + 1/6;}
      }
     }
    }
    resRel = (double) resAbs/(NX*NY*NZ);
    result[0] = resRel;
  }
  
  /* This function is a deasater as well, I have to optmize it later!!!!
   * It might be even so, that this is numerical not stable, when N is getting big!!!
   */
  //void varU(double *x, double *y, int *nx, int *ny, double *result ){
    double varU(double *x, double *y, int *nx, int *ny){
    int i;
    int N = nx[0] + ny[0];
    //double joined [N];
    std::vector<double> joinedVec(N);
     
    for(i=0;i<nx[0];i++) joinedVec[i] = x[i];
    for(i=0;i<ny[0];i++) joinedVec[i+nx[0]] = y[i];
  
    // Write the vector into an array - so the old code still works...
    double* joined = &joinedVec[0];
    
    //Sorting the joined array
    int elements = sizeof(joined) / sizeof(joined[0]);
    sort(joined,joined+elements);
   
    int diffValues = 0;
    // We assume, that all values are different and that we do not have any ties:
    //int mult [N];
    std::vector<double> multVec(N);
    // Write the vector into an array - so the old code still works...
    double* mult = &multVec[0];

    for(i=0;i<N;i++) mult[i] = 0.0;
    mult[0] = 1.0;
    
    for(i=1;i<N;i++){
     if(joined[i]!=joined[i-1]){
	diffValues += 1;
	mult[diffValues] = 1.0;    
     } else {
	mult[diffValues] += 1.0;
      }
    }
    // calculate the correction term
    double innerSum;
    innerSum = 0;
    for(i=0;i<N;i++){
      innerSum += pow(mult[i],3)-mult[i];
    }
    // Remember here: In the reference we use the Variance for W, and here we use proportions!!!
    innerSum = innerSum/(12*N*(N-1)*nx[0]*ny[0]);
    
    double totalTerm = (double) (N + 1) /( 12* nx[0] * ny[0]);
    
    totalTerm -= innerSum;
    //result[0] = totalTerm;
    return(totalTerm);
  }
  
   void varUR(double *x, double *y, int *nx, int *ny, double *result ){
   //double varU(double *x, double *y, int *nx, int *ny){
    int i;
    int N = nx[0] + ny[0];
    //double joined [N];
    std::vector<double> joinedVec(N);
    for(i=0;i<nx[0];i++) joinedVec[i] = x[i];
    for(i=0;i<ny[0];i++) joinedVec[i+nx[0]] = y[i];
  
     // Write the vector into an array - so the old code still works...
    double* joined = &joinedVec[0];
    
    //Sorting the joined array
    int elements = sizeof(joined) / sizeof(joined[0]);
    sort(joined,joined+elements);
   
    int diffValues = 0;
    // We assume, that all values are different and that we do not have any ties:
    //int mult [N];
    std::vector<double> multVec(N);
    double* mult = &multVec[0];
    
    for(i=0;i<N;i++) mult[i] = 0;
    mult[0] = 1;
    
    for(i=1;i<N;i++){
     if(joined[i]!=joined[i-1]){
	diffValues += 1;
	mult[diffValues] = 1;    
     } else {
	mult[diffValues] += 1;
      }
    }
    // calculate the correction term
    double innerSum;
    innerSum = 0;
    for(i=0;i<N;i++){
      innerSum += pow(mult[i],3)-mult[i];
    }
    // Remember here: In the reference we use the Variance for W, and here we use proportions!!!
    innerSum = innerSum/(12*N*(N-1)*nx[0]*ny[0]);
    
    double totalTerm = (double) (N + 1) /( 12* nx[0] * ny[0]);
    
    totalTerm -= innerSum;
    result[0] = totalTerm;
    //return(totalTerm);
  }
  
  double getT(double *x, double *y, int *nx, int *ny){
    double output;
    output = (getP(x,y,nx,ny) - 0.5)/(sqrt(varU(x,y,nx,ny)));
    return(output);
  }
  
  double getRho(int nx, int ny, int nz){
    double nxD = (double) nx;
    double nyD = (double) ny;
    double nzD = (double) nz;
    //double nom = (double) sqrt(nx * ny);
    double nom = sqrt(nxD * nyD);
    double den = sqrt((nxD + nzD + 1) * (nyD + nzD + 1));
    double output = nom / den;
    return(output);
  }

  void getRhoR(int *nx, int *ny, int *nz, double *result){
    double nxD = (double) nx[0];
    double nyD = (double) ny[0];
    double nzD = (double) nz[0];
    
    //double nom = (double) sqrt(nx[0] * ny[0]);
    //double den = (double) sqrt((nx[0] + nz[0] + 1) * (ny[0] + nz[0] + 1));
    double nom = sqrt(nxD * nyD);
    double den = sqrt((nxD + nzD + 1) * (nyD + nzD + 1));

    double output = nom / den;
    result[0] = output;
    result[1] = nom;
    result[2] = den;
    result[3] = nx[0];
    result[4] = ny[0];
    result[5] = nz[0];
  }
  
  int determineTS(double ts1, double ts2, double rho){
    int output;
    
    if((ts1<=(rho*ts2)) && (ts2<=(rho*ts1))){
       output = 3;
     } else if( (ts1<=0) && (ts2 > (rho*ts1))){
       output = 2;
     } else if((ts2<=0) && (ts1 > (rho*ts2))){
       output = 4;
     } else output = 1;    
   
    return(output);
  }
 
  /* This function calculates the test statistic for the Union-Intersection test
   */
  void uitR(double *x, double *y, double *z, int *nx, int *ny, int *nz, double *result ){
    
    double ts1 = getT(x,z,nx,nz);
    double ts2 = getT(y,z,ny,nz);
    double rho = getRho(nx[0],ny[0],nz[0]);
   // double rho = 1/(12*nz[0]);  
    int quadrant = (int) determineTS(ts1,ts2,rho);

    double output;
    
    if(quadrant==3){
	output = 0;
    } else if(quadrant == 2){
	output = pow((ts2-rho*ts1),2)/(1-pow(rho,2));
    } else if(quadrant == 4){
	output = pow((ts1-rho*ts2),2)/(1-pow(rho,2));
    } else {
	/* WARNING!!!! I wrote first this line, but I think it's wrong and it has to be a sum, NOT a difference!!!
	 */
	  output = ts1*(ts1-ts2*rho) + ts2*(ts2 - ts1 * rho);
    }
  
  result[0] = output;
  result[1] = quadrant;
  result[2] = ts1;
  result[3] = ts2;
  result[4] = rho;
 
  }

  double uit(double *x, double *y, double *z, int *nx, int *ny, int *nz){
    
    double ts1 = getT(x,z,nx,nz);
    double ts2 = getT(y,z,ny,nz);
    double rho = getRho(nx[0],ny[0],nz[0]);  
    //double rho = 1/(12 * nz[0]);
    int quadrant = (int) determineTS(ts1,ts2,rho);

    double output;
    
    if(quadrant==3){
	output = 0;
    } else if(quadrant == 2){
	output = pow((ts2-rho*ts1),2)/(1-pow(rho,2));
    } else if(quadrant == 4){
	output = pow((ts1-rho*ts2),2)/(1-pow(rho,2));
    } else {
	/* WARNING!!!! I wrote first this line, but I think it's wrong and it has to be a sum, NOT a difference!!!
	 */
	 output = ts1*(ts1-ts2*rho) + ts2*(ts2 - ts1 * rho);
	 
	//output = ts1*(ts1+ts2*rho) + ts2*(ts2 + ts1 * rho);
    }
  
  return(output); 
  }

  
  // ADD HERE STILL THE OPTION FOR A SEED!!!!
  // CHECK IF THE MEMORY IS FREED!!!
  double *permObs(double *x, double *y, double *z,int *nx, int *ny, int *nz){
  
      int i;
      int N = nx[0] + ny[0] + nz[0];
  
      // join the three arrays
      double *joined = new double [N];  
      
      for(i=0;i<nx[0];i++) joined[i] = x[i];
      for(i=0;i<ny[0];i++) joined[i+nx[0]] = y[i];
      for(i=0;i<nz[0];i++) joined[i+nx[0]+ny[0]] = z[i];
  
      GetRNGstate();
      
      // Shuffle the joined array:
      for (int i=0; i<(N-1); i++) {
         double rv = unif_rand();
	 int r = i + floor(rv * (N-i));
         //int r = i + (rand() % (N-i)); // Random remaining position.
         double temp = joined[i]; 
	 joined[i] = joined[r];
	 joined[r] = temp;
        }

      PutRNGstate();
      
      return(joined);
  }

   // ADD HERE STILL THE OPTION FOR A SEED!!!!
  void permObsR(double *x, double *y, double *z,int *nx, int *ny, int *nz,double *result){

      int i;
      int N = nx[0] + ny[0] + nz[0];
  
      // join the three arrays
      //double joined [N];
      std::vector<double> joined(N);
      
      for(i=0;i<nx[0];i++) joined[i] = x[i];
      for(i=0;i<ny[0];i++) joined[i+nx[0]] = y[i];
      for(i=0;i<nz[0];i++) joined[i+nx[0]+ny[0]] = z[i];
  
      double temp;
      
      GetRNGstate();

      // Shuffle the joined array:
      for (int i=0; i<(N-1); i++) {
	double rv = unif_rand();
         int r = i + floor(rv * (N-i));
         //int r = i + (rand() % (N-i)); // Random remaining position.
         temp = joined[i]; 
	 joined[i] = joined[r];
	 joined[r] = temp;
        }
        
      PutRNGstate();      
     for(i=0;i<N;i++) result[i] = joined[i];
  }
  
  void permObs2R(double *x, double *y,int *nx, int *ny,double *result){

      int i;
      int N = nx[0] + ny[0];
  
      // join the three arrays
      //double joined [N];
      std::vector<double> joined(N);
      
      for(i=0;i<nx[0];i++) joined[i] = x[i];
      for(i=0;i<ny[0];i++) joined[i+nx[0]] = y[i];
  
      double temp;
      
      GetRNGstate();
      // Shuffle the joined array:
      for (int i=0; i<(N-1); i++) {
	 double rv = unif_rand();
	 int r = i + floor(rv * (N-i));
         // int r = i + (rand() % (N-i)); // Random remaining position.
         temp = joined[i]; 
	 joined[i] = joined[r];
	 joined[r] = temp;
        }
     PutRNGstate();     
     for(i=0;i<N;i++) result[i] = joined[i];
  }
  
  void permUIT(double *x, double *y, double *z,int *nx, int *ny, int *nz, int *nper, double *result){
    int N=nx[0]+ny[0]+nz[0]; 
    //double output[nper[0]];
    std::vector<double> output(N);
    //double tempValues [N];
    std::vector<double> tempValues(N);
    
    //double tempX [nx[0]];
    std::vector<double> tempXVec(nx[0]);
    double* tempX = &tempXVec[0];

    //double tempY [ny[0]];
    std::vector<double> tempYVec(ny[0]);
    double* tempY = &tempYVec[0];

    //double tempZ [nz[0]];
    std::vector<double> tempZVec(nz[0]);
    double* tempZ = &tempZVec[0];

    
    int i;
    int j;
    int k;
    
    
    
    for(i=0;i<nper[0];i++){
      for(j=0;j<N;j++){
      tempValues[j]= permObs(x,y,z,nx,ny,nz)[j];

      for(k=0;k<nx[0];k++) tempX[k] = tempValues[k];
      for(k=0;k<ny[0];k++) tempY[k] = tempValues[k + nx[0]];
      for(k=0;k<nz[0];k++) tempZ[k] = tempValues[k + nx[0] + ny[0]];
      
      output[i] = uit(tempX,tempY,tempZ,nx,ny,nz);

      }
    }
    for(i=0;i<nper[0];i++) result[i] = output[i];    
  }    
  
  int sgn(double d){
    return d<0? -1 : d>0;
    }

  
  
}
