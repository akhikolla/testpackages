// MultiTrainWithInventory.cpp file
 /*
 * Author: David Silksorth
 *         (c) 2014-2018 OpenReliability.org
 */


#include <R.h>
#include <Rcpp.h>
#include "stosim.h"


 SEXP MultiTrainWithInventoryCPP(SEXP arg1, SEXP arg2, SEXP arg3, SEXP arg4){

using namespace Rcpp;
//	src <- '
// set dataframe arguments from R into Rcpp vector classes for use in C++											
	Rcpp::NumericVector inTime(arg1);										
	Rcpp::NumericVector inDuration(arg2);										
	Rcpp::NumericVector GenRate(arg3);										
	Rcpp::NumericVector Constraints(arg4);										
	double CapacityHrs = Constraints[0];										
	double ReserveHrs = Constraints[1];										
	double RefillTime = Constraints[2];										
	double DischargeCap = Constraints[3];										
	double TurndownLimit = Constraints[4];										
	double TurndownTime = Constraints[5];										
	double numTrains = Constraints[6];										
// initialize the output object vectors											
	int insize = inTime.size();										
	int outsize = insize;										
	Rcpp::NumericVector outTime(outsize);										
	Rcpp::NumericVector outDuration(outsize);										
	Rcpp::NumericVector outProd(outsize);										
	Rcpp::IntegerVector DisCapEx(outsize);										
	Rcpp::IntegerVector RunOut(outsize);										
	Rcpp::IntegerVector EmptyOnD(outsize);										
	Rcpp::NumericVector EndInv(insize);										
	Rcpp::IntegerVector LogicPath(insize);										
// make first entries into out vectors											
	outTime[0]=0.0;										
	outDuration[0]=0.0;										
	outProd[0]=GenRate[0];										
	DisCapEx[0]=0;										
	RunOut[0]=0;										
	EmptyOnD[0]=0;										
	EndInv[0]=CapacityHrs;										
	LogicPath[0]=100;										
// separate pointer to the out vectors											
	int y=1;										
// declare variables											
	double LastProdRate;										
	double StartInventory;										
	double DemandRate;										
	double EndInventory;										
	double DurationToReserve;										
	double TurndownDemandRate;										
	int path_adder=0;										
											
// begin a line by line loop starting at line 2 of OpLineDetail											
	for(int x=1;  x < insize; x++)  {										
											
// establish a few terms for ease of use											
		LastProdRate =  outProd[y-1];									
		StartInventory = EndInv[x-1];									
											
											
// At Failure											
		if(GenRate[x] < GenRate[x-1])  {									
// exclude evaluation of non-operational system											
			if(LastProdRate > 0)  {								
				DemandRate = LastProdRate-GenRate[x];							
// exclude evaluation of excess demand											
				if( (DischargeCap- DemandRate*numTrains + 1e-9)>0 )  {							
// exclude evaluation of no inventory											
					if(StartInventory > 0)  {						
// establish whether above or below Reserve with more than a turndown margin											
// Inventory_consumption_before_turndown <- DemandRate*numTrains*TurndownTime											
						if(StartInventory > ReserveHrs + DemandRate*numTrains*TurndownTime)  {					
// Determine EndInventory											
							EndInventory = StartInventory-DemandRate*numTrains*inDuration[x];				
							if(EndInventory > ReserveHrs)  {				
// successful backup from inventory, no change to ProdRate											
								EndInv[x] = EndInventory;			
								LogicPath[x] = 104;			
							}else{				
// determine time at which Reserve level was hit											
								DurationToReserve =  (StartInventory - ReserveHrs)/(DemandRate*numTrains);			
// This case will end at Reserve inventory or below											
								if((GenRate[x] -TurndownLimit+1e-9)>0)  {			
// fill elements of the output vectors											
									outTime[y] = inTime[x] + DurationToReserve;		
//  the previous outDuration is known											
									outDuration[y-1] = outTime[y] - outTime[y-1];		
									outProd[y] = GenRate[x];		
									DisCapEx[y] = 0;		
									RunOut[y] = 0;		
									EmptyOnD[y] = 0;		
											
											
											
									y=y+1;		
// this is an end of the train state											
									EndInv[x] = ReserveHrs;		
									LogicPath[x] = 105;		
								}else{			
// fill elements of the output vectors											
									outTime[y] = inTime[x] + DurationToReserve;		
//  the previous outDuration is known											
									outDuration[y-1] = outTime[y] - outTime[y-1];		
									outProd[y] = TurndownLimit;		
									DisCapEx[y] = 0;		
									RunOut[y] = 0;		
									EmptyOnD[y] = 0;		
											
											
											
									y=y+1;		
											
// the train state case is not necessarily over yet											
									TurndownDemandRate = TurndownLimit - GenRate[x];		
									EndInventory = ReserveHrs - TurndownDemandRate*numTrains*(inDuration[x] - DurationToReserve);		
									if(EndInventory > 0)  {		
// this is an end of the train state											
											EndInv[x] = EndInventory;
											LogicPath[x] = 106;
									}else{		
// this case will runout based on TurndownDemandRate											
// fill elements of the output vectors											
// use the negative EndInventory to figure out when the runout occurred											
										outTime[y] = inTime[x] + inDuration[x] + EndInventory/(TurndownDemandRate*numTrains);	
//  the previous outDuration is known											
										outDuration[y-1] = outTime[y] - outTime[y-1];	
										outProd[y] = 0;	
										DisCapEx[y] = 0;	
										RunOut[y] = 1;	
										EmptyOnD[y] = 0;	
											
											
											
										y=y+1;	
// this is an end of the train state											
										EndInv[x] = 0;	
										LogicPath[x] = 107;	
									}		
								}			
							}				
						}else{					
// this case starts within Reserve, it is permitted to run for TurndownTime											
// Will duration of event permit full turndown?
											
							if(( inDuration[x] -TurndownTime -  1e-9) > 0)  {				
// Determine EndInventory											
								EndInventory = StartInventory-DemandRate*numTrains*TurndownTime;			
								if(EndInventory > 0)  {			
									if( GenRate[x] >= TurndownLimit) {		
											
										outTime[y] = inTime[x] + TurndownTime;	
//  the previous outDuration is known											
										outDuration[y-1] = outTime[y] - outTime[y-1];	
										outProd[y] = GenRate[x];	
										DisCapEx[y] = 0;	
										RunOut[y] = 0;	
										EmptyOnD[y] = 0;	
											
											
											
										y=y+1;	
											
// this is an end of the train state											
										EndInv[x] = EndInventory;	
										LogicPath[x] = 109;	
									}else{		
// there will still be demand since GenRate is below TurndownLimit											
// after setting a new set of conditions for wid entry the EndInventory needs to be redetermined											
											
											
										outTime[y] = inTime[x] + TurndownTime;	
											
										outDuration[y-1] = outTime[y] - outTime[y-1];	
										outProd[y] = TurndownLimit;	
										DisCapEx[y] = 0;	
										RunOut[y] = 0;	
										EmptyOnD[y] = 0;	
											
											
										y=y+1;	
											
// now re-determine the EndInventory based on remaining Duration											
// Determine EndInventory											
										TurndownDemandRate = TurndownLimit-GenRate[x];	
										EndInventory = EndInventory-TurndownDemandRate*numTrains*(inDuration[x]-TurndownTime);	
										if(EndInventory > 0)  {	
// successful backup from inventory, no further change to ProdRate											
											EndInv[x] = EndInventory;
											LogicPath[x] = 110;
										}else{	
// This case has been a run out at TurndownLimit											
// use the negative EndInventory to figure out when the runout occurred											
											
											outTime[y] = inTime[x] + inDuration[x] + EndInventory/(TurndownDemandRate*numTrains);
											
											outDuration[y-1] = outTime[y] - outTime[y-1];
											outProd[y] = 0;
											DisCapEx[y] = 0;
											RunOut[y] = 1;
											EmptyOnD[y] = 0;
											
											
											
											y=y+1;
// this is an end of the train state											
											EndInv[x] = 0;
											LogicPath[x] = 111;
										}	
									}		
								}else{			
// inventory runs out during turndown time											
											
											
									outTime[y] = inTime[x] + inDuration[x] + EndInventory/(DemandRate*numTrains);		
											
									outDuration[y-1] = outTime[y] - outTime[y-1];		
									outProd[y] = 0;		
									DisCapEx[y] = 0;		
									RunOut[y] = 1;		
									EmptyOnD[y] = 0;		
											
											
											
									y=y+1;		
											
									EndInv[x] = 0;		
									LogicPath[x] = 108;		
								}			
							}else{				
// Debugging revealed improper handling of incomplete turndown cases											
// Determine EndInventory											
								EndInventory = StartInventory-DemandRate*numTrains*inDuration[x];			
								if(EndInventory > 0)  {			
// in this case there is no change in ProdRate, no entry to wid											
									EndInv[x] = EndInventory;		
									LogicPath[x] = 112;		
								}else{			
// inventory runs out during turndown time											
// use the negative EndInventory to figure out when the runout occurred											
											
									outTime[y] = inTime[x] + inDuration[x] +  EndInventory/(DemandRate*numTrains);		
									outDuration[y-1] = outTime[y] - outTime[y-1];		
									outProd[y] = 0;		
									DisCapEx[y] = 0;		
									RunOut[y] = 1;		
									EmptyOnD[y] = 0;		
											
											
											
									y=y+1;		
											
									EndInv[x] = 0;		
											
									LogicPath[x] = 113;		
								}			
							}				
						}					
					}else{						
// A special case that may only occur if ReserveHrs=0, this is the EmptyOnDemand case											
											
						outTime[y] = inTime[x];					
											
						outDuration[y-1] = outTime[y] - outTime[y-1];					
						outProd[y] = 0;					
						DisCapEx[y] = 0;					
						RunOut[y] = 0;					
						EmptyOnD[y] = 1;					
											
											
											
						y=y+1;					
											
						EndInv[x] = 0;					
						LogicPath[x] = 103;					
					}						
				}else{							
// Demand has exceeded Discharge Capacity, sudden failure of downstream process											
											
					outTime[y] = inTime[x];						
											
					outDuration[y-1] = outTime[y] - outTime[y-1];						
					outProd[y] = 0;						
					DisCapEx[y] = 1;						
					RunOut[y] = 0;						
					EmptyOnD[y] = 0;						
											
											
											
					y=y+1;						
											
					EndInv[x] = StartInventory;						
					LogicPath[x] = 102;						
				}							
			}else{								
// ProdRate is already at zero, merely end the train state with same inventory											
				EndInv[x] = StartInventory;							
				LogicPath[x] = 101;							
			}								
											
											
		}else{									
// At Repair											
											
											
			path_adder=0;								
//  Most often this is a return to full GenLevel, if EndInventory will be based on refill											
			if(GenRate[x] == 1.0)   {								
// there should be no wid entry if ProdRate is unchanged (as successful backup)											
				if(LastProdRate<1.0)  {							
											
					outTime[y] = inTime[x];						
											
					outDuration[y-1] = outTime[y] - outTime[y-1];						
					outProd[y] = 1.0;						
					DisCapEx[y] = 0;						
					RunOut[y] = 0;						
					EmptyOnD[y] = 0;						
											
											
											
					y=y+1;						
					path_adder = 2;						
				}							
// this is an end of the train state											
// EndLevel is lessor of StartLevel + duration*CapacityHrs/RefillHrs, or CapacityHrs											
				if(StartInventory + inDuration[x]*CapacityHrs/RefillTime  <  CapacityHrs)   {							
				EndInv[x] = StartInventory + inDuration[x]*CapacityHrs/RefillTime;							
				LogicPath[x] = 201+path_adder;							
				}else{							
				EndInv[x] = CapacityHrs;							
				LogicPath[x] = 202+path_adder;							
				}							
			}else{								
				if(LastProdRate>GenRate[x])  {							
// there is still a demand from storage, treat similarly to fail events											
// an unusual case that can occur with large inventory											
					DemandRate = LastProdRate-GenRate[x];						
					if(StartInventory > ReserveHrs + DemandRate*numTrains*TurndownTime)  {						
// Determine EndInventory											
						EndInventory = StartInventory-DemandRate*numTrains*inDuration[x];					
						if(EndInventory > ReserveHrs)  {					
// successful backup from inventory, no change to ProdRate											
							EndInv[x] = EndInventory;				
							LogicPath[x] = 205;				
						}else{					
// determine time at which Reserve level was hit											
							DurationToReserve =  (StartInventory - ReserveHrs)/(DemandRate*numTrains);				
// This case will end at Reserve inventory or below											
							if(GenRate[x] >TurndownLimit)  {				
											
								outTime[y] = inTime[x] + DurationToReserve;			
											
								outDuration[y-1] = outTime[y] - outTime[y-1];			
								outProd[y] = GenRate[x];			
								DisCapEx[y] = 0;			
								RunOut[y] = 0;			
								EmptyOnD[y] = 0;			
											
											
											
								y=y+1;			
// this is an end of the train state											
								EndInv[x] = ReserveHrs;			
								LogicPath[x] = 206;			
							}else{				
											
								outTime[y] = inTime[x] + DurationToReserve;			
											
								outDuration[y-1] = outTime[y] - outTime[y-1];			
								outProd[y] = TurndownLimit;			
								DisCapEx[y] = 0;			
								RunOut[y] = 0;			
								EmptyOnD[y] = 0;			
											
											
											
								y=y+1;			
											
// the train state case is not necessarily over yet											
								TurndownDemandRate = TurndownLimit-GenRate[x];			
								EndInventory  =  ReserveHrs - TurndownDemandRate*numTrains*(inDuration[x] - DurationToReserve);			
								if(EndInventory > 0)  {			
// this is an end of the train state											
									EndInv[x] = EndInventory;		
									LogicPath[x] = 207;		
								}else{			
// this case will runout based on TurndownDemandRate											
											
											
									outTime[y] = inTime[x] + inDuration[x] + EndInventory/(TurndownDemandRate*numTrains);		
											
									outDuration[y-1] = outTime[y] - outTime[y-1];		
									outProd[y] = 0;		
									DisCapEx[y] = 0;		
									RunOut[y] = 1;		
									EmptyOnD[y] = 0;		
											
											
											
									y=y+1;		
											
									EndInv[x] = 0;		
									LogicPath[x] = 208;		
								}			
							}				
						}					
					}else{						
// This is a very unusual case -- There have already been two aborted TurndownTime events											
// Without further checking this case will be simplified for immediate turndown as appropriate											
						if((GenRate[x] -TurndownLimit+1e-9)>0)  {					
											
							outTime[y] = inTime[x];				
											
							outDuration[y-1] = outTime[y] - outTime[y-1];				
							outProd[y] = GenRate[x];				
							DisCapEx[y] = 0;				
							RunOut[y] = 0;				
							EmptyOnD[y] = 0;				
											
											
											
							y=y+1;				
											
							EndInv[x] = StartInventory;				
							LogicPath[x] = 209;				
						}else{					
// there will still be demand since GenRate is below TurndownLimit											
// after setting a new set of conditions for wid entry the EndInventory needs to be redetermined											
											
							outTime[y] = inTime[x];				
											
							outDuration[y-1] = outTime[y] - outTime[y-1];				
							outProd[y] = TurndownLimit;				
							DisCapEx[y] = 0;				
							RunOut[y] = 0;				
							EmptyOnD[y] = 0;				
											
											
											
							y=y+1;				
											
// now re-determine the EndInventory based on remaining Duration											
// Determine EndInventory											
							TurndownDemandRate = TurndownLimit-GenRate[x];				
							EndInventory = StartInventory-TurndownDemandRate*numTrains*inDuration[x];				
							if(EndInventory > 0)  {				
// successful backup from inventory, no further change to ProdRate											
								EndInv[x] = EndInventory;			
								LogicPath[x] = 210;			
							}else{				
// This case has been a run out at TurndownLimit											
											
											
								outTime[y] = inTime[x] + inDuration[x] + EndInventory/(TurndownDemandRate*numTrains);			
											
								outDuration[y-1] = outTime[y] - outTime[y-1];			
								outProd[y] = 0;			
								DisCapEx[y] = 0;			
								RunOut[y] = 1;			
								EmptyOnD[y] = 0;			
											
											
											
								y=y+1;			
											
								EndInv[x] = 0;			
								LogicPath[x] = 211;			
							}				
						}					
					}						
				}else{							
// GenRate must now be greater than or equal to LastProdRate, but less than full											
// there should be no wid entry if ProdRate is unchanged (as successful backup)											
					if(GenRate[x]-TurndownLimit+1e-9>0)  {						
//  the case where LastProdRate was zero, but system is still below TurndownLimit has been excluded											
						if(LastProdRate==GenRate[x])  {					
// this is an end of the train state											
							EndInv[x] = StartInventory;				
							LogicPath[x] =212;				
						}else{					
// GenRate must be > LastProdRate, so there is a new ProdRate entry in wid											
											
							outTime[y] = inTime[x];				
											
							outDuration[y-1] = outTime[y] - outTime[y-1];				
							outProd[y] = GenRate[x];				
							DisCapEx[y] = 0;				
							RunOut[y] = 0;				
							EmptyOnD[y] = 0;				
											
											
											
							y=y+1;				
// this is an end of the train state											
							EndInv[x] = StartInventory;				
							LogicPath[x] =213;				
						}					
											
					}else{						
// this is an end of the train state											
						EndInv[x] = StartInventory;					
						LogicPath[x] =214;					
					}						
											
				}							
			}								
		} 									
	if(y==insize)  break;										
// close the loop											
 	}										
// Make a final entry for last duration in wid											
	outDuration[y-1]=inTime[insize-1]+ inDuration[insize-1]-outTime[y-1];										
											
	Rcpp::DataFrame NDF1 =										
	Rcpp::DataFrame::create(Rcpp::Named("Time")=outTime,										
	                  Rcpp::Named("Duration")=outDuration,										
	                  Rcpp::Named("ProdRate")=outProd,										
	                  Rcpp::Named("DisCapEx")=DisCapEx,										
	                  Rcpp::Named("RunOut")=RunOut,										
	                  Rcpp::Named("EmptyOnD")=EmptyOnD);										
											
	Rcpp::DataFrame NDF2 =										
	Rcpp::DataFrame::create(										
	                  Rcpp::Named("EndInv")=EndInv,										
	                  Rcpp::Named("LogicPath")=LogicPath);										
	Rcpp::List L=										
	Rcpp::List::create(NDF1,NDF2);										
											
	return(L);										

//	'
}
