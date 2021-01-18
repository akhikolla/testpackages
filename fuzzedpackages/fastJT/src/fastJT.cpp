					/********************************/
					/*  Author: Jiaxing Lin 	    */
					/*  Email: Jiaxing.Lin@duke.edu */
					/*  Joncheere-Terpstra Test     */
					/*  parallel version with openMP*/
					/*	03/22/2016                  */
					/********************************/

/*
		     Title: Rcpp routine for Jonckheere Terpstra test (with ties)
					parallel version based on openMP
*/

#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include "pairCompare.h"
#include "JTscore.h"

using namespace Rcpp;
using namespace std;

    /*******************************************************************************
        This is a c++ implemented Jonckheere-Terpstra test routine. The routine con-
        tains two major loops for both markers and gene snips.
    *******************************************************************************/

// [[Rcpp::export]]
Rcpp::List fastJT(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, 
                  bool outTopNFlag = true, int outItemNo = 15, bool standardized = true ) 
{	
	int 	gTypes 		= 0;
	long  	gSnips 		= Y.ncol();	// number of gene snips
	long  	numPatWithNA= X.nrow();	// number of patients
	int 	numMarker 	= X.ncol();	// number of makers(genes) same as geneSnips

	Rcpp::NumericMatrix gSnipID(outItemNo,numMarker);
    
	if(outItemNo > gSnips)
		outItemNo = gSnips;
	if(outTopNFlag == false)
		outItemNo = gSnips;

	Rcpp::NumericMatrix JTJStar(outItemNo,numMarker);
    Rcpp::NumericMatrix jStatistic(outItemNo,numMarker);

	//List result that return to R. 
	Rcpp::List JTtestRes; 	
	std::list<JTscore>::iterator it;
	//Variable used in the loop defined here:
	double termA1 = 0.00;
	double termA2 = 0.00;
	double termA3 = 0.00;

	double termB1 = 0.00;
	double termB2 = 0.00;

	double termC1 = 0.00;
	double termC2 = 0.00;

	double termA  = 0.00; 
	double termB  = 0.00;
	double termC  = 0.00;

	bool tieState;
	long itor = 0;

	long long nonEqCount = 0;
	long long eqCount    = 0;
	long tmpItor         = 0;

	double mean = 0.00;
	double var  = 0.00;
	double    n   = 0;
	
	double Jstar = 0.00;		
	
	// Outer loop for markers.
	for(int curMarker = 0; curMarker < numMarker; curMarker++){ 
		
		//(1) obtain the tie informtion: number of tie group, and number of patients
		//    in each tied group. Here we use unordered_map for this purpose.
		// 	  Tie group information only depends on the curMarker, therefore handled
		//    in the outer loop. We precompute the ties and them for each SNP, we w-
		//    ill kick out those SNPs are NA. 
		std::unordered_map<double, long> tieGroups_pre; // patient number for tie groups
		for(long pid = 0; pid < numPatWithNA; pid++){
			if(!R_IsNA(X(pid, curMarker)))		
				tieGroups_pre[X(pid,curMarker)]
				= tieGroups_pre[X(pid,curMarker)] + 1;
		}
		
		// Output size contorl for each marker
		int outcount = 1;	
        std::list<JTscore> jTScoreTopN;// list to store top N scores for each marker

		// Inner loop for gene snips
		for(long curGSnip = 0; curGSnip < gSnips; curGSnip++ ){	
			
			int NA_count =0;	
			gTypes 	= 0;
			
			std::unordered_map<double, long> tieGroups = tieGroups_pre;	
			std::set<double> gTypesSet;//use this map the check how many type geno.

			for(int i =0; i < numPatWithNA; i++){
				if(!R_IsNA(Y(i, curGSnip))&&!R_IsNA(X(i,curMarker))){
					if (gTypesSet.find(Y(i,curGSnip)) == gTypesSet.end()){
						gTypesSet.insert(Y(i,curGSnip));
						gTypes ++;
      				}
				}
				else{
					if(!R_IsNA(X(i,curMarker)))
						tieGroups[X(i,curMarker)] -= 1;
					NA_count ++; 
				}
			}
			
			std::map<double, int> gTypesMap;
			std::set<double>::iterator iter;

			// sort the SNP id and put into a map.
			int SNPmap = 0;
			for(iter = gTypesSet.begin(); iter!= gTypesSet.end(); ++iter)
			{
				gTypesMap[*iter] = SNPmap;
				SNPmap++;
			}
			
			// adjust the patient number to exclude those with NA for their SNP type and marker level. 
			int	numPat = numPatWithNA-NA_count;
		
			// Terms in computing the varience that base on the tie. For all the gene sn-
			// ips, this will be the same. It only depend on the marker. Therefore,
			// compute in the outer loop.
			termA3 = 0.00;
			termB2 = 0.00;
			termC2 = 0.00;

			for(const auto &tieItor:tieGroups){
				if(tieItor.second > 1){
					termA3 += tieItor.second
							*(tieItor.second - 1)
							*(2*tieItor.second +5);

					termB2 +=  tieItor.second
							*(tieItor.second - 1)
							*(tieItor.second - 2);

					termC2 +=  tieItor.second
							*(tieItor.second - 1);
				}	
			}
	
			// Tie state
			tieState  = termA3 > 0.00? true:false;
 
			//(2) determine the number of a gene types for current gene snip this proces-
			//    sing can be expensive, since it is in the geno snip loop that looped m-
			//    any times. My current testing shows that this part cost 2% of the over-
			//    all computation cost.
			
			
			// only on type of geno is not able for performing JT test.  in the eval-
			// uating this geno snip is skipped.
			if(gTypes < 2)	
			{
				if(!outTopNFlag)
				{
					JTJStar(curGSnip,curMarker)   = NA_REAL;
                	jStatistic(curGSnip,curMarker) = NA_REAL;
				}
				continue;
			}
			//(3) determine number of patients associated with each type of geno for the 
			//    current gene snip 
			vector<long> numPatGType(gTypes,0);
			for(long i = 0; i < numPatWithNA; i++){
				if(!R_IsNA(Y(i,curGSnip))&&!R_IsNA(X(i, curMarker)))
					numPatGType[gTypesMap[Y(i,curGSnip)]]++;
    		}
			
			// Test if certain geno type for this gene snip has unsignificant patient 
			// number, if so this geno snip is skipped.
			//for(unsigned int i = 0; i < numPatGType.size(); i++)
			//{
			//	if(numPatGType[i]/numPat < 0.05||numPatGType[i]/numPat > 0.95)
			//		continue;	
			//}
	        
			//(4) create a index array to specify the patient ID that associated certain
			//    geno type. For a specified geno type, the patient ID is not sorted or 
			//    ordered by any rule. The algorithm need gTypes pass of all the patieni-
			//    ts.The cost is O(gTypes*number of patient). There is one pass algorithm 
			//    that needs a extra map to store and update locate current indexed patie-
			//    nt ID, however very complicated to implement, and the overhead could be 
			//    significant. Only the gTypes is large. Then you will need one pass algo-
			//    rithm.
			vector<long> gTypePatInd(numPat,0);
			itor = 0; 	//iterator for asign values to the index vector
			
			for(int genoType = 0; genoType < gTypes; genoType++){
				for (long j = 0; j < numPatWithNA; j++){
					if(!R_IsNA(X(j,curMarker)) && !R_IsNA(Y(j,curGSnip))){
						if( gTypesMap[Y(j,curGSnip)] == genoType ){ 	
							gTypePatInd[itor] = j;  		//asign the index
							itor++; 		  				//update asigning iterator
        				}
					}
     			}
    		}
    		
			//(5) Scan and implicitly construct for the Mann-Whitney count.
			//	  Here we do not directly perform comparision for the Mann-Whitney count.
			//	  Indeed, the marker level for each geno types are sorted before compari-
			// 	  sion, such that the cost of algorithm complexity is reduced from O(N^2) 
			//    to O(Nlog(N)), 1000 time faster can be obtained with one million patie-
			//    nts. Current computation cost is 1 sec in single core per marker per g-
			//    ene snip for one million patient.

			nonEqCount		= 0; 	// count for the unequal elements
			eqCount			= 0; 	// count for the tie elements.
			tmpItor		 	= 0; 	// tempary itorator used in the for loops.
    	
			// Testing of only passing the current marker column of X_copy to the compa-
			// ring subroutine.
			
			// perform the comprision
			vector<long long> pairCom(2,0); // vector for non equal and equal count.
			for(int genoType = 0; genoType < gTypes-1; genoType++){ 
					// loop through each geno type
					pairCom = pairCompare(X, curMarker,tmpItor, 
							tmpItor + numPatGType[genoType]-1, 
							tmpItor + numPatGType[genoType], 
							numPat-1, gTypePatInd);
					nonEqCount += pairCom[0];
					eqCount += pairCom[1];
					tmpItor = tmpItor + numPatGType[genoType];
			}
			
			// mean value ,and varience of J
            mean = 0.00;
			var  = 0.00;
			n   = numPat;//for concise expression

            mean =  numPat*numPat -
                    inner_product(
                        numPatGType.begin(),
                        numPatGType.end(),
                        numPatGType.begin(),0);
            mean = mean /4.0;
					
			// If there is no tie, then the no tie JT test score expression is  used.
			if(!tieState)
			{
				termA1 = n*n*(2*n+3);
				termA2 = 0.00;
				for(int i = 0; i < gTypes; i++){
					termA2 +=  numPatGType[i]
							*  numPatGType[i]
							* (2*numPatGType[i]+3);
				}
				var = (termA1 - termA2)/72; 
			}else{
				// If there is tie, the tie case JT test score expression is used.
				//termA, termB, termC;
				termA2 = 0.00, termB1 = 0.00, termC1 = 0.00;

				termA1 = n*(n-1)*(2*n+5);
				for(int i = 0; i < gTypes; i++){
					termA2 +=  numPatGType[i]
						* (numPatGType[i]-1)
						* (2*numPatGType[i]+5);
				
					termB1 +=  numPatGType[i]
						* (numPatGType[i] - 1)
						* (numPatGType[i] - 2);
      			
					termC1 +=  numPatGType[i]
      					* (numPatGType[i] - 1);
    			}
			
    			termA = termA1 - termA2 - termA3;
				termB = termB1*termB2;
				termC = termC1*termC2;
			
				// Assemble for variance of J from term A, B, C.	
				var   =  termA /72.00
        				+termB /(36.00*n*(n-1)*(n-2))
        				+termC /(8.00*n*(n-1));
			}
			
			// Compute the statistics Jstar.
			Jstar = ( nonEqCount + 0.5*eqCount - mean)/sqrt(var);
			
			// Store the Jstat value and in euqal count for current marker and current g-
			// ene snip
//			JTJStar[curGSnip][curMarker] = Jstar;
//			jStatistic[curGSnip][curMarker] = nonEqCount;	
	
			if(outTopNFlag==false)
			{
				JTJStar(curGSnip,curMarker)   = Jstar;
				jStatistic(curGSnip,curMarker) = nonEqCount + 0.5*eqCount;
				continue;
			}
			JTscore js(Jstar, nonEqCount+0.5*eqCount, curGSnip);			
           	if (outcount <= outItemNo){ // populate up the result list before droping out
				outcount++;
               	it = jTScoreTopN.begin(); // reset the index to the beginning of the list
				if (jTScoreTopN.empty())
				{
					jTScoreTopN.insert(it,js); //when result list is empty, skip the find
										       //insert position step.
					continue;
				}else{
               		//while((abs(Jstar) <= abs((*it).getJStar())) && (it!= jTScoreTopN.end()))
                   	//{it++;}   			// find the insert postion
					//*****************************//
					// fix to remove the stack overflow error
					while(it!= jTScoreTopN.end())
					{
						if(abs(Jstar) <= abs((*it).getJStar()))
							it++;
						else
							break;
					}
					jTScoreTopN.insert(it,js); //insert the new result item
					continue;
				}
			}
					
			if(abs(Jstar) <= abs(jTScoreTopN.back().getJStar())) // updating the result list.
       		{	
				continue;          // new result item is smaller than smallest element in the result list
							   // no update needed for this new result item
       		}else{ 
				it = jTScoreTopN.begin(); // reset the index to the beginning of the list
				while(abs(Jstar) <= abs((*it).getJStar()))
				{it++;}               // find the insert position
				jTScoreTopN.insert(it,js);  // insort the new result item
				jTScoreTopN.pop_back();		// kick out the smallest element
			}
			 	
	    }		
	   		
		if(outTopNFlag == false)
			continue;	
		// Strip the result from the class to the Rcpp vector.
   		int i = 0;
   		for(it = jTScoreTopN.begin(); it !=  jTScoreTopN.end(); it++)
   		{
       		JTJStar(i,curMarker)      = (*it).getJStar();
       		jStatistic(i,curMarker)   = (*it).getNonEqCount();
 			gSnipID(i,curMarker)      = (*it).getGSnipID();      		
			i++ ;
   		}
	}
    // Import the results to a master R list object to return to R.
    if(standardized)
		JTtestRes["J"]  	= JTJStar;
	else
    	JTtestRes["J"]  	= jStatistic;
	
	if(outTopNFlag == true)
		JTtestRes["XIDs"]= gSnipID;
	
	return(JTtestRes);
}

