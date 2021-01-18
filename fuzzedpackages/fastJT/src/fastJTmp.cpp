					/********************************/
					/*  Author: Jiaxing Lin 	    */
					/*  Email: Jiaxing.Lin@duke.edu */
					/*  Joncheere-Terpstra Test     */
					/*  parallel version with openMP*/
					/*	03/28/2016                  */
					/********************************/
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include "pairCompareMP.h"
#include "JTscoreMP.h"

using namespace Rcpp;
using namespace std;

    /*******************************************************************************
        This is a c++ implemented Jonckheere-Terpstra test routine. The routine con-
        tains two major loops for both markers and gene snips.
    *******************************************************************************/

// [[Rcpp::export]]
Rcpp::List fastJTmp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, 
					bool outTopNFlag = true, int numThreads = 8, int outItemNo = 15,
					bool standardized = true) 
{
	
	long  	gSnips 		= Y.ncol();	// number of gene snips
	long  	numPatWithNA= X.nrow();	// number of patients
	int 	numMarker 	= X.ncol();	// number of makers(genes) same as geneSnips
   
	Rcpp::NumericMatrix gSnipID(outItemNo,numMarker); 
	
	if(outItemNo > gSnips)
		outItemNo = gSnips;
	if(outTopNFlag == false)
		outItemNo = gSnips;
	
    Rcpp::NumericMatrix JTJStar(outItemNo, numMarker);
    Rcpp::NumericMatrix jStatistic(outItemNo, numMarker);

	//List result that return to R. 
	Rcpp::List JTtestRes; 	
			
	// Outer loop for markers.
	#ifdef _OPENMP
	omp_set_num_threads(numThreads); // define number of threads.	
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for 
	#endif
	for(int curMarker = 0; curMarker < numMarker; curMarker++){ 

		int 	NA_count 	= 0;
		int 	gTypes 		= 0;
		int 	gTypesMax	= 0;
		int		outcount	= 0;   		// output size control

		std::list<JTscoreMP> jTScoreTopN;
		std::list<JTscoreMP>::iterator it;
		
		bool tieState;

		long itor = 0; 	//iterator for asign values to the index vector
		
		long long nonEqCount = 0; 	// count for the unequal elements
		long long eqCount    = 0; 	// count for the tie elements.
		long tmpItor  		 = 0; 	// tempary itorator
	
        double mean = 0.00;
		double var  = 0.00;
	
		double termA, termA1, termA2, termA3;
		double termB, termB1, termB2;
		double termC, termC1, termC2;
					
		//(1) obtain the tie informtion: number of tie group, and number of patients
		//    in each tied group. Here we use unordered_map for this purpose.
		// 	  Tie group information only depends on the curMarker, therefore handled
		//    in the outer loop.
		std::unordered_map<double, long> tieGroups_pre; // patient number for tie groups
		for(long pid = 0; pid < numPatWithNA; pid++){
			if(!R_IsNA(X(pid,curMarker)))
			tieGroups_pre[X(pid,curMarker)]
			= tieGroups_pre[X(pid, curMarker)] + 1;
        }
			
		for(long curGSnip = 0; curGSnip < gSnips; curGSnip++ ){	
				
			// Inner loop for gene snips
			//(2) determine the number of a gene types for current gene snip this proces-
			//    sing can be expensive, since it is in the geno snip loop that looped m-
			//    any times. My current testing shows that this part cost 2% of the over-
			//    all computation cost.
			NA_count 	= 0;
			gTypes  	= 0;

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

			// adjust the patient number to exclude those with NA for their SNP type
			int numPat = numPatWithNA-NA_count;

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
 
			//(3) determine number of patients associated with each type of geno for the 
			//    current gene snip 
			vector<long> numPatGType(gTypes,0);
			for(long i = 0; i < numPatWithNA; i++){
				if(!R_IsNA(Y(i,curGSnip))&& !R_IsNA(X(i,curMarker)) )
					numPatGType[gTypesMap[Y(i, curGSnip)]]++;
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
						if(gTypesMap[Y(j, curGSnip)] == genoType){ 
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

			nonEqCount	= 0; 	// reset count for the unequal elements
			eqCount		= 0; 	// reset count for the tie elements.
			tmpItor		= 0; 	// reset tempary itorator used in the for loops.
    	
			// perform the comprision
			vector<long long> pairCom(2,0); // vector for non equal and equal count.
			try
			{
				for(int genoType = 0; genoType < gTypes-1; genoType++){ 
					// loop through each geno type
					pairCom = pairCompareMP(X, curMarker,tmpItor, 
							tmpItor + numPatGType[genoType]-1, 
							tmpItor + numPatGType[genoType], 
							numPat-1, gTypePatInd);
					nonEqCount += pairCom[0];
					eqCount += pairCom[1];
					tmpItor = tmpItor + numPatGType[genoType];
				}
			}catch(std::bad_alloc &){}

			// mean value ,and varience of J
            mean = 0.00;
			var  = 0.00;
			long    n   = numPat;//for concise expression

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
				var = (termA1 - termA2)/72; }else{
					// If there is tie, the tie case JT test score expression is used.
				termA2 = 0.00;
				termB1 = 0.00; 
				termC1 = 0.00;

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
        				+termB /(36*n*(n-1)*(n-2)) 
        				+termC /(8*n*(n-1));
			}
			
			// Compute the statistics Jstar.
			double Jstar = 0.00;
  			Jstar =  ( nonEqCount + 0.5*eqCount - mean)/sqrt(var);
				
			if(outTopNFlag == false)
			{
				JTJStar(curGSnip, curMarker)   = Jstar;
				jStatistic(curGSnip, curMarker) = nonEqCount + 0.5*eqCount; 		
			}else{
			// Store the Jstat value and in euqal count for current marker and current g-
			// ene snip
			JTscoreMP js(Jstar, nonEqCount+0.5*eqCount, curGSnip);		
			#ifdef _OPENMP
			#pragma omp critical
			#endif
			{            
				if (outcount < outItemNo){	
                	// sort and store the top N results JT test in a list.
                	it = jTScoreTopN.begin();
                	//while(abs(Jstar) <= abs((*it).getJStar()) && it!= jTScoreTopN.end())
                   	//	it++;
					//************************//
					// fix to remove the stack overflow error
					while(it!= jTScoreTopN.end())
					{
						if(abs(Jstar) <= abs((*it).getJStar()))
							it++;
						else
							break;
					}
                	jTScoreTopN.insert(it,js);
					outcount++;
            	}else{
                	// Dynamicly update the top N significant results list with sorting.
                	if(abs(Jstar) <= abs(jTScoreTopN.back().getJStar())){} 
			        	  // smaller than smallest element, no update needed.
                	else{
						it = jTScoreTopN.begin();
                		while(abs(Jstar) <= abs((*it).getJStar()) && it!= jTScoreTopN.end())
                	    	it++;               // find the insert position
               	 		//Insert qulified JT test results for current gene snip and marker
						jTScoreTopN.insert(it,js);
                		jTScoreTopN.pop_back();
					}
            	}
			}
			}
			// end of selecting top N significant JT score for out put.
		} 	
		// strip the values from the JTscoreMP class to numericvectors
    	if (outTopNFlag != false)
		{	int i = 0;
    		for(it = jTScoreTopN.begin(); it !=  jTScoreTopN.end(); it++)
    		{
        		JTJStar(i,curMarker)    =(*it).getJStar();
        		jStatistic(i,curMarker) =(*it).getNonEqCount();
        		gSnipID(i,curMarker)	=(*it).getGSnipID();
				i++ ;
    		}
		}

    }

    // Import the results to a master R list object to return to R.
	if(standardized)
    	JTtestRes["J"]  = JTJStar;
    else
		JTtestRes["J"]  = jStatistic;
	if(outTopNFlag!= false)
		JTtestRes["XIDs"] = gSnipID;		
//	JTtestRes["p.value"]	= JTJStar;
	return(JTtestRes);
}

