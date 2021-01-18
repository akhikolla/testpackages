RCPP_MODULE(sub1){
  using namespace Rcpp ;
                      
  class_<Sub1Design>( "Sub1Design" )
	
	  // expose the default constructor
	  .constructor()   
    
    .property("miniMaxPos", &Sub1Design::getMinimaxPos, "Returns the position of the minimax design under the found designs.")
    .property("optimalPos", &Sub1Design::getOptimalPos, "Returns the position of the optimal design under the found designs.")
	  
		.method( "setAlpha", &Sub1Design::setAlpha , "Set the maximal type I error rate" ) 
    .method( "getAlpha", &Sub1Design::getAlpha , "Returns the maximal type I error rate" )
  	.method( "setBeta", &Sub1Design::setBeta     , "Set the maximla type II error rate" )
    .method( "setPc0", &Sub1Design::setPc0 , "Set the response rate for an event in the subset endpoint under the null hypothesis" )
    .method( "setPt0", &Sub1Design::setPt0 , "Set the response rate for an event in the superset endpoint under the null hypothesis" )
    .method( "setPc1", &Sub1Design::setPc1 , "Set the response rate for an event in the subset endpoint under the alternative hypothesis" )
    .method( "setPt1", &Sub1Design::setPt1 , "Set the response rate for an event in the superset endpoint under the alternative hypothesis" )
    .method( "aproximateMaxN", &Sub1Design::aproximateMaxN     , "Approximates and sets the maximal value for \"n\"" )  
    .method( "calculateStudySolutions", &Sub1Design::calculateStudySolutions, "Searches for possible designs for the parameter set (alpha, beta, pc0, pt0, pc1, pt1)")
    .method( "getResultsForR", &Sub1Design::getResultsForR, "Returns all found designs, which were found through calling 'calculateStudySolutions' in an data.frame.")
    .method( "binsum", &Sub1Design::binsum, "Implements the binomial distribution function.")
    .method( "get_p_exact", &Sub1Design::get_p_exact, "Returns the exact p value.")
    .method( "calculateSC", &Sub1Design::calculateSC, "Estimates the effect of (non-)stochastic curtailment.")
    .method( "getSolutionCount", &Sub1Design::getSolutionCount, "Returns the number of found designs/solutions (after calling \"calculateStudySolutions\").")    
    .method( "get_conditionalPower", &Sub1Design::get_conditionalPower, "Returns the conditional power.")
	;
}
