RCPP_MODULE(simon){
  using namespace Rcpp ;
                      
	class_<SimonDesign>( "SimonDesign" )
	
	  // expose the default constructor
	  .constructor()   
    
    .property("maxN", &SimonDesign::getMaxN, "Returns the set maximal value for \"n\"")
    .property("alpha", &SimonDesign::getAlpha, &SimonDesign::setAlpha, "Maximum value for type I error")
    .property("beta", &SimonDesign::getBeta, &SimonDesign::setBeta, "Maximum value for type II error")
    .property("p0", &SimonDesign::getP0, &SimonDesign::setP0, "The probability of an event under the null hypothesis")
    .property("p1", &SimonDesign::getP1, &SimonDesign::setP1, "The probability of an event under the alternative hypothesis")
    .property("miniMaxPos", &SimonDesign::getMinimaxPos, "Returns the ID of the minimax design under the found designs.")
    .property("optimalPos", &SimonDesign::getOptimalPos, "Returns the ID of the optimal design under the found designs.")
	  
    .method( "getSolutionCount", &SimonDesign::getSolutionCount, "Returns the number of found Solutions (after calling \"calculateStudySolutions\")")
    .method( "getMaxN", &SimonDesign::getMaxN,"Returns the set maximal value for \"n\"")
    .method( "setMaxN", &SimonDesign::setMaxN, "Sets maxN (but should be calculated automaticly with \"aproximateMaxN\").")
		.method( "setAlpha", &SimonDesign::setAlpha , "Sets the maximal type I error rate" )
  	.method( "setBeta", &SimonDesign::setBeta     , "Sets the maximla type II error rate" )
    .method( "setP0", &SimonDesign::setP0 , "Sets the response rate for an event under the null hypothesis" )
  	.method( "setP1", &SimonDesign::setP1     , "Sets the response rate for an event under the alternative hypothesis" )
    .method( "aproximateMaxN", &SimonDesign::aproximateMaxN     , "Approximate the maximal value of \"n\"" )
    .method( "calculateStudySolutions", &SimonDesign::calculateStudySolutions, "Searches for possible designs for the parameter set (alpha, beta, p0, p1)")
    .method( "getResultsForR", &SimonDesign::getResultsForR, "Returns all found designs, which were found through calling 'calculateStudySolutions' in an data.frame.")
    .method( "calculateSC", &SimonDesign::calculateSC, "Estimates the effect of (non-)stochastic curtailment")
    .method( "getCurResultForR", &SimonDesign::getCurResultForR, "Returns the estimated effect of (non-)stochastic curtailment if it was simulated through \"calculateSC\"")
    .method( "getConditionalPower", &SimonDesign::getConditionalPower, "Returns the conditional power.")
    .method( "calcAlpha", &SimonDesign::calcAlpha, "Returns the actual alpha level.")
	;
}
