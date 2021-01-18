
#include "exportFile.h"
using namespace Rcpp;
using namespace std;
using namespace arma;

/********************************************
*  Constructors ~ Model initialization      *
********************************************/

//  Initialization of a hidden Markov model with categorical observations
//  It has n - hidden states
//  It has m - possible emissions
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP initHMM(SEXP n, SEXP m)
{
    //  error handling
    try {
        // Model initialization
        HMM hmm(as<unsigned int>(n),
            as<unsigned int>(m));                        
        //  Return Model parameters to R environment as a list instead of an object  
        return hmm.toList();        
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("C++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);  
}

//  Initialization of a hidden Markov model with continuous observations
//  It has n - hidden states and m - dimensions
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::export]]
SEXP initGHMM(SEXP n, SEXP m)
{
    //  error handling
    try {
        // Model initialization
        MultiGHMM hmm(as<unsigned int>(n), as<unsigned int>(m));            
        //  Return Model parameters to R environment as a list instead of an object  
        return hmm.toList();       
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);  
}

//  Initialization of a hidden Markov model with categorical observations
//  It has n - hidden states
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::export]]
SEXP initPHMM(SEXP n)
{
    //  error handling
    try {
        // Model initialization
        HMMpoisson hmm(as<unsigned int>(n));    
        //  Return Model parameters to R environment as a list instead of an object
        return hmm.toList();       
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);  
}

//  Function used to verify that all the parameters fulfill the model requirements
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP verifyModel(SEXP model)
{       
    //  error handling
    try {
        //  We get what kind of observation the model uses
        List init_hmm(model);
        string model = as<string>(init_hmm["Model"]); 
        
        //  We init the model        

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }                
        else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }    
        else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }     
        else
            Rf_error("That model is not supported.");

    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);            
}

/********************************************
*                  Setters                  *
********************************************/

//  Set names: State names and Emission names
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP setNames(SEXP hmm, SEXP names)
{
    //  error handling
    try {
        //  We get what kind of observation the model uses
        List init_names(names);
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]); 

        //  We init the model                   
        if(model == "HMM")
        {
            HMM newHMM( init_names["StateNames"] , init_names["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }        
        else if(model == "PHMM")   
        {
            HMMpoisson newHMM( init_names["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }
        else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_names["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }                 
        else
            Rf_error("That model is not supported.");

    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);                  
}

//  Set parameters: A: Transition matrix. B: Emission matrix. Pi: Initial probability vector
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP setParameters(SEXP hmm, SEXP params)
{  
    //  error handling  
    try {
        //  We get what kind of observation the model uses
        List init_params(params);
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);     

        //  We init the model
        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_params["A"], init_params["B"], init_params["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }                                
        else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_params["A"], init_params["B"], init_params["Pi"]);
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        } 
        else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_params["A"]), as<mat>(init_params["Mu"]), as<cube>(init_params["Sigma"]), as<rowvec>(init_params["Pi"]));
            //  Return Model parameters to R environment as a list instead of an object
            return newHMM.toList();
        }           
        else
            Rf_error("That model is not supported.");
                               
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);     
}

/********************************************
*           Sequence Evaluation             *
********************************************/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP evaluation(SEXP hmm, SEXP sequence, SEXP method)
{    
    //  error handling
    try {
        //  We get which evauation algorithm is to be used
        char init_method = as<char>(method);
        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]); 

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return wrap(newHMM.evaluation(sequence, init_method));
        }
        else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return wrap(newHMM.evaluation(sequence, init_method));            
        }
        else if(model == "GHMM")   
        {
            //  We init the model
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            return wrap(newHMM.evaluation(as<mat>(sequence), init_method));
        }
        else
            Rf_error("That model is not supported.");
        
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);  
}

/********************************************
*          Sequence Decodification          *
********************************************/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP viterbi(SEXP hmm, SEXP sequence)
{    
    //  error handling
    try {
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.viterbi(sequence);
        }else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.viterbi(sequence);
        }
        else if(model == "GHMM")   
        {
            //  We init the model
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            return newHMM.viterbi(as<mat>(sequence));
        }
        else
            Rf_error("That model is not supported.");

    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);     
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP forwardBackward(SEXP hmm, SEXP sequence)
{    
    //  error handling
    try {
        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.forwardBackward(sequence);
        }else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.forwardBackward(sequence);
        }
        else if(model == "GHMM")   
        {
            //  We init the model
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            return newHMM.forwardBackward(as<mat>(sequence) );
        }
        else
            Rf_error("That model is not supported.");

    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  Return Model parameters to R environment as a list instead of an object
    return wrap(NA_REAL);      
}

/********************************************
*                 Learning                  *
********************************************/

//  Function used to determine the loglikelihood of a sequence set given a specific model
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP loglikelihood(SEXP hmm, SEXP sequences)
{    
    //  error handling
    try {     
        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return wrap(newHMM.loglikelihood(sequences));
        }else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return wrap(newHMM.loglikelihood(sequences));
        }
        else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            return wrap(newHMM.loglikelihood(as<cube>(sequences)) );
        }
        else
            Rf_error("That model is not supported.");

    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);    
}

//  Function used for parameter estimation given two or more observation sequences
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 
SEXP learnEM(SEXP hmm, SEXP sequences, SEXP iter, SEXP delta, SEXP pseudo, SEXP print)
{
    //  error handling
    try {
        //  How many iterations until the function stops
        unsigned short int init_iter = as<unsigned short int>(iter);
        //  Convergence parameter
        double init_delta = as<double>(delta);
        //  Pseudo counts used if there are not enough observations of a parameter
        unsigned char init_pseudo = as<unsigned char>(pseudo);
        //  Value that determine if the method prints the error value in each iteration 
        bool init_print = as<bool>(print);

        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {            
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);        
            newHMM.learnEM(sequences, init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList();
        }else if(model == "PHMM")
        {                      
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            newHMM.learnEM(sequences, init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList();        
        }else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            newHMM.learnEM(as<cube>(sequences), init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList(); 
        }
        else
            Rf_error("That model is not supported.");
         
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);    
}

//  Function used for parameter estimation given only one observation sequence
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP learnBW(SEXP hmm, SEXP sequences, SEXP iter, SEXP delta, SEXP pseudo, SEXP print)
{
    //  error handling
    try {
        //  How many iterations until the function stops
        unsigned short int init_iter = as<unsigned short int>(iter);
        //  Convergence parameter
        double init_delta = as<double>(delta);
        //  Pseudo counts used if there are not enough observations of a parameter
        unsigned char init_pseudo = as<unsigned char>(pseudo);
        //  Value that determine if the method prints the error value in each iteration
        bool init_print = as<bool>(print);

        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {            
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);        
            newHMM.learnBW(sequences, init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList();
        }else if(model == "PHMM")
        {                      
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            newHMM.learnBW(sequences, init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList();         
        }else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            newHMM.learnBW(as<mat>(sequences), init_iter, init_delta, init_pseudo, init_print);
            return newHMM.toList(); 
        }else
            Rf_error("That model is not supported.");
         
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);    
}


/********************************************
*                Simulation                 *
********************************************/

// Given a specific model, it generates an observation sequence and the hidden states that originated them.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 
SEXP generateObservations(SEXP hmm, SEXP length)
{
    //  error handling
    try {
        //  Length of the sequence
        unsigned short int init_length = as<unsigned short int>(length);

        //  We get what kind of observation the model uses
        List init_hmm(hmm);
        string model = as<string>(init_hmm["Model"]);        

        if(model == "HMM")
        {
            HMM newHMM( init_hmm["StateNames"] , init_hmm["ObservationNames"], init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.generateObservations(init_length);
        }else if(model == "PHMM")
        {
            HMMpoisson newHMM( init_hmm["StateNames"] , init_hmm["A"], init_hmm["B"], init_hmm["Pi"]);
            return newHMM.generateObservations(init_length);
        }else if(model == "GHMM")   
        {
            MultiGHMM newHMM( init_hmm["StateNames"] , as<mat>(init_hmm["A"]), as<mat>(init_hmm["Mu"]),as<cube>(init_hmm["Sigma"]), as<rowvec>(init_hmm["Pi"]));
            return newHMM.generateObservations(init_length);
        }else
            Rf_error("That model is not supported.");
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    Rf_error("c++ exception (unknown reason)"); 
    }
    //  A necessary return for the function, although it is never reached
    return wrap(NA_REAL);     
}


