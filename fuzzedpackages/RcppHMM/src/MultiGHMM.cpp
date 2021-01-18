// included dependencies
#include "MultiGHMM.h"

//System libraries
#include <iostream>
#include <algorithm>    /* find */
#include <math.h>       /* log */

//namespaces

using namespace std;
using namespace Rcpp;
using namespace arma;

//--------------------------------------------------------------------------
// CONSTRUCTORS & DESTRUCTOR:
//--------------------------------------------------------------------------

//  First constructor
MultiGHMM::MultiGHMM(unsigned short int  numberStates, unsigned short int numberObservations)
{
    //Validate the values
    if(numberStates < 2)
        Rf_error("The number of states must be bigger or equal to 2.");

    //  Set known values 
    m_M = numberObservations;
    m_N = numberStates;    
    m_StateNames = CharacterVector(m_N);

    // Memory allocation for parameters
    m_Pi    = rowvec(m_N,fill::zeros);
    m_A     = mat(m_N,m_N,fill::zeros);        
    m_mu    = mat(m_M,m_N,fill::zeros);
    m_sigma = cube(m_M,m_M,m_N,fill::zeros);

    /**********************************************************/
    //  Proposed state names
    for(int i = 1; i <= m_N; i++ )
        m_StateNames[i-1] = "x" + to_string(i);         

    //  Parameter random initialization               
    randomInit(-10,10);
}

//  Second constructor used for model validation 
MultiGHMM::MultiGHMM(CharacterVector stateNames, mat A, mat mu, cube sigma, rowvec Pi)
{
    //Validate the size
    if(stateNames.size() < 2) 
        Rf_error("The number of states must be bigger or equal to 2.");           
    
    //Validate the size
    if(stateNames.size() != A.n_rows|| stateNames.size() != A.n_cols)            
        Rf_error("The number of states must be the same as the transition matrix column and row size");

    //Validate the size
    if(sigma.n_cols != mu.n_rows || sigma.n_rows != mu.n_rows || sigma.n_slices != mu.n_cols)
        Rf_error("The size of the rows and columns in the Sigma matrix must be the same as the rows size in Mu. Also, the depth of Sigma must be equal to the cols size in Mu");        
    
    //Validate the size
    if(stateNames.size() != Pi.size())            
        Rf_error("The number of states must be the same as the initial probability vector size");

    //Validate the parameters
    if(verifyVector(Pi) == false)
        Rf_error("The initial probability vector is not normalized");
    
    //Validate the parameters
    if(verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized");

     //Validate the parameters
     bool flag = true;
    for(unsigned int i = 0; i < sigma.n_slices; i++)
    {
        mat slice = sigma.slice(i);
        if(isPositiveDefinite(slice) == false)
            Rf_error("All the Sigma slices must be positive definite."); 

        if(flag == true && mu.n_rows > 1 && det(slice) < 1.0 / pow((2 * 3.14159), mu.n_rows) ) 
        {
            Rf_warning("It is recommended to have a covariance matrix with a determinant bigger than 1/ ((2*PI)^k) .");                  
            flag = false;
        }else if( flag == true && mu.n_rows == 1 && slice(0,0)  < 1.0 / (2 * 3.14159))
        {
            Rf_warning("The variance is recommended to be bigger than 1/(2*PI)");
            flag = false;
        }
            
    }  
    
    //  If all the paremeters have been validated, then they are used.
    m_M = mu.n_rows;
    m_N = stateNames.size();
    m_StateNames = stateNames ;
                        
    m_A     = A;  
    m_mu    = mu;
    m_Pi    = Pi;
    m_sigma = sigma;
}

//  Destructor
MultiGHMM::~MultiGHMM(void){}

//--------------------------------------------------------------------------
//  PUBLIC
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//  GETTERS
//--------------------------------------------------------------------------

CharacterVector MultiGHMM::getStateNames(void) const
{
    return m_StateNames;
}

mat MultiGHMM::getA(void) const
{    
    return m_A;
}

mat MultiGHMM::getMu(void) const
{
    return m_mu;
}

cube MultiGHMM::getSigma(void) const
{
    return m_sigma;
}

rowvec MultiGHMM::getPi(void) const
{
    return m_Pi;
}

//--------------------------------------------------------------------------
//  SETTERS
//
//  Every method has an error handler
//--------------------------------------------------------------------------

void MultiGHMM::setStateNames(CharacterVector stateNames)
{
    if( stateNames.size() != m_N)
        Rf_error("The number of state names does not coincide with the one declared.");
    m_StateNames = CharacterVector(clone(stateNames));
}

void MultiGHMM::setA(mat A)
{    
    if(A.n_rows != m_N || A.n_cols != m_N || verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized or the size is wrong");
    m_A = A;
}

void MultiGHMM::setMu(mat mu)
{
    if(mu.n_rows != m_M || mu.n_cols != m_N)
        // Rf_error("The mu columns must be the same as the number of states and the rows must have the same dimensionality of the observation vector");
        Rf_error("The mu matrix size is wrong");
    m_mu = mu;
}

void MultiGHMM::setSigma(cube sigma)
{
    if(sigma.n_rows != m_M || sigma.n_cols != m_M || sigma.n_slices != m_N)
        Rf_error("The covariance matrix size is wrong");

    //Validate the parameters
    for(unsigned int i = 0; i < sigma.n_slices; i++)
        if(isPositiveDefinite(sigma.slice(i)) == false)
            Rf_error("All the Sigma slices must be positive definite.");

    m_sigma = sigma;
}

void MultiGHMM::setPi(rowvec Pi)
{
    if(Pi.size() != m_N || verifyVector(Pi) == false)        
        Rf_error("The initial probability vector is not normalized or the size is wrong");
    m_Pi = Pi;
}

void MultiGHMM::setParameters(mat A, mat mu, cube sigma, rowvec Pi)
{     
    //Validate the size
    if(sigma.n_cols != mu.n_rows || sigma.n_rows != mu.n_rows || sigma.n_slices != mu.n_cols)
        Rf_error("The size of the rows and columns in the Sigma matrix must be the same as the rows size in Mu. Also, the depth of Sigma must be equal to the cols size in Mu");        
    
    //Validate the size
    if( Pi.size() != m_N)            
        Rf_error("The number of states must be the same as the initial probability vector size");

    //Validate the size
    if(A.n_rows != m_N || A.n_cols != m_N)            
        Rf_error("The number of states must be the same as the transition matrix column and row size");

    //Validate the parameters
    if(verifyVector(Pi) == false)
        Rf_error("The initial probability vector is not normalized");
    
    //Validate the parameters
    if(verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized");

    //Validate the parameters
    for(unsigned int i = 0; i < sigma.n_slices; i++)
    {
        mat slice = sigma.slice(i);
        Rcout << slice(0,0) ;
        if(isPositiveDefinite(slice) == false)
            Rf_error("All the Sigma slices must be positive definite.");                 
        if(m_M > 1 && det(slice) < 1/pow((2 * 3.14159), m_M) ) 
            Rf_warning("It is recommended to have a covariance matrix with a determinant bigger than 1/ ((2*PI)^k) .");                  
        else if(m_M == 1 && slice(0,0)  < MINSD)
            Rf_warning("The standard deviation is recommended to be bigger than 1/sqrt(2*PI)");
    }
        
        
    m_M     = mu.n_rows;
    m_A     = A;
    m_mu    = mu;
    m_Pi    = Pi;
    m_sigma = sigma;
}

//--------------------------------------------------------------------------
//  EVALUATION
//--------------------------------------------------------------------------

double MultiGHMM::evaluation(mat sequences, char method)
{
    double eval = 0.0;
    unsigned int i, length;
    
    length = sequences.n_cols;    

    // Memory allocation
    rowvec scaled(length,fill::zeros);
    mat matrix(m_N, length,fill::zeros);
    scaledMatrixM eva = {scaled, matrix};

    //  Selected evaluation method
    switch(method)
    {
        case 'f':
            forwardMatrix( sequences, length , eva );             
            break;
        case 'b':
            backwardMatrix( sequences, length , eva );                
            break;
    }

    //  double variables lose precision if multiplication is used
    //  Therfore, a sum of log values is used
    for(i = 0; i < length ; i++ )
                eval+= log(eva.scaling[i]);
    
    return eval;
}

//--------------------------------------------------------------------------
// DECODING
//--------------------------------------------------------------------------

CharacterVector MultiGHMM::forwardBackward(mat sequences)
{
    unsigned i, j;
    unsigned int length = sequences.n_cols;

    //  We get P(State|Data) for all states and observations
    mat gamma = forwardBackwardGamma(sequences);

    //  Most probable state-path traveled
    IntegerVector best(length);

    //  Temp vector to store a column 
    rowvec temp(m_N, fill::zeros ); 
    for(j = 0; j < length; j++)
    {
        for(i = 0; i < m_N ; i++)
            temp[i] = gamma(i,j);
        best[j] = distance(temp.begin(), max_element(temp.begin(), temp.end()));
    }

    //  Change discrete data into categorical values
    return toName(best, 'S');   
}

CharacterVector MultiGHMM::viterbi(mat sequences)
{
    unsigned int length, i, j, k;    
    
    length = sequences.n_cols;     

    //  Most probable state-path traveled 
    IntegerVector best(length);

    // Memory allocation
    mat phi(m_N, length, fill::zeros);
    mat delta(m_N, length, fill::zeros);        

    // Used to find the max value
    rowvec temp(m_N, fill::zeros);   

    //  log-transform to avoid double precision-loss
    mat A(m_N, m_N, fill::zeros);      
    rowvec Pi(m_N, fill::zeros);

    for(i = 0 ; i < m_N; i++)
    {
        Pi[i] = log(m_Pi[i]);
        for(j = 0; j < m_N; j++)
            A(i,j) = log(m_A(i,j)); 
    }

    //  Init step
    for(i = 0; i < m_N; i++)
        delta(i,0) = Pi[i] + dmvnormSingle(sequences.col(0) , m_mu.col(i), m_sigma.slice(i), true) ;

    //  Recursion step
    for(j = 1; j < length; j++)
        for(i = 0; i < m_N; i++)
        {
            for(k = 0; k < m_N; k++)
                temp[k] = delta(k,j-1) + A(k,i);
            // The auto keyword is simply asking the compiler to deduce the type of the variable from the initialization                
            auto maximum = max_element(temp.begin(), temp.end());             
            delta(i,j) = (*maximum) + dmvnormSingle(sequences.col(j), m_mu.col(i), m_sigma.slice(i), true) ;                    
            phi(i,j) = distance(temp.begin(), maximum);             
        }

    // Termination step
    for(k = 0; k < m_N; k++)
        temp[k] = delta(k,length-1);

    auto maximum = max_element(temp.begin(), temp.end());
    best[length-1] = distance(temp.begin(), maximum);

    for(j = length-1; j > 0; j--)
        best[j-1] = phi(best[j],j);                 
    
    //  Change discrete data into categorical values
    return toName(best, 'S');
}

//--------------------------------------------------------------------------
// LEARNING
//--------------------------------------------------------------------------

//  Loglikelihood of a sequence set
double MultiGHMM::loglikelihood(cube sequences)
{
    double ll = 0.0;
    unsigned int i, seqLen;

    seqLen = sequences.n_slices;

    for(i = 0; i < seqLen; i++)
        ll+=evaluation(sequences.slice(i));
    return ll;
}

//  Parameter estimation using a Expectation Maximization approach for multiple sequences
void MultiGHMM::learnEM(cube sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{
    
    double newLL, error;
    double lastLL = loglikelihood(sequences);
    unsigned int counter = 0; 
    double max = 0.0, min = 0.0;    
    mat seqRow;

    // This flags goes false when a slice of the covariance matrix is not positive definite
    bool EM = true;

    //  We search the min and max value, if a reinitialization is needed 
    auto tmin =  min_element ( sequences.begin(), sequences.end() );
    auto tmax =  max_element ( sequences.begin(), sequences.end() );    
    min = *tmin;                
    max = *tmax;    
    
    //  Parameter estimation
    do{
        //  If the error is nan, it may be a big error. 
        //  A new parameter initialization is recomended
        EM = expectationMaximization(sequences, pseudo);
        newLL = loglikelihood(sequences) ;

        if(std::isnan(newLL) || !EM)
        {
            if(print)
                Rcout << "Convergence error, ll: " << newLL << ", EM: " << EM << " new initialization needed\n";                              
            randomInit(min, max);
            lastLL = loglikelihood(sequences);
            counter++;
            error = 1e10;
            continue;
        }  

        error = fabs(newLL - lastLL);        
        lastLL = newLL;
        counter++;
        
        if(print)
            Rcout << "Iteration: " << counter << " Error: " << error  << "\n";            
    } while(counter < iter && error > delta);   // Convergence criteria
    
    Rcout << "Finished at Iteration: " << counter << " with Error: " << error  << "\n";
}

//  Parameter estimation using a Baum Welch approach for a single sequence
void MultiGHMM::learnBW(mat sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{
    double newLL, error;
    double lastLL = evaluation(sequences);
    unsigned int counter = 0; 
    double max = 0.0, min = 0.0;    

    bool EM = true;

    //  We search the min and max value, if a reinitialization is needed
    auto tmin =  min_element ( sequences.begin(), sequences.end() );
    auto tmax =  max_element ( sequences.begin(), sequences.end() );    
    min = *tmin;            
    max = *tmax;

    //  Parameter estimation
    do{
        //  If the error is nan, it may be a big error. 
        //  A new parameter initialization is recomended
        EM =  BaumWelch(sequences, pseudo);
        newLL = evaluation(sequences);
        if(std::isnan(newLL)  || !EM)
        {
            if(print)
                Rcout << "Convergence error, new initialization needed\n";               
            randomInit(min,max);
            lastLL = evaluation(sequences);
            counter++;
            error = 1e10;
            continue;
        }
        error = fabs(newLL - lastLL);        
        lastLL = newLL;
        counter++;
        
        if(print)
            Rcout << "Iteration: " << counter << " Error: " << error  << "\n";            
    } while(counter < iter && error > delta); // Convergence criteria
    Rcout << "Finished at Iteration: " << counter << " with Error: " << error  << "\n";
}

//--------------------------------------------------------------------------
// SIMULATION
//--------------------------------------------------------------------------

//  Funtion used to generate observations given a model
List MultiGHMM::generateObservations( unsigned short int length)
{

    unsigned int i,j;
    double x;    

    IntegerVector X(length, 0);
    mat Y(m_M, length, fill::zeros);

    //  Used for set.seed compatibility
    RNGScope scope;

    //  Matrix rearrangement to use a uniform distribution to generate the hidden states and its corresponding observation
    mat A(m_N, m_N);     
    rowvec Pi(m_N);

    double tempPi = 0.0, tempA1 = 0.0;

    //  We  fill each value with its corresponding new one
    for(i=0; i < m_N; i++)
    {
        //  We fill first the initial probability vector
        tempPi += m_Pi[i];
        Pi[i] = tempPi;                
        //  Then, we fill the transition matrix
        tempA1 = 0;
        for(j = 0; j < m_N; j++)
        {
            tempA1 += m_A(i,j);
            A(i,j) = tempA1;            
        }                     
    }

    //  Random variable generation based in the rearranged matrices    
    rowvec tempA;

    x = as<double>(runif(1));       
    X[0] = lower_bound (Pi.begin(), Pi.end(), x) - Pi.begin();
    Y.col(0) = rmvnormSingle(m_mu.col(X[0]), m_sigma.slice(X[0])) ; 

    for(j = 1; j < length; j++)
    {
        x = as<double>(runif(1));

        tempA = A.row(X[j-1]);     
        X[j] = lower_bound (tempA.begin(), tempA.end(), x) - tempA.begin();                
        Y.col(j) = rmvnormSingle(m_mu.col(X[j]), m_sigma.slice(X[j])) ;
    }

    //  Returns the hidden state path 'X' and its emissions 'Y'
    return List::create(
                Named("X", toName(X, 'S')),
                Named("Y", Y)
            );

}

//--------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------

//  Funtion to return all the model parameters as an R List
List MultiGHMM::toList(void) const
{
    return List::create(            
            Named("Model", "GHMM" ),
            Named("StateNames", getStateNames() ),
            Named("A",          getA()),
            Named("Mu",         getMu()),
            Named("Sigma",      getSigma()),
            Named("Pi",         getPi() )            
        );
}

//--------------------------------------------------------------------------
// PROTECTED
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// EVALUATION
//--------------------------------------------------------------------------

//  Forward
void MultiGHMM::forwardMatrix( mat sequences, unsigned int length , scaledMatrixM & forward)
{  
    unsigned int i, j, k;   

    //  Base case and forward matrix initialization
    i = 0;    
    
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0) = dmvnormSingle(sequences.col(0) , m_mu.col(i), m_sigma.slice(i), false) * m_Pi[i];                                    
        forward.scaling[0] += forward.matrix(i,0);

        
    }  
    //  First column normalization (sum = 1)
    for( i = 0; i < m_N; i++)
        forward.matrix(i,0) /=  forward.scaling[0];
    
    //  Recursive step
    //  It starts at the second observation
    for(j = 1; j < length; j++)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
                forward.matrix(i,j) += m_A(k,i) * forward.matrix(k,j-1);
            forward.matrix(i,j) *= dmvnormSingle(sequences.col(j) , m_mu.col(i), m_sigma.slice(i), false);            
            forward.scaling[j] += forward.matrix(i,j);        
        }

        //  Factor normalization
        for(i = 0; i < m_N; i++ )
            forward.matrix(i,j) /=  forward.scaling[j];
    }  

}

//  Backward
void MultiGHMM::backwardMatrix( mat sequences, unsigned int length , scaledMatrixM & backward)
{  
    unsigned int i, j, k;

    //  Base case and backward matrix initialization
    for( i = 0; i < m_N; i++)          
        backward.matrix(i, length-1)= 1;        

    //  Recursive step 
    for(j = length - 1 ; j > 0  ; j--)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
                backward.matrix(i,j-1) +=  dmvnormSingle(sequences.col(j) , m_mu.col(k), m_sigma.slice(k), false)*m_A(i,k) * backward.matrix(k,j);            
            backward.scaling[j] += backward.matrix(i,j-1);      
        }

        //  Factor normalization             
        for(i = 0; i < m_N; i++ )
            backward.matrix(i,j-1) /=  backward.scaling[j];
    } 

    //  Last step
    for(i = 0; i < m_N; i++ )
        backward.scaling[0] += m_Pi[i] * dmvnormSingle(sequences.col(0) , m_mu.col(i), m_sigma.slice(i), false) * backward.matrix(i,0); 
}

//--------------------------------------------------------------------------
// DECODING
//--------------------------------------------------------------------------

//  Function dedicated to memory allocation for the forward backward algorithm
mat MultiGHMM::forwardBackwardGamma(mat sequences)
{           
    unsigned int length = sequences.n_cols;

    //  scaling factors for the forward and backward matrices
    rowvec scaledf(length, fill::zeros);
    rowvec scaledb(length + 1, fill::zeros);  //length+1 given the prior used at the beginning of the algorithm
    scaledb[length] = 0; //log(1) = 0

    //  Memory reserved for each matrix.s
    mat matrix(m_N, length, fill::zeros);
    scaledMatrixM forward = {scaledf, matrix};
    scaledMatrixM backward = {scaledb,  matrix};    

    //  Algorithm call
    forwardBackwardGamma(sequences, forward, backward, scaledf, scaledb, matrix, length);

    //  Gamma matrix
    return matrix;
}

void MultiGHMM::forwardBackwardGamma(mat sequences, scaledMatrixM & forward, scaledMatrixM & backward,  rowvec & scaledf, rowvec & scaledb, mat & matrix, unsigned int length)
{
    unsigned int i, j, k;
    double eval;              

    //  Initial step   
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0) = dmvnormSingle(sequences.col(0) , m_mu.col(i), m_sigma.slice(i), false) * m_Pi[i];
        forward.scaling[0] += forward.matrix(i,0);
        backward.matrix(i, length-1)= 1;   
    } 
        
    //  First factor normalization (sum = 1)
    for( i = 0; i < m_N; i++)
        forward.matrix(i,0) /=  forward.scaling[0];

    //  Recursive step   
    for(j = 1; j < length; j++)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
            {
                forward.matrix(i,j) += m_A(k,i)*forward.matrix(k,j-1);
                backward.matrix(i,length-j-1) += dmvnormSingle(sequences.col(length -j ) , m_mu.col(k), m_sigma.slice(k), false) * m_A(i,k) * backward.matrix(k,length-j);
            }                
            forward.matrix(i,j) *=  dmvnormSingle(sequences.col(j) , m_mu.col(i), m_sigma.slice(i), false);                        
            
            //  Scaling factor calculation
            forward.scaling[j] += forward.matrix(i,j);
            backward.scaling[length - j] += backward.matrix(i,length-j-1);        
        }

        //  Factor normalization
        for(i = 0; i < m_N; i++ )
        {
            forward.matrix(i,j) /=  forward.scaling[j];
            backward.matrix(i,length-j-1) /=  backward.scaling[length - j];
        }
    } 
    
    //  Last step
    for(i = 0; i < m_N; i++ )
        backward.scaling[0] += m_Pi[i] * dmvnormSingle(sequences.col(0) , m_mu.col(i), m_sigma.slice(i), false) * backward.matrix(i,0);

    //  After the matrices calculation it is needed to get P(X(t) , data(1,2,...,t)) and P(data(t+1, t+2,...,T) | X(t))  
    //  To do it, a sum of logarithms was used
    scaledf[0] = log(forward.scaling[0]);
    scaledb[length-1] = log(backward.scaling[length-1]);

    for(i = 1; i < length; i++ )
    {
        scaledf[i] = scaledf[i-1] + log(forward.scaling[i]);;
        scaledb[length - 1 - i] = scaledb[length - i] + log(backward.scaling[length - 1 - i]);               
    }

    //  We get the value P(data)
    eval = scaledf[length -1];      
    
    //  We get the value P(X(t)|data) = P(X(t),data) / P(data)
    //  P(X(t)|data)  = log(P(X(t) , data(1,2,...,t))) + log(P(data(t+1, t+2,...,T) | X(t))) - log(P(data))
    for(j = 0; j < length; j++)
        for(i = 0; i < m_N; i++)
            matrix(i,j) = exp( // exponential needed to return to a probability value [0,1]
                                log(forward.matrix(i,j)) +
                                scaledf[j] +
                                log(backward.matrix(i,j)) +
                                scaledb[j+1]  // The +1 shift is done because the saled factors in the backward matrix are shifted
                                - eval );
}


//--------------------------------------------------------------------------
// LEARNING
//--------------------------------------------------------------------------

//  Function used for parameter estimation given a set of sequences
bool MultiGHMM::expectationMaximization( cube sequences, unsigned int pseudo)
{
    unsigned int seqLen, length, i, j, k, s;   
    double temp;

    seqLen = sequences.n_slices;
    length = sequences.n_cols;      

    //  Memory allocation
    mat Amean(m_N, m_N, fill::zeros);
    mat muMean(m_M, m_N, fill::zeros);
    cube sigmaMean(m_M, m_M, m_N, fill::zeros);
    rowvec Pimean(m_N, fill::zeros);

    //  Normalizing factors
    rowvec denomA(m_N, fill::zeros);
    rowvec denomB(m_N, fill::zeros);                   
         
    //  Analysis per sequence  
    for(s = 0; s < seqLen ; s++)
    {        
        //  Expectation step
        //  Memory allocation for expectation step
        rowvec scaledf(length, fill::zeros);
        rowvec scaledb(length + 1, fill::zeros);

        mat matrix(m_N, length, fill::zeros); //  Gamma matrix         
        scaledMatrixM forward = {scaledf, matrix};
        scaledMatrixM backward = {scaledb,  matrix};

        mat observation = sequences.slice(s);
        
        forwardBackwardGamma(observation, forward, backward, scaledf, scaledb, matrix, length);

        //  Maximization         
        //  Hidden state X(t)
        for( i = 0; i < m_N; i++)
        {
            Pimean[i] += matrix(i,0);                                                           

            //  Each observation in the sequence 's' 
            for(j = 0; j < (length-1); j++)
            {                
                //  Hidden state X(t+1)
                for(k = 0; k < m_N ; k++ )
                {
                    temp = (matrix(i,j)*m_A(i,k)*dmvnormSingle(observation.col(j+1) , m_mu.col(k), m_sigma.slice(k), false)*backward.matrix(k,j+1))/
                        (backward.matrix(i,j)* backward.scaling[j+1] ); // j+1 is the shift in the backward scaling vector 
                    Amean(i,k) += temp;                    
                    denomA[i] += temp;
                }
                
                //  For the emission matrix, we can use P(X|data)
                muMean.col(i) = muMean.col(i) + (matrix(i,j)*observation.col(j));                   //  Means
                colvec varTemp = (observation.col(j) - m_mu.col(i));                                //  Difference between the mean and the observation  
                sigmaMean.slice(i) = sigmaMean.slice(i) + (matrix(i,j) * (varTemp * varTemp.t()));  //  Variance
                denomB[i] += matrix(i,j);
            }
            //  The loop ends at 'length-1', thus it is necessary to sum the las value             
            muMean.col(i) = muMean.col(i) + (matrix(i,j)*observation.col(length-1));                //  Means
            colvec varTemp = (observation.col(length-1) - m_mu.col(i));                             //  Difference between the mean and the observation  
            sigmaMean.slice(i) = sigmaMean.slice(i) + (matrix(i,j) * (varTemp * varTemp.t()));      //  Variance
            denomB[i] += matrix(i,length-1); 
        }                        
    } 
           
    //  We use the normalizing factor and also use the pseudo count value
    for(i = 0; i < m_N; i++)
    {
        //  Initial vector normalization
        m_Pi[i] = (Pimean[i] + pseudo) / (seqLen + m_N*pseudo);
        //  Transition matrix normalization
        for(k = 0; k < m_N; k++)
            m_A(i,k) = (Amean(i,k) + pseudo) / (denomA[i] + m_N*pseudo);             
        //  Emission matrix normalization
        m_mu.col(i) = muMean.col(i) / denomB[i];                          //  Normalized means
        if( isPositiveDefinite(sigmaMean.slice(i)) )
            //m_sigma.slice(i) = sqrtmat_sympd(sigmaMean.slice(i) / denomB[i]); //  standard deviation       
            m_sigma.slice(i) = sigmaMean.slice(i) / denomB[i]; //  Covariance matrix
        else{
            return false;
        }
    }   
    return true;
}

//  Function used for parameter estimation given a single sequence.
//  The initial probability vector is not estimated
bool MultiGHMM::BaumWelch( mat sequences, unsigned int pseudo)
{
    unsigned int length, i, j, k;   
    double temp;
    length = sequences.n_cols;      

    //  Memory allocation
    mat Amean(m_N, m_N, fill::zeros);
    mat muMean(m_M, m_N, fill::zeros);
    cube sigmaMean(m_M, m_M, m_N, fill::zeros);                   

    //  Normalizing factors
    rowvec denomA(m_N, fill::zeros);
    rowvec denomB(m_N, fill::zeros);                      
         
    //  Expectation step
    //  Memory allocation for expectation step
    rowvec scaledf(length, fill::zeros);
    rowvec scaledb(length+ 1, fill::zeros);

    mat matrix(m_N, length, fill::zeros); //  Gamma matrix        
    scaledMatrixM forward = {scaledf, matrix};
    scaledMatrixM backward = {scaledb,  matrix};
        
    forwardBackwardGamma(sequences, forward, backward, scaledf, scaledb, matrix, length);
    
    //  Maximization         
    //  Hidden state X(t)
    for( i = 0; i < m_N; i++)
    {                                                             
        //  Each observation in the sequence
        for(j = 0; j < (length-1); j++)
        {                
            //  Hidden state X(t+1)
            for(k = 0; k < m_N ; k++)
            {                
                temp = (matrix(i,j) * m_A(i,k) * dmvnormSingle(sequences.col(j+1) , m_mu.col(k), m_sigma.slice(k), false) * backward.matrix(k,j+1))/
                        (backward.matrix(i,j) * backward.scaling[j+1] ); // j+1 is the shift in the backward scaling vector
                Amean(i,k) += temp;                
                denomA[i] += temp;
            }
                
            //  Emission matrix estimation
            muMean.col(i) = muMean.col(i) + (matrix(i,j)*sequences.col(j));                     //  Means
            colvec varTemp = (sequences.col(j) - m_mu.col(i));                                  //  Difference between the mean and the observation  
            sigmaMean.slice(i) = sigmaMean.slice(i) + (matrix(i,j) * (varTemp * varTemp.t()));  //  Variance
            denomB[i] += matrix(i,j);
        }
        //  The loop ends at 'length-1', thus it is necessary to sum the las value             
        muMean.col(i) = muMean.col(i) + (matrix(i,j)*sequences.col(length-1));                //  Means
        colvec varTemp = (sequences.col(length-1) - m_mu.col(i));                             //  Difference between the mean and the observation  
        sigmaMean.slice(i) = sigmaMean.slice(i) + (matrix(i,j) * (varTemp * varTemp.t()));      //  Variance
        denomB[i] += matrix(i,length-1); 
    }                            
           
    //  We use the normalizing factor and also use the pseudo count value
    for(i = 0; i < m_N; i++)
    {
        //  Transition matrix normalization
        for(k = 0; k < m_N; k++)
            m_A(i,k) = (Amean(i,k) + pseudo) / (denomA[i] + m_N*pseudo);             

        //  Emission matrix normalization
        m_mu.col(i) = muMean.col(i) / denomB[i];                          //  Normalized means
        if( isPositiveDefinite(sigmaMean.slice(i)) )
            //m_sigma.slice(i) = sqrtmat_sympd(sigmaMean.slice(i) / denomB[i]); //  standard deviation    
            m_sigma.slice(i) = sigmaMean.slice(i) / denomB[i]; //  Covariance matrix    
        else{
            return false;
        }
    }
    return true;
}

//--------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------

//  Function used to set random model parameters
void MultiGHMM::randomInit(double min, double max)
{

    //  Used for set.seed compatibility
    RNGScope scope;

    // Sigma variables
    mat temp;

    // Random values:
    m_A.randu();
    m_Pi.randu();
    m_mu.randu();    

    // We set the real range to the means    
    m_mu = m_mu*(max - min)  + min;
    
    //  Normalizing factors
    double maxPi = accu(m_Pi);  // Sum all values
    colvec maxA = sum(m_A,1);   // Sum over the each row

    m_Pi = m_Pi / maxPi;     
    for(int i = 0; i < m_N; i++)
    {
        m_A.row(i) = m_A.row(i)/maxA(i);
        // Each slice in sigma must be positive definite
        // Since a symmetric, positive definite matrix is uniquely determined by its Cholesky decomposition, 
        // you could just randomly choose a lower triangular matrix LL with positive diagonal entries and 
        // obtain your matrix as L* t(L).
        do{
            temp = (randu<mat>(m_M, m_M) * (max - min) ); //Range := [0, max - min ]      
            temp = trimatl(temp);
            temp = temp * temp.t();
            temp.diag() += 1 ;//+ randu<vec>(m_M); // To make sure that the variance of each dimension is positive                      
        }while(isPositiveDefinite(temp) == false);
        m_sigma.slice(i) = temp;          
    }
}

//  Verify if it is a stochastic matrix (row elements sum equals to 1)
bool MultiGHMM::verifyMatrix(mat matrix)
{
    bool flag = true;  
    double normalized;
    unsigned int i;
    
    //  Sum of each row must be 1
    for(i = 0 ; i < matrix.n_rows ; i++)
    {
        normalized = sum(matrix.row(i)); 
        if(normalized < 1 - EPSILON || normalized > 1 + EPSILON)
        {
            flag = false;
            break;
        } 
    }
    return flag ;    
}

//  Verify if the sum of all its elements sum 1
bool MultiGHMM::verifyVector(rowvec vector)
{
    bool flag = true;       

    double normalized = sum(vector); 
    if(normalized < 1 - EPSILON || normalized > 1 + EPSILON)
            flag = false;
 
    return flag;
}

//  All the data expressed as discrete {0,1,2,...} is returned to categorical
//  S: State, O: Observation
Rcpp::CharacterVector MultiGHMM::toName( Rcpp::IntegerVector index, char vectorName)
{
    int length = 0, i = 0;

    length = index.size(); 
    CharacterVector names(length);       

    switch(vectorName)
    {
        case 'S':
            for(i = 0; i < length; i++)
                names[i] = m_StateNames[index[i]];
        break;        
    }

    return names;
}



