// included dependencies

#include "HMM.h"

//System libraries
#include <iostream>
#include <string>
#include <algorithm>    /* find */
#include <math.h>       /* log */
//#include <omp.h>      /* multithead */

//namespaces
using namespace std;
using namespace Rcpp;

//--------------------------------------------------------------------------
// CONSTRUCTORS & DESTRUCTOR:
//--------------------------------------------------------------------------

//  First constructor
HMM::HMM(unsigned short int  numberStates, unsigned short int  numberEmissions)
{
    //  Validate the values
    if(numberStates < 2 || numberEmissions < 2)    
        Rf_error("The number of states and possible symbols must be bigger or equal to 2.");

    //  Set known values 
    m_N = numberStates;
    m_M = numberEmissions;
    m_StateNames = CharacterVector(m_N);
    m_ObservationNames = CharacterVector(m_M);

    // Memory allocation for parameters
    m_A = NumericMatrix(m_N,m_N);
    m_B = NumericMatrix(m_N,m_M);
    m_Pi = NumericVector(m_N);
    /**********************************************************/
    //  Proposed state and emission names
    for(int i = 1; i <= m_N; i++ )
        m_StateNames[i-1] = "x" + to_string(i); 

    for(int i = 1; i <= m_M; i++ )
        m_ObservationNames[i-1] = "y" + to_string(i);           
    
    //  Parameter random initialization
    randomInit();
}

//  Second constructor
HMM::HMM(CharacterVector stateNames, CharacterVector emissionNames)
{
    //Validate the values
    if(stateNames.size() < 2 || emissionNames.size() < 2)    
        Rf_error("The number of states and possible symbols must be bigger or equal to 2.");    

    //  Set known values 
    m_N = stateNames.size();
    m_M =  emissionNames.size();
    m_StateNames = stateNames ;
    m_ObservationNames = emissionNames ;

    // Memory allocation for parameters
    m_A = NumericMatrix(m_N,m_N);
    m_B = NumericMatrix(m_N,m_M);
    m_Pi = NumericVector(m_N); 

    //  Parameter random initialization
    randomInit();
}

//  Third constructor used for model validation 
HMM::HMM(CharacterVector stateNames, CharacterVector emissionNames, NumericMatrix A, NumericMatrix B, NumericVector Pi)
{
    //Validate the values
    if(stateNames.size() < 2 || emissionNames.size() < 2)    
        Rf_error("The number of states and possible symbols must be bigger or equal to 2.");
    
    //Validate the values
    if(stateNames.size() != A.ncol() || stateNames.size() != A.nrow())    
        Rf_error("The number of states must be the same as the transition matrix column and row size");

    //Validate the values
    if(emissionNames.size() != B.ncol() || stateNames.size() != B.nrow())    
        Rf_error("The number of symbols must be the same as the emission matrix column size and the number of states must be the same as the row size");
    
    //Validate the values
    if(stateNames.size() != Pi.size())    
        Rf_error("The number of states must be the same as the initial probability vector size");

    //  If all the paremeters have been validated, then they are used.
    m_N = stateNames.size();
    m_M =  emissionNames.size();    
    m_StateNames = stateNames ;
    m_ObservationNames = emissionNames ;

    setParameters(A,B,Pi);
}

//  Destructor
HMM::~HMM(void){}

//--------------------------------------------------------------------------
//  PUBLIC
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//  GETTERS
//--------------------------------------------------------------------------

CharacterVector HMM::getEmissionNames(void) const
{

    return m_ObservationNames;
}

NumericMatrix HMM::getB(void) const
{
    return m_B;
}

//--------------------------------------------------------------------------
//  SETTERS
//--------------------------------------------------------------------------

void HMM::setEmissionNames(CharacterVector emissionNames)
{
    if( emissionNames.size() != m_M)
        Rf_error("The number of state names does not coincide with the one declared.");
    m_ObservationNames = CharacterVector(clone(emissionNames));
}

void HMM::setB(NumericMatrix B)
{
    if(B.ncol() != m_M || B.nrow() != m_N)
        Rf_error("The emission matrix size is wrong");
    if(verifyMatrix(B) == false)
        Rf_error("The emission matrix is not normalized");
    m_B = NumericMatrix(clone(B));
}

void HMM::setParameters(NumericMatrix A, NumericMatrix B, NumericVector Pi)
{
    if(Pi.size() != m_N)
        Rf_error("The initial probability vector size is wrong");
    if(verifyVector(Pi) == false)
        Rf_error("The initial probability vector is not normalized");
    
    if(A.ncol() != m_N || A.nrow() != m_N)
        Rf_error("The transition matrix size is wrong");
    if(verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized");  

    setB(B);            
    m_Pi = NumericVector(clone(Pi));
    m_A = NumericMatrix(clone(A));    
}

//--------------------------------------------------------------------------
//  EVALUATION
//--------------------------------------------------------------------------

double HMM::evaluation(CharacterVector sequence, char method)
{    
    unsigned int i, length;
    double eval = 0.0;
    
    //  Change categorical data into discrete values
    length = sequence.size();
    IntegerVector index = toIndex(sequence);    

    // Memory allocation
    NumericVector scaled(length, 0);
    NumericMatrix matrix(m_N, length);
    scaledMatrix eva = {scaled, matrix};

    //  Selected evaluation method
    switch(method)
    {  
        case 'b':
            backwardMatrix( index, length , eva );                
            break;
        default:            
            forwardMatrix( index, length , eva ); 
            break;
    }

    //  double variables lose precision if multiplication is used
    //  Therfore, a sum of log values is used
    for(i = 0; i < length; i++ )
                eval+= log(eva.scaling[i]);
    
    return eval;
}

//--------------------------------------------------------------------------
// DECODING
//--------------------------------------------------------------------------

CharacterVector HMM::forwardBackward(CharacterVector sequence)
{
    unsigned j;
    unsigned int length = sequence.size();

    //  We get P(State|Data) for all states and observations
    NumericMatrix gamma = forwardBackwardGamma(sequence);
    
    //  Most probable state-path traveled
    IntegerVector best(length,0);

    //  Temp vector to store a column 
    NumericVector temp(m_N, 0); 
    for(j = 0; j < length; j++)
    {
        temp = gamma(_,j);
        best[j] = distance(temp.begin(), max_element(temp.begin(), temp.end()));
    }

    //  Change discrete data into categorical values
    return toName(best, 'S');  
}

CharacterVector HMM::viterbi(CharacterVector sequence)
{
    unsigned int length, i, j, k;    
    
    length = sequence.size();      

    //  Most probable state-path traveled 
    IntegerVector best(length, 0);

    // Memory allocation
    NumericMatrix phi(m_N, length);
    NumericMatrix delta(m_N, length); 
    
    //  Change categorical data into discrete values
    IntegerVector index = toIndex(sequence); 

    //  log-transform to avoid double precision loss
    NumericMatrix A(m_N, m_N);
    NumericMatrix B(m_N, m_M);
    NumericVector Pi = log(m_Pi);  
    
    for(i = 0 ; i < m_N; i++)
    {                
        A(i,_)= log(m_A(i,_)); 
        B(i,_)= log(m_B(i,_));
    }

    // Used to find the max value
    NumericVector temp(m_N); 

    //  Init step
    for(i = 0; i < m_N; i++)
        delta(i,0)= Pi[i] +  B(i,index[0]);

    //  Recursion step
    for(j = 1; j < length; j++)
        for(i = 0; i < m_N; i++)
        {
            for(k = 0; k < m_N; k++)
                temp[k] = delta(k,j-1) + A(k,i);
            // The auto keyword is simply asking the compiler to deduce the type of the variable from the initialization                
            auto maximum = max_element(temp.begin(), temp.end());   
            //  (*maximum) to access the value of the pointer
            delta(i,j)= (*maximum) + B(i,index[j]);
            phi(i,j) = distance(temp.begin(), maximum);             
        }

    // Termination step
    temp = delta(_,length-1);            

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
double HMM::loglikelihood(CharacterMatrix sequences)
{
    double ll = 0.0;
    unsigned int i, seqLen;

    seqLen = sequences.nrow();

    for(i = 0; i < seqLen; i++)
        ll+=evaluation(sequences.row(i));
    return ll;
}

//  Parameter estimation using a Expectation Maximization approach for multiple sequences
void HMM::learnEM(CharacterMatrix sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{
    double newLL, error;
    double lastLL = loglikelihood(sequences);
    unsigned int counter = 0;       
    
    //  Parameter estimation
    do{
        //  If the error is nan, it may be a big error. 
        //  A new parameter initialization is recomended
        expectationMaximization(sequences, pseudo);
        newLL = loglikelihood(sequences) ;
        if(std::isnan(newLL))
        {
            if(print)
                Rcout << "Convergence error, new initialization needed\n";                       
            randomInit();
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
    } while(counter < iter && error > delta); // Convergence criteria
    
    Rcout << "Finished at Iteration: " << counter << " with Error: " << error  << "\n";
}

//  Parameter estimation using a Baum Welch approach for a single sequence
void HMM::learnBW(CharacterVector sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{
    double newLL, error;
    double lastLL = evaluation(sequences);
    unsigned int counter = 0;       
    
    //  Parameter estimation
    do{
        //  If the error is nan, it may be a big error. 
        //  A new parameter initialization is recomended
        BaumWelch(sequences, pseudo);
        newLL = evaluation(sequences);
        if(std::isnan(newLL))
        {
            if(print)
                Rcout << "Convergence error, new initialization needed\n";            
            randomInit();
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
    } while(counter < iter && error > delta);
    Rcout << "Finished at Iteration: " << counter << " with Error: " << error  << "\n";
}

//--------------------------------------------------------------------------
// SIMULATION
//--------------------------------------------------------------------------

//  Funtion used to generate observations given a model
List HMM::generateObservations( unsigned short int length)
{
    unsigned int i,j;
    double x, y;    
    IntegerVector X(length, 0);
    IntegerVector Y(length, 0);

    //  Used for set.seed compatibility
    RNGScope scope;

    //  Matrix rearrangement to use a uniform distribution to generate the hidden states and its corresponding observation
    NumericMatrix A(m_N, m_N);
    NumericMatrix B(m_N, m_M);
    NumericVector Pi(m_N);

    double tempPi = 0.0, tempAB = 0.0;

    //  We  fill each value with its corresponding new one
    for(i=0; i < m_N; i++)
    {
        //  We fill first the initial probability vector
        tempPi += m_Pi[i];
        Pi[i] = tempPi;        
        
        //  Then, we fill the transition matrix
        tempAB = 0;
        for(j = 0; j < m_N; j++)
        {
            tempAB += m_A(i,j);
            A(i,j)= tempAB;
        }
        
        //  Finally, we fill the emission matrix
        tempAB = 0;
        for(j = 0; j < m_M; j++)
        {
            tempAB += m_B(i,j);
            B(i,j)= tempAB;                     
        }
    }

    //  Random variable feneration based in the rearranged matrices
    x = as<double>(runif(1));
    y = as<double>(runif(1));

    NumericVector tempA;
    NumericVector tempB = B.row(X[0]);
    
    X[0] = lower_bound (Pi.begin(), Pi.end(), x) - Pi.begin();
    Y[0] = lower_bound (tempB.begin(), tempB.end(), y) - tempB.begin();  

    for(j = 1; j < length; j++)
    {
        x = as<double>(runif(1));
        y = as<double>(runif(1));

        tempA = A.row(X[j-1]);
        X[j] = lower_bound (tempA.begin(), tempA.end(), x) - tempA.begin();
        tempB = B.row(X[j]);
        Y[j] = lower_bound (tempB.begin(), tempB.end(), y) - tempB.begin();
    }
    //  Returns the hidden state path 'X' and its emissions 'Y'
    return List::create(
                Named("X", toName(X, 'S')),
                Named("Y", toName(Y, 'O'))
            );
}

//--------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------

//  Funtion to return all the model parameters as an R List
List HMM::toList(void) const
{
    return List::create(
            Named("Model", "HMM"),
            Named("StateNames", getStateNames() ),
            Named("ObservationNames", getEmissionNames() ),
            Named("A", getA()),
            Named("B", getB()),
            Named("Pi", getPi() )
    );
}

//--------------------------------------------------------------------------
// PROTECTED
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// EVALUATION
//--------------------------------------------------------------------------

//  Forward
void HMM::forwardMatrix( IntegerVector sequence, unsigned int length , scaledMatrix & forward)
{  
    unsigned int i, j, k;   

    //  Base case and forward matrix initialization
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0) = m_B(i,sequence[0])* m_Pi[i];
        //  A scaling factor is needed to avoid a double variable precision loss
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
                forward.matrix(i,j) += m_A(k,i)*forward.matrix(k,j-1);
            forward.matrix(i,j) *= m_B(i,sequence[j]);
            //  Scaling factor calculation
            forward.scaling[j] += forward.matrix(i,j);                  
        }

        //  Factor normalization
        for(i = 0; i < m_N; i++ )
            forward.matrix(i,j) /=  forward.scaling[j];            
    }   
}

//  Backward
void HMM::backwardMatrix( IntegerVector sequence, unsigned int length , scaledMatrix & backward)
{  
    unsigned int i, j, k;

    //  Base case and backward matrix initialization
    for( i = 0; i < m_N; i++)          
        backward.matrix(i,length - 1) = 1;    

    //  Recursive step    
    for(j = length - 1 ; j > 0  ; j--)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
                backward.matrix(i,j-1) += m_B(k,sequence[j]) * m_A(i,k) * backward.matrix(k,j);  
            //  Scaling factor calculation                        
            backward.scaling[j] += backward.matrix(i,j-1);      
        }

        //  Factor normalization        
        for(i = 0; i < m_N; i++ )
            backward.matrix(i,j-1) /=  backward.scaling[j];
    } 

    //  Last step
    for(i = 0; i < m_N; i++ )
        backward.scaling[0] += m_Pi[i] * m_B(i,sequence[0])* backward.matrix(i,0); 
}

//--------------------------------------------------------------------------
// DECODING
//--------------------------------------------------------------------------

//  Function dedicated to memory allocation for the forward backward algorithm
NumericMatrix HMM::forwardBackwardGamma(CharacterVector sequence)
{           
    unsigned int length = sequence.size();     

    //  scaling factors for the forward and backward matrices
    NumericVector scaledf(length, 0);
    NumericVector scaledb(length + 1, 0);  //length+1 given the prior used at the beginning of the algorithm
    scaledb[length] = 0; //log(1) = 0

    //  Memory reserved for each matrix. The matrices are cloned to make a "safe" memory management
    NumericMatrix matrix(m_N, length);
    scaledMatrix forward = {clone(scaledf), clone(matrix)};
    scaledMatrix backward = {clone(scaledb),  clone(matrix)};

    // Change categorical data to discrete values    
    IntegerVector index = toIndex(sequence); 

    //  Algorithm call
    forwardBackwardGamma(index, forward, backward, scaledf, scaledb, matrix, length);

    //  Gamma matrix
    return matrix;
}

void HMM::forwardBackwardGamma(IntegerVector index, scaledMatrix & forward, scaledMatrix & backward,  NumericVector & scaledf, NumericVector & scaledb, NumericMatrix & matrix, unsigned int length)
{
    unsigned int i, j, k;
    double eval;              
    
    //  Initial step    
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0) = m_B(i,index[0]) * m_Pi[i];
        forward.scaling[0] += forward.matrix(i,0);
        backward.matrix(i,length - 1) = 1;   
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
                backward.matrix(i,length-j-1) += 
                    m_B(k, index[length-j]) * m_A(i,k) * backward.matrix(k,length - j);
            }                
            forward.matrix(i,j) *= m_B(i,index[j]);
            
            //  Scaling factor calculation              
            forward.scaling[j] += forward.matrix(i,j);
            backward.scaling[length - j] += backward.matrix(i,length-j-1);        
        }

        //  Factor normalization
        for(i = 0; i < m_N; i++ )
        {
            forward.matrix(i,j) /=  forward.scaling[j];
            backward.matrix(i,length - j - 1) /=  backward.scaling[length-j];   
        }
    } 
    
    //  Last step
    for(i = 0; i < m_N; i++ )
        backward.scaling[0] += m_Pi[i] * m_B(i,index[0]) * backward.matrix(i,0);

    //  After the matrices calculation it is needed to get P(X(t) , data(1,2,...,t)) and P(data(t+1, t+2,...,T) | X(t))  
    //  To do it, a sum of logarithms was used
    scaledf[0] = log(forward.scaling[0]);
    scaledb[length-1] = log(backward.scaling[length-1]);

    for(j = 1; j < length; j++ )
    {
        scaledf[j] = scaledf[j-1] + log(forward.scaling[j]);;
        scaledb[length - 1 - j] = scaledb[length - j] + log(backward.scaling[length - 1 - j]);               
    }
    
    //  We get the value P(data)
    eval = scaledf[length - 1];      

    //  We get the value P(X(t)|data) = P(X(t),data) / P(data)
    //  P(X(t)|data)  = log(P(X(t) , data(1,2,...,t))) + log(P(data(t+1, t+2,...,T) | X(t))) - log(P(data))
    for(j = 0; j < length; j++)
        for(i = 0; i < m_N; i++)          
            matrix(i,j) = exp( // exponential needed to return to a probability value [0,1]
                                log(forward.matrix(i,j)) +
                                scaledf[j] +
                                log(backward.matrix(i,j)) +
                                scaledb[j+1]  - // The +1 shift is done because the saled factors in the backward matrix are shifted
                                eval );
     
}

//--------------------------------------------------------------------------
// LEARNING
//--------------------------------------------------------------------------

//  Function used for parameter estimation given a set of sequences
void HMM::expectationMaximization( CharacterMatrix sequences, unsigned int pseudo)
{
    unsigned int seqLen, length, i, j, k, s;
    IntegerVector index;

    double temp;    

    seqLen = sequences.nrow();
    length = sequences.ncol();    

    //  Memory allocation
    NumericMatrix Amean(m_N, m_N);
    NumericMatrix Bmean(m_N, m_M);            
    NumericVector Pimean(m_N);

    //  Normalizing factors
    NumericVector denomA(m_N);
    NumericVector denomB(m_N);                      
         
    //  Analysis per sequence  
    for(s = 0; s < seqLen ; s++)
    {       
        //  Change categorical data into discrete values
        index = toIndex(sequences.row(s)); 
        //  Expectation step
        //  Memory allocation for expectation step
        NumericVector scaledf(length);
        NumericVector scaledb(length + 1);

        NumericMatrix matrix(m_N, length); //  Gamma matrix       
        scaledMatrix forward = {clone(scaledf), clone(matrix)};
        scaledMatrix backward = {clone(scaledb),  clone(matrix)};
        
        forwardBackwardGamma(index, forward, backward, scaledf, scaledb, matrix, length);

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
                    temp = (matrix(i,j) * m_A(i,k) * m_B(k,index[j+1]) * backward.matrix(k,j+1))/
                                        (backward.matrix(i,j) * backward.scaling[j+1] ); // j+1 is the shift in the backward scaling vector
                    Amean(i,k) += temp;                    
                    denomA[i] += temp;                                       
                }
                
                //  For the emission matrix, we can use P(X|data)
                Bmean(i,index[j]) += matrix(i,j);
                denomB[i] += matrix(i,j);
            }
            //  The loop ends at 'length-1', thus it is necessary to sum the las value
            Bmean(i,index[length-1]) += matrix(i,length-1);
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
        for(j = 0; j < m_M; j++)
            m_B(i,j) = (Bmean(i,j) + pseudo) / (denomB[i] + m_M*pseudo );
    }
}


//  Function used for parameter estimation given a single sequence.
//  The initial probability vector is not estimated
void HMM::BaumWelch( CharacterVector sequence, unsigned int pseudo)
{
    unsigned int length, i, j, k;
    IntegerVector index;

    double temp;      
    length = sequence.size();    

    //  Memory allocation
    NumericMatrix Amean(m_N, m_N);
    NumericMatrix Bmean(m_N, m_M);               

    //  Normalizing factors
    NumericVector denomA(m_N);
    NumericVector denomB(m_N);                      
         
    //  Change categorical data into discrete values
    index = toIndex(sequence); 
    //  Expectation step
    //  Memory allocation for expectation step
    NumericVector scaledf(length);
    NumericVector scaledb(length + 1);

    NumericMatrix matrix(m_N, length); //  Gamma matrix             
    scaledMatrix forward = {clone(scaledf), clone(matrix)};
    scaledMatrix backward = {clone(scaledb),  clone(matrix)};
        
    forwardBackwardGamma(index, forward, backward, scaledf, scaledb, matrix, length);

    //  Maximization         
    //  Hidden state X(t)
    for( i = 0; i < m_N; i++)
    {                                                              
        //  Each observation in the sequence
        for(j = 0; j < (length-1); j++)
        {                
            //  Hidden state X(t+1)
            for(k = 0; k < m_N ; k++ )
            {
                temp = (matrix(i,j) * m_A(i,k) * m_B(k,index[j+1]) * backward.matrix(k,j+1))/
                    (backward.matrix(i,j) * backward.scaling[j+1] ); // j+1 is the shift in the backward scaling vector
                Amean(i,k) += temp;                
                denomA[i] += temp;                                       
            }              
            //  For the emission matrix, we can use P(X|data)
            Bmean(i,index[j]) += matrix(i,j);
            denomB[i] += matrix(i,j);
        }
        //  The loop ends at 'length-1', thus it is necessary to sum the las value
        Bmean(i,index[length-1]) += matrix(i,length-1);
        denomB[i] += matrix(i,length-1);
    }                   

    //  We use the normalizing factor and also use the pseudo count value
    for(i = 0; i < m_N; i++)
    {        
        //  Transition matrix normalization
        for(k = 0; k < m_N; k++)
            m_A(i,k) = (Amean(i,k) + pseudo) / (denomA[i] + m_N*pseudo);             

        //  Emission matrix normalization
        for(j = 0; j < m_M; j++)
            m_B(i,j) = (Bmean(i,j) + pseudo) / (denomB[i] + m_M*pseudo );
    }
}

//--------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------

//  Function used to set random model parameters
void HMM::randomInit()
{    
    //  Used for set.seed compatibility
    RNGScope scope;

    //  Normalizing factors
    double maxPi = 0.0;
    NumericVector maxA(m_N, 0.0);
    NumericVector maxB(m_N, 0.0);
 

    for (int i=0; i< m_N ; i++) 
    {
        //  Initial probability vector
        m_Pi[i] = as<double>(runif(1));
        maxPi += m_Pi[i] ;
        maxA[i] = 0;
        maxB[i] = 0;
        //  Transition matrix 
        for(int j = 0; j < m_N; j++)
        {
            m_A(i,j) = as<double>(runif(1));
            maxA[i] += m_A(i,j);
        }
        //  Emission matrix
        for(int j = 0; j < m_M; j++)
        {
            m_B(i,j) = as<double>(runif(1));
            maxB[i] += m_B(i,j);
        }
    }

    //  Parameter normalization
    for(int i = 0; i < m_N; i++)
    {
        m_Pi[i] /= maxPi;        
        for(int j = 0; j < m_N; j++)
            m_A(i,j) /= maxA[i];        
        for(int j = 0; j < m_M; j++)
            m_B(i,j) /= maxB[i];
    }  
}

//  Function used to change the categorial data into discrete values  
//  Function used for the provided sequences
IntegerVector HMM::toIndex( CharacterVector observations)
{
    int length = observations.size();
    IntegerVector indices(length);

    int pos;

    for(int i = 0; i < length; i++)
    {
        //  First it is necessary to confirm that all values provided by the user are in the observation names set.
        pos = find(m_ObservationNames.begin(), m_ObservationNames.end(), observations[i]) - m_ObservationNames.begin();
        if(pos < m_ObservationNames.size())
            indices[i] = pos;
        else
        {
            Rcout << "Error in "  <<  observations[i] <<  " , " << i << endl;        
            Rf_error("The values must exist in the possible observations of the model");
        }
            
    }     
    return indices;
}


//  Function used to change the discrete data into categorical values
//  S: State, O: Observation
CharacterVector HMM::toName( IntegerVector index, char vectorName)
{
    unsigned int length, i; 
    length = index.size();
    CharacterVector names(length);

    switch(vectorName)
    {
        case 'S':
            for(i = 0; i < length; i++)
                names[i] = m_StateNames[index[i]];
        break;
        case 'O':
            for(i = 0; i < length; i++)
                names[i] = m_ObservationNames[index[i]];
        break;
    }

    return names;
}

