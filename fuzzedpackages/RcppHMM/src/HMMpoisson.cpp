// included dependencies

#include "HMMpoisson.h"

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
HMMpoisson::HMMpoisson(unsigned short int  numberStates)
{
    //Validate the values
    if(numberStates < 2)
        Rf_error("The number of states must be bigger or equal to 2.");        

    //  Set known values 
    m_N = numberStates;    
    m_StateNames = CharacterVector(m_N);

    // Memory allocation for parameters
    m_A = NumericMatrix(m_N,m_N);
    m_B = NumericVector(m_N);
    m_Pi = NumericVector(m_N);

    /**********************************************************/
    //  Proposed state names
    for(int i = 1; i <= m_N; i++ )
        m_StateNames[i-1] = "x" + to_string(i); 
    
    //  Parameter random initialization                                    
    randomInit(1,10);
}

//  Second constructor
HMMpoisson::HMMpoisson(CharacterVector stateNames)
{
    //Validate the values
    if(stateNames.size() < 2)
        Rf_error("The number of states must be bigger or equal to 2.");

    //  Set known values 
    m_N = stateNames.size();    
    m_StateNames = stateNames ;    

    // Memory allocation for parameters
    m_A = NumericMatrix(m_N,m_N);
    m_B = NumericVector(m_N);
    m_Pi = NumericVector(m_N);                  
                     
    //  Parameter random initialization                     
    randomInit(1,10);
}

//  Third constructor used for model validation 
HMMpoisson::HMMpoisson(CharacterVector stateNames, NumericMatrix A, NumericVector B, NumericVector Pi)
{
    //Validate the values
    if(stateNames.size() < 2) 
        Rf_error("The number of states must be bigger or equal to 2.");           
    
    //Validate the values
    if(stateNames.size() != A.ncol() || stateNames.size() != A.nrow())            
        Rf_error("The number of states must be the same as the transition matrix column and row size");

    //Validate the values
    if(stateNames.size() != B.size())              
        Rf_error("The number of parameters in the emission matrix must be the number of states");        
    
    //Validate the values
    if(stateNames.size() != Pi.size())            
        Rf_error("The number of states must be the same as the initial probability vector size");

    //  If all the paremeters have been validated, then they are used.
    m_N = stateNames.size();
    m_StateNames = stateNames ;

    m_A = NumericMatrix(m_N,m_N);
    m_B = NumericVector(m_N);
    m_Pi = NumericVector(m_N);                         
       
    setParameters(A,B,Pi);  
}

//  Destructor
HMMpoisson::~HMMpoisson(void){}

//--------------------------------------------------------------------------
//  PUBLIC
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//  GETTERS
//--------------------------------------------------------------------------

NumericVector HMMpoisson::getB(void) const
{    
    return m_B;
}

//--------------------------------------------------------------------------
//  SETTERS
//--------------------------------------------------------------------------

void HMMpoisson::setB(NumericVector B)
{        
    if(B.size() != m_N)
        Rf_error("The emission vector size is wrong");

    for(int i = 0; i < m_N; i++)    
        if(B[i] <= 0)
            throw invalid_argument("Lambda must be greater than zero"); 
    
    m_B = NumericVector(clone(B));           
}

void HMMpoisson::setParameters(NumericMatrix A, NumericVector B, NumericVector Pi)
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

double HMMpoisson::evaluation(IntegerVector sequence, char method)
{
    //Validate the sequence
    if(sum(sequence) < sum(abs(sequence))) 
        Rf_error("All values in the sequence must be positive");

    double eval = 0.0;
    unsigned int i, length;
    
    length = sequence.size();    

    // Memory allocation
    NumericVector scaled(length);
    NumericMatrix matrix(m_N, length);
    scaledMatrix eva = {scaled, matrix};

    //  Selected evaluation method
    switch(method)
    {
        case 'f':
            forwardMatrix( sequence, length , eva );             
            break;
        case 'b':
            backwardMatrix( sequence, length , eva );                
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

CharacterVector HMMpoisson::forwardBackward(IntegerVector sequence)
{
    //Validate the sequence
    if(sum(sequence) < sum(abs(sequence))) 
        Rf_error("All values in the sequence must be positive");

    unsigned i, j;
    unsigned int length = sequence.size();

    //  We get P(State|Data) for all states and observations
    NumericMatrix gamma = forwardBackwardGamma(sequence);

    //  Most probable state-path traveled
    IntegerVector best(length);

    //  Temp vector to store a column 
    NumericVector temp(m_N); 
    for(j = 0; j < length; j++)
    {
        for(i = 0; i < m_N ; i++)
            temp[i] = gamma(i,j);
        best[j] = distance(temp.begin(), max_element(temp.begin(), temp.end()));
    }

    //  Change discrete data into categorical values
    return toName(best, 'S');  
}

CharacterVector HMMpoisson::viterbi(IntegerVector sequence)
{
    //Validate the sequence
    if(sum(sequence) < sum(abs(sequence))) 
        Rf_error("All values in the sequnce must be positive");

    unsigned int length, i, j, k;    
    
    length = sequence.size();       

    //  Most probable state-path traveled 
    IntegerVector best(length);

    // Memory allocation
    NumericMatrix phi(m_N, length);
    NumericMatrix delta(m_N, length);    
        
    //  log-transform to avoid double precision loss
    NumericMatrix A(m_N, m_N);    
    NumericVector Pi(m_N);

    // Used to find the max value
    NumericVector temp(m_N);    

    for(i = 0 ; i < m_N; i++)
    {
        Pi[i] = log(m_Pi[i]);
        for(j = 0; j < m_N; j++)
            A(i,j) = log(m_A(i,j)); 
    }

    //  Init step
    for(i = 0; i < m_N; i++)
        delta(i,0)= Pi[i] + R::dpois(sequence[0], m_B[i], true);


    //  Recursion step
    for(j = 1; j < length; j++)
        for(i = 0; i < m_N; i++)
        {
            for(k = 0; k < m_N; k++)
                temp[k] = delta(k,j-1) + A(k,i);
            // The auto keyword is simply asking the compiler to deduce the type of the variable from the initialization                
            auto maximum = max_element(temp.begin(), temp.end());             
            delta(i,j) = (*maximum) + 
                    R::dpois(sequence[j], m_B[i], true);            
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
double HMMpoisson::loglikelihood(IntegerMatrix sequences)
{
    double ll = 0.0;
    unsigned int i, seqLen;

    seqLen = sequences.nrow();

    for(i = 0; i < seqLen; i++)
        ll+=evaluation(sequences.row(i));
    return ll;
}

//  Parameter estimation using a Expectation Maximization approach for multiple sequences
void HMMpoisson::learnEM(IntegerMatrix sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{

    double newLL, error;
    double lastLL = loglikelihood(sequences);
    unsigned int counter = 0; 
    double max = 0.0, min = 0.0;    
    NumericVector seqRow;

    //  We search the min and max value, if a reinitialization is needed
    for(int i = 0; i < sequences.nrow(); i++)
    {
        seqRow = sequences.row(i);
        auto tmin =  min_element ( seqRow.begin(), seqRow.end() );
        auto tmax =  max_element ( seqRow.begin(), seqRow.end() );
        if(*tmin < min)
            min = *tmin;            
        if(*tmax > max)
            max = *tmax;
    }
    
    //Validate the sequence
    if(min < 0) 
        Rf_error("All values in the sequnce must be positive");
    
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
    } while(counter < iter && error > delta);
    Rcout << "Finished at Iteration: " << counter << " with Error: " << error  << "\n";
}

//  Parameter estimation using a Baum Welch approach for a single sequence
void HMMpoisson::learnBW(IntegerVector sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print )
{
    double newLL, error;
    double lastLL = evaluation(sequences);
    unsigned int counter = 0; 
    double max = 0.0, min = 0.0;    

    //  We search the min and max value, if a reinitialization is needed
    auto tmin =  min_element ( sequences.begin(), sequences.end() );
    auto tmax =  max_element ( sequences.begin(), sequences.end() );
    if(*tmin < min)
        min = *tmin;            
    if(*tmax > max)
        max = *tmax;

    //Validate the sequence
    if(min < 0) 
        Rf_error("All values in the sequnce must be positive");        

        
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
List HMMpoisson::generateObservations( unsigned short int length)
{
    unsigned int i,j;
    double x;    

    IntegerVector X(length, 0);
    IntegerVector Y(length, 0);

    //  Used for set.seed compatibility
    RNGScope scope; 

    //  Matrix rearrangement to use a uniform distribution to generate the hidden states and its corresponding observation
    NumericMatrix A(m_N, m_N);     
    NumericVector Pi(m_N);

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

    //  Random variable feneration based in the rearranged matrices
    x = as<double>(runif(1));       
    X[0] = lower_bound (Pi.begin(), Pi.end(), x) - Pi.begin();
    Y[0] = as<double>(rpois(1, m_B[X[0]]));     
      
    NumericVector tempA;

    for(j = 1; j < length; j++)
    {
        x = as<double>(runif(1));

        tempA = A.row(X[j-1]);     
        X[j] = lower_bound (tempA.begin(), tempA.end(), x) - tempA.begin();
        Y[j] = as<double>(rpois(1, m_B[X[j]]));          
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
List HMMpoisson::toList(void) const
{
    return List::create(
            Named("Model", "PHMM"),
            Named("StateNames", getStateNames() ),
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
void HMMpoisson::forwardMatrix( IntegerVector sequence, unsigned int length , scaledMatrix & forward)
{  
    unsigned int i, j, k;   

    //  Base case and forward matrix initialization
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0)= R::dpois(sequence[0], m_B[i], false) * m_Pi[i];
        forward.scaling[0] += forward.matrix(i,0);
    } 
        
    //  First column normalization (sum = 1)
    for( i = 0; i < m_N; i++)
        forward.matrix(i,0)/=  forward.scaling[0];
    //  Recursive step
    //  It starts at the second observation
    for(j = 1; j < length; j++)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
                forward.matrix(i,j) += m_A(k,i)*forward.matrix(k,j-1);
            forward.matrix(i,j) *= R::dpois(sequence[j], m_B[i], false);
            forward.scaling[j] += forward.matrix(i,j);        
        }

        //  Factor normalization
        for(i = 0; i < m_N; i++ )
            forward.matrix(i,j) /=  forward.scaling[j];
    }   
}

//  Backward
void HMMpoisson::backwardMatrix( IntegerVector sequence, unsigned int length , scaledMatrix & backward)
{  
    unsigned int i, j, k;

    //  Base case and backward matrix initialization
    for( i = 0; i < m_N; i++)          
        backward.matrix(i,length-1) = 1;    

    //  Recursive step 
    for(j = length - 1 ; j > 0  ; j--)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
                backward.matrix(i,j-1) += R::dpois(sequence[j], m_B[k], false) * m_A(i,k) * backward.matrix(k,j);            
            backward.scaling[j] += backward.matrix(i,j-1);      
        }

        //  Factor normalization                 
        for(i = 0; i < m_N; i++ )
            backward.matrix(i,j-1) /=  backward.scaling[j];
    } 
    
    //  Last step
    for(i = 0; i < m_N; i++ )
        backward.scaling[0] += m_Pi[i] * R::dpois(sequence[0], m_B[i], false) * backward.matrix(i,0); 
}

//--------------------------------------------------------------------------
// DECODING
//--------------------------------------------------------------------------

//  Function dedicated to memory allocation for the forward backward algorithm
NumericMatrix HMMpoisson::forwardBackwardGamma(IntegerVector sequence)
{           
    unsigned int length = sequence.size();       

    //  scaling factors for the forward and backward matrices
    NumericVector scaledf(length, 0);
    NumericVector scaledb(length + 1, 0);  //+1 para tener el valor prior de 1
    scaledb[length] = 0; //log(1) = 0

    //  Memory reserved for each matrix. The matrices are cloned to make a "safe" memory management
    NumericMatrix matrix(m_N, length);
    scaledMatrix forward = {clone(scaledf), clone(matrix)};
    scaledMatrix backward = {clone(scaledb),  clone(matrix)};    

    //  Algorithm call
    forwardBackwardGamma(sequence, forward, backward, scaledf, scaledb, matrix, length);

    //  Gamma matrix
    return matrix;
}

void HMMpoisson::forwardBackwardGamma(IntegerVector sequence, scaledMatrix & forward, scaledMatrix & backward,  NumericVector & scaledf, NumericVector & scaledb, NumericMatrix & matrix, unsigned int length)
{
    unsigned int i, j, k;
    double eval;              

    //  Initial step   
    for( i = 0; i < m_N; i++)
    {           
        forward.matrix(i,0)= R::dpois(sequence[0], m_B[i], false) * m_Pi[i];
        forward.scaling[0] += forward.matrix(i,0);
        backward.matrix(i,length-1) = 1;   
    } 
        
    //  First factor normalization (sum = 1)
    for( i = 0; i < m_N; i++)
        forward.matrix(i,0)/=  forward.scaling[0];

    //  Recursive step   
    for(j = 1; j < length; j++)  
    {        
        for(i = 0; i < m_N; i++ )
        {
            for(k = 0; k < m_N; k++)
            {
                forward.matrix(i,j) += m_A(k,i)*forward.matrix(k,j-1);
                backward.matrix(i,length-j-1) += 
                    R::dpois(sequence[length - j], m_B[k], false) * m_A(i,k) * backward.matrix(k,length-j);
            }                
            forward.matrix(i,j) *=  R::dpois(sequence[j], m_B[i], false);
            
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
        backward.scaling[0] += m_Pi[i] * R::dpois(sequence[0], m_B[i], false) * backward.matrix(i,0);

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
void HMMpoisson::expectationMaximization( IntegerMatrix sequences, unsigned int pseudo)
{
    unsigned int seqLen, length, i, j, k, s;   
    double temp;

    seqLen = sequences.nrow();
    length = sequences.ncol();       

    //  Memory allocation
    NumericMatrix Amean(m_N, m_N);
    NumericVector Bmean(m_N);                    
    NumericVector Pimean(m_N);

    //  Normalizing factors
    NumericVector denomA(m_N);
    NumericVector denomB(m_N);                      
         
    //  Analysis per sequence  
    for(s = 0; s < seqLen ; s++)
    {        
        //  Expectation step
        //  Memory allocation for expectation step
        NumericVector scaledf(length);
        NumericVector scaledb(length+ 1);

        NumericMatrix matrix(m_N, length); //  Gamma matrix       
        scaledMatrix forward = {clone(scaledf), clone(matrix)};
        scaledMatrix backward = {clone(scaledb),  clone(matrix)};
        
        forwardBackwardGamma(sequences.row(s), forward, backward, scaledf, scaledb, matrix, length);

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
                    temp = (matrix(i,j)*m_A(i,k)*
                        R::dpois(sequences(s,j+1), m_B[k], false) *
                        backward.matrix(k,j+1))/
                        (backward.matrix(i,j)* backward.scaling[j+1] );  // j+1 is the shift in the backward scaling vector
                    Amean(i,k) += temp;                    
                    denomA[i] += temp;
                }
                
                //  For the emission matrix, we can use P(X|data)
                Bmean[i] += (matrix(i,j)*sequences(s,j));  //  Lambda                
                denomB[i] += matrix(i,j);
            }
            //  The loop ends at 'length-1', thus it is necessary to sum the las value
            Bmean[i] += (matrix(i,length-1)*sequences(s,length-1));  //Lambda
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
        m_B[i] = Bmean[i] / denomB[i];  //  Normalized mean  
    }   
}

//  Function used for parameter estimation given a single sequence.
//  The initial probability vector is not estimated
void HMMpoisson::BaumWelch( IntegerVector sequences, unsigned int pseudo)
{
    unsigned int length, i, j, k;   
    double temp;
    length = sequences.size();      

    //  Memory allocation
    NumericMatrix Amean(m_N, m_N);
    NumericVector Bmean(m_N);                        

    //  Normalizing factors
    NumericVector denomA(m_N);
    NumericVector denomB(m_N);                     
         
    //  Expectation step
    //  Memory allocation for expectation step
    NumericVector scaledf(length);
    NumericVector scaledb(length+ 1);

    NumericMatrix matrix(m_N, length); //  Gamma matrix        
    scaledMatrix forward = {clone(scaledf), clone(matrix)};
    scaledMatrix backward = {clone(scaledb),  clone(matrix)};
        
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
                temp = (matrix(i,j)*
                        m_A(i,k)*
                        R::dpois(sequences[j+1], m_B[k], false)*
                        backward.matrix(k,j+1))/
                        (backward.matrix(i,j) * backward.scaling[j+1] ); // j+1 is the shift in the backward scaling vector
                Amean(i,k) += temp;                
                denomA[i] += temp;
            }
                
            //  Emission matrix estimation
            Bmean[i] += (matrix(i,j)*sequences[j]);  //  Lambda      
            denomB[i] += matrix(i,j);
        }
        //  The loop ends at 'length-1', thus it is necessary to sum the las value
        Bmean[i] += (matrix(i,length-1)*sequences[length-1]);  // Lambda
        denomB[i] += matrix(i,length-1);  
    }                            
           
    //  We use the normalizing factor and also use the pseudo count value
    for(i = 0; i < m_N; i++)
    {
        //  Transition matrix normalization
        for(k = 0; k < m_N; k++)
            m_A(i,k) = (Amean(i,k) + pseudo) / (denomA[i] + m_N*pseudo);             

        //  Emission matrix normalization
        m_B[i] = Bmean[i] / denomB[i];      //  Normalized mean        
    }
}

//--------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------

//  Function used to set random model parameters
void HMMpoisson::randomInit(double min, double max)
{

    //  Used for set.seed compatibility
    RNGScope scope;  

    //  Normalizing factors
    double maxPi = 0.0;
    NumericVector maxA(m_N);  

    for (int i=0; i< m_N ; i++) 
    {
        //  Initial probability vector
        m_Pi[i] = as<double>(runif(1));
        maxPi += m_Pi[i] ;
        maxA[i] = 0;        
        //  Transition matrix 
        for(int j = 0; j < m_N; j++)
        {
            m_A(i,j) = as<double>(runif(1));
            maxA[i] += m_A(i,j);
        }
        //  Emission matrix
        // One parameter: lambda
        m_B[i] = as<double>(runif(1, min, max));
    }

    //  Parameter normalization
    for(int i = 0; i < m_N; i++)
    {
        m_Pi[i] /= maxPi;
        for(int j = 0; j < m_N; j++)
            m_A(i,j) /= maxA[i];
    }  
}

