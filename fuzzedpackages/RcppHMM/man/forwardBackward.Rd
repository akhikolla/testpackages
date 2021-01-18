\name{forwardBackward}
\alias{forwardBackward}

\title{
Forward-backward algortihm for hidden state decoding
}

\description{
Function used to get the most likely hidden states at each observation in the provided sequence.
}

\usage{
forwardBackward(hmm, sequence)
}

\arguments{
  \item{hmm}{a list with the necessary variables to define a hidden Markov model.}
  \item{sequence}{sequence of observations to be decoded. HMM and PHMM use a vector. GHMM uses a matrix.}
}

\details{
GHMM uses a matrix with the variables as rows and consecutive observations in the columns. 
}

\value{
A vector of hidden states in the traveled path of observations.
}

\references{
Cited references are listed on the \link{RcppHMM} manual page.
}

\seealso{
\code{\link{generateObservations}} , \code{\link{verifyModel}} , \code{\link{viterbi}}
}

\examples{
## Values for a hidden Markov model with categorical observations
# Set the model parameters
n <- c("First","Second")
m <- c("A","T","C","G")
A <- matrix(c(0.8,0.2,
              0.1,0.9),
            nrow = 2,
            byrow = TRUE)

B <- matrix(c(0.2, 0.2, 0.3, 0.3,
              0.4, 0.4, 0.1, 0.1),
            nrow = 2,
            byrow = TRUE)

Pi <- c(0.5, 0.5)

params <- list( "Model" = "HMM",
                "StateNames" = n,
                "ObservationNames" = m,
                "A" = A,
                "B" = B,
                "Pi" = Pi)

HMM <- verifyModel(params)

# Data simulation
set.seed(100)
length <- 100
observationSequence <- generateObservations(HMM, length)

#Sequence decoding
hiddenStates <- forwardBackward(HMM, observationSequence$Y)
print(hiddenStates)

## Values for a hidden Markov model with discrete observations

n <- c("Low","Normal","High")

A <- matrix(c(0.5, 0.3,0.2,
              0.2, 0.6, 0.2,
              0.1, 0.3, 0.6),
            ncol=length(n), byrow=TRUE)

B <- c(2600,  # First distribution with mean 2600
       2700,  # Second distribution with mean 2700
       2800)  # Third distribution with mean 2800

Pi <- rep(1/length(n), length(n))

HMM.discrete <- verifyModel(list("Model"="PHMM", "StateNames" = n, "A" = A, "B" = B, "Pi" = Pi))

# Data simulation
set.seed(100)
length <- 100
observationSequence <- generateObservations(HMM.discrete, length)

#Sequence decoding
hiddenStates <- forwardBackward(HMM.discrete, observationSequence$Y)
print(hiddenStates)

## Values for a hidden Markov model with continuous observations                          
# Number of hidden states = 3
# Univariate gaussian mixture model

N = c("Low","Normal", "High")
A <- matrix(c(0.5, 0.3,0.2,
              0.2, 0.6, 0.2,
              0.1, 0.3, 0.6),
            ncol= length(N), byrow = TRUE)

Mu <- matrix(c(0, 50, 100), ncol = length(N))
Sigma <- array(c(144, 400, 100), dim = c(1,1,length(N)))
Pi <- rep(1/length(N), length(N))

HMM.cont.univariate <- verifyModel(list( "Model"="GHMM", 
                                         "StateNames" = N,
                                         "A" = A, 
                                         "Mu" = Mu, 
                                         "Sigma" = Sigma, 
                                         "Pi" = Pi))

# Data simulation
set.seed(100)
length <- 100
observationSequence <- generateObservations(HMM.cont.univariate, length)

#Sequence decoding
hiddenStates <- forwardBackward(HMM.cont.univariate, observationSequence$Y)
print(hiddenStates)

## Values for a hidden Markov model with continuous observations                          
# Number of hidden states = 2
# Multivariate gaussian mixture model
# Observed vector with dimensionality of 3
N = c("X1","X2")
M <- 3

# Same number of dimensions
Sigma <- array(0, dim =c(M,M,length(N)))
Sigma[,,1] <- matrix(c(1.0,0.8,0.8,
                       0.8,1.0,0.8,
                       0.8,0.8,1.0), ncol = M,  
                     byrow = TRUE)
Sigma[,,2] <- matrix(c(1.0,0.4,0.6,
                       0.4,1.0,0.8,
                       0.6,0.8,1.0), ncol = M,
                     byrow = TRUE)
Mu <- matrix(c(0, 5, 
               10, 0, 
               5, 10), 
             nrow = M, 
             byrow = TRUE)

A <- matrix(c(0.6, 0.4,
              0.3, 0.7), 
            ncol = length(N),
            byrow = TRUE)
Pi <- c(0.5, 0.5)

HMM.cont.multi <- verifyModel(list( "Model" = "GHMM",
                                    "StateNames" = N,
                                    "A" = A, 
                                    "Mu" = Mu, 
                                    "Sigma" = Sigma, 
                                    "Pi" = Pi))

# Data simulation
set.seed(100)
length <- 100
observationSequence <- generateObservations(HMM.cont.multi, length)

#Sequence decoding
hiddenStates <- forwardBackward(HMM.cont.multi, observationSequence$Y)
print(hiddenStates)

}

\keyword{methods}
