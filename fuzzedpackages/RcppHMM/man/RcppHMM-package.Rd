\name{RcppHMM-package}
\alias{RcppHMM-package}
\alias{RcppHMM}

\title{
  Overview of Package RcppHMM
}

\description{
  This package can model observations based on hidden Markov models. The observations can be considered to be emitted by a multinomial distribution, A mixture of Gaussians or a mixture of Poissons. It can be used for inference, parameter estimation and simulation.  
}

\details{
  The package can be used to represent a discrete-time hidden Markov model. The states can generate categorical (labeled), continuous or discrete observations. The hidden state transition and observations can be randomly generated based on fixed parameters. Also, the inference methods can be used to evaluate sequences or decode the hidden states that generated the observations. Finally, the model parameters can be estimated by a single or multiple observed sequences.  
}

\author{
Roberto A. Cardenas-Ovando, Julieta Noguez and Claudia Rangel-Escareno

Maintainer: Roberto A. Cardenas-Ovando <robalecarova@gmail.com>
}

\references{

Bilmes, J.E. (1998). A Gentle Tutorial of the EM Algorithm and its Application to Parameter Estimation for Gaussian Mixture and Hidden Markov Models. \emph{International Computer Science Institute}.

Ibe, O. (2009). Markov processes for stochastic modeling. \emph{Oxford}.

Rabiner, L.R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. \emph{Proceedings of the IEEE}.

Rabiner L.; Juang, B.H. (1993) Fundamentals of Speech Recognition. \emph{Prentice Hall Signal Processing Series}.

}

\examples{
# Multinomial case
# Set the model parameters to be estimated
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

# Model parameters validation

HMM <- verifyModel(params)

# Data simulation
# Multiple sequences

set.seed(100)
length <- 100
seqs <- 100
observationSequences<- c()
for(i in 1:seqs){
  Y <- generateObservations(HMM , length)$Y
  observationSequences <- rbind(observationSequences , Y)
}

# New model random initialization

set.seed(1000)
newModel <- initHMM(2,4) 
n = c("X1","X2")
m = c("A","T","C","G")

# Change model names

newModel <- setNames(newModel,
                        list( "StateNames" = n,
                              "ObservationNames" = m) )

# Model parameters estimation

newModel <- learnEM(newModel,
        observationSequences,
        iter=300, 
        delta = 1E-8,
        pseudo = 0,
        print = TRUE)

# New sequence simulation to compare the new model
# Data simulation

# Single sequence
Y <- generateObservations(HMM , length)$Y

# Evaluation

evaluation(newModel, Y, "f")
evaluation(newModel, Y, "b")

# Hidden state decoding

hiddenStatesViterbi <- viterbi(newModel, Y)
hiddenStatesFB <- forwardBackward( newModel, Y)
}

