\name{deeplearning-methods}
\alias{Rda}
\alias{Rsda}
\alias{Rrbm}
\alias{Rdbn}
\alias{train}
\alias{pretrain}
\alias{finetune}
\alias{predict}
\alias{reconstruct}

\title{
  Create \code{deeplearning} objects from training set.
}
\usage{
 Rda(x)
 
 Rsda(x, y, hidden)
 
 Rrbm(x)
 
 Rdbn(x, y, hidden)
 
 train(object)
 
 pretrain(object)
 
 finetune(object)
 
 reconstruct(object, test)
 
 predict(object, test)
}
\arguments{
  \item{x}{The training dataset.}
  \item{y}{The labels for training dataset.}
  \item{test}{The testing dataset.}
  \item{hidden}{The number of hidden representation in each layer.}
  \item{object}{ An instantiated \code{deeplearning} object. }
}

\description{
  \code{Rda}, \code{Rsda}, \code{Rrbm} and \code{Rdbn} will return an instantiated \code{deeplearning} object
  for denoising autoencoder, stacked denoising autoencoder, restricted Boltzmann machine and deep belief net.
  \code{train} and \code{reconstruct} are for training and reconstructing from denoising autoencoder and restricted Boltzmann machine;
  \code{pretrain}, \code{finetune} and \code{predict} are used for pretraining, finetuning and predicting
  using stacked denoising autoencoder and deep belief net.
}

\author{
  Qiang Kou
}

\examples{
 data(test)
 dbn_test <- Rdbn(train_X, train_Y, hidden)
 summary(dbn_test)
 LearningRate(dbn_test)
 pretrain(dbn_test)
 finetune(dbn_test)
 predict(dbn_test, test_X)
}
