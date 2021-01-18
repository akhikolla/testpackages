\name{metadata}
\alias{PretrainLearningRate}
\alias{CorruptionLevel}
\alias{PretrainingEpochs}
\alias{FinetuneLearningRate}
\alias{FinetuneEpochs}
\alias{Step}
\alias{HiddenRepresentation}
\alias{LearningRate}
\alias{TrainingEpochs}
\alias{setCorruptionLevel}
\alias{setFinetuneEpochs}
\alias{setFinetuneLearningRate}
\alias{setHiddenRepresentation}
\alias{setLearningRate}
\alias{setPretrainEpochs}
\alias{setPretrainLearningRate}
\alias{setStep}
\alias{setTrainingEpochs}

\title{
  Access the metadata from a \code{deeplearning} object.
}

\usage{
LearningRate(object)
PretrainLearningRate(object)
CorruptionLevel(object)
TrainingEpochs(object)
PretrainingEpochs(object)
FinetuneLearningRate(object)
FinetuneEpochs(object)
Step(object)
HiddenRepresentation(object)
setCorruptionLevel(object, x)
setFinetuneEpochs(object, x)
setFinetuneLearningRate(object, x)
setHiddenRepresentation(object, x)
setLearningRate(object, x)
setPretrainEpochs(object, x)
setPretrainLearningRate(object, x)
setStep(object, x)
setTrainingEpochs(object, x)
}
\arguments{
  \item{object}{An instantiated \code{deeplearning} object.}
  \item{x}{A numeric number to set.}
}

\description{
  Accessors to the metadata of a \code{deeplearning} object.
  \code{summary} will show a summary of the paramters as a named list,
  and other functions will return corresponding features.
}

\author{
  Qiang Kou
}
