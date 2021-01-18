\name{sepsis_dataset_test}
\alias{sepsis_dataset_test}
\title{Sepsis dataset test}
\usage{
sepsis_dataset_test(induce_biomarker_decision_tree = TRUE,
  threshold = 75, debug = FALSE, hsc = NULL)
}
\arguments{
\item{induce_biomarker_decision_tree}{(logical) - If FALSE, "Biomarker assessment" is one activity ignoring biomarker values.
If TRUE, based on the biomarker values, several distinct "Biomarker assessment" activities are inferred.}

\item{threshold}{(numeric) - Projection threshold.}

\item{debug}{(logical) - Switch for debug printout.}

\item{hsc}{(HybridSequenceClassifier) - An existing Sequence Detector that should be used instead of creating a new one.}
}
\value{
None
}
\description{
\link[eventdataR]{sepsis} dataset is taken from the package \emph{eventdataR} and used to test the Sequence Detector. The results are
available in [1].
}
\references{
[1] D. Krleža, B. Vrdoljak, and M. Brčić, Latent Process Discovery using Evolving Tokenized Transducer, \emph{IEEE Access}, vol. 7, pp. 169657 - 169676, Dec. 2019
}
