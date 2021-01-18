\name{synthetic_test_agenda}
\alias{synthetic_test_agenda}
\title{Synthetic process test}
\usage{
synthetic_test_agenda(label_aspect=1)
}
\arguments{
\item{label_aspect}{(numeric) - A vector of all relevant consolidated data stream fields.}
}
\description{
A synthetic process that was introduced in the process mining agenda [1].
The original event log introduced in [1] did not comprise any timestamps, and a process discovery algorithm was intended to infer this based on the event position in the log.
ETT and new process discovery algorithms require events to have at least some sort of timing, and this was added for this test.
It is worth noticing that the given event log has some parallel activities, which should be detected by the process discovery algorithm.
The final results of this test are described in [2].
}
\references{
[1] W. M. P. van der Aalst and A. J. M. M. Weijters, Process mining: a research agenda, \emph{Computers in Industry}, vol. 53, no. 2, pp. 231--244, Apr. 2004\cr 
[2] D. Krleža, B. Vrdoljak, and M. Brčić, Latent Process Discovery using Evolving Tokenized Transducer, \emph{IEEE Access}, vol. 7, pp. 169657 - 169676, Dec. 2019
}
