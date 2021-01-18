\name{ISOpure.calculate.tac}
\alias{ISOpure.calculate.tac}
\title{Perform calculation for Tumour Adjacent Cell (TAC) profiles}
\description{Performs the mathematical calculations taking bulk tumor data and deconvolved profiles
	and returning deconvolved tumour adjacent cell profiles.}
\usage{
ISOpure.calculate.tac(tumor.profiles, deconvolved.profiles, purity.estimates)
}
\arguments{

	\item{tumor.profiles}{a GxD matrix representing gene expression profiles of heterogeneous (mixed) 
		tumor samples, where G is the number of genes, D is the number of tumor samples.}

	\item{deconvolved.profiles}{a GxD matrix representing gene expression profiles of purified (ISOpure output)
		tumor samples, where G is the number of genes, D is the number of tumor samples.}

	\item{purity.estimates}{a vector D representing the purity estimates (output from ISOpure)}

}
\value{a GxD matrix representing gene expression profiles of purified (ISOpure output)
	tumor adjacent cell signal, where G is the number of genes, D is the number of tumor samples. 
    }
\author{Natalie Fox}
