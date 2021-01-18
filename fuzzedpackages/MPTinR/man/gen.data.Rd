\name{gen.data}
\alias{gen.data}
\alias{sample.data}
\alias{gen.predictions}
\title{
Generate or bootstrap data and get predictions from a model specified in a model file (or connection).
}
\description{
\code{gen.data} generates random dataset(s) from given paramater values and model (specified via model file or textConnection) for paramteric bootstrap.\cr
\code{sample.data} generates random dataset(s) from given data for nonparametric bootstrap.\cr
\code{gen.predictions} generates response probabilities or predicted responses from given paramater values and model (specified via model file or textConnection). \cr
}
\usage{
gen.data(parameter.values, samples,
	model.filename, 
	data = NULL, n.per.item.type = NULL,
	restrictions.filename = NULL, model.type = c("easy", "eqn", "eqn2"), 
	reparam.ineq = TRUE, check.model = TRUE)

sample.data(data, samples,
	model.filename = NULL, categories.per.type = NULL,
	model.type = c("easy", "eqn", "eqn2"), check.model = TRUE)	

gen.predictions(parameter.values, 
	model.filename, 
	restrictions.filename = NULL, 
	n.per.item.type = NULL, 
	model.type = c("easy", "eqn", "eqn2"), 
	reparam.ineq = TRUE, check.model = TRUE)
}


\arguments{
  \item{parameter.values}{
\code{vector} of paramater values. Either named then order is irrelevant or unnamed then must follow the alphabetical order of paramaters (\code{\link{check.mpt}} returns the alphabetical order of paramater names).
}
  \item{samples}{
Number of random datasets to be generated from a given set of paramater values.
}
  \item{n.per.item.type}{
\code{vector} of length equal to number of item types (or trees) specifying how many item each item type has. Default is \code{NULL}. See Details.
} 
  \item{data}{
data \code{vector}. See Details.
}
  \item{categories.per.type}{
numeric vector indicating how many response categories each item type has.
} 
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file, pssoble a \code{\link{textConnection}}. 
}
  \item{restrictions.filename}{
\code{NULL} or a character \code{vector} or a \code{list} of characters. The default is \code{NULL} which corresponds to no restrictions. A character \code{vector} specifies the location or name of the restrictions file. A \code{list} of characters contains the restrictions directly. See \code{\link{fit.mpt}} for Details and Examples.
}
  \item{model.type}{
Character vector specifying whether the model file is formatted in the easy way (\code{"easy"}; i.e., each line represents all branches corresponding to a response category) or the traditional EQN syntax (\code{"eqn"} or \code{"eqn2"}; see Details and e.g., Stahl & Klauer, 2007). If \code{model.filename} ends with .eqn or .EQN, \code{model.type} is automatically set to \code{"eqn"}. Default is \code{"easy"}.
}
  \item{reparam.ineq}{
Should inequality restrictions be applied (i.e., the model reparametrized)? Default is \code{TRUE}.
}
  \item{check.model}{
logical. Should model be chekced with random values whether or not the expected values sum to one per tree? Default is \code{TRUE}.
}
}

\details{
\code{gen.data} and \code{sample.data} are basically wrapper for \code{\link{rmultinom}} (called multiple times, if there is more than one item type). The \code{prob} argument of \code{rmultinom} is obtained differently for the two functions. For \code{gen.data} it corresponds to the predicted response proportions as returned by \code{get.predictions} (which is actually called by \code{gen.data}). For \code{sample.data} it is the proprtion of responses for each item type.

\code{gen.data} needs to know how big the n for each item type is. This can either be specified via the \code{data} or the \code{n.per.item.type} argument (i.e., one of those needs to be non-\code{NULL}). See the examples.

\code{sample.data} needs to know which response categories correspond to each item type. This can either be specified via the \code{model.filename} or the \code{categories.per.type} argument (i.e., one of those needs to be non-\code{NULL}). See the examples.

}
\value{
Either a vector or matrix containing the generated data (for \code{gen.data} and \code{sample.data}) or a vector containing the predictions (for \code{gen.predictions}).
}

\author{
Henrik Singmann and David Kellen
}

\seealso{
\code{\link{fit.mpt}} or \code{\link{fit.model}} for functions that will fit the generated data. Note that it is probably a very good idea to set \code{fit.aggregated = FALSE} when fitting larger sets of generated data.
}
\examples{
#### using the model and data from Broeder & Schuetz:
data(d.broeder, package = "MPTinR")
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
m.sdt <- "pkg/MPTinR/inst/extdata/broeder.sdt.model"

m.sdt <- system.file("extdata", "broeder.sdt.model", package = "MPTinR")

# fit the 2HTM
br.2htm <- fit.mpt(colSums(d.broeder), m.2htm)

# fit the SDT model
br.sdt <- fit.model(colSums(d.broeder), m.sdt, lower.bound = c(rep(-Inf, 5), 0, 1),
	upper.bound = Inf)

# get one random dataset using the paramater values obtained (i.e., parametric bootstrap) 
# and the data argument.
gen.data(br.2htm[["parameters"]][,1], 1, m.2htm, data = colSums(d.broeder))

gen.data(br.sdt[["parameters"]][,1], 1, m.sdt, data = colSums(d.broeder))

# get one random dataset using the paramater values obtained (i.e., parametric bootstrap) 
# and the n.per.item.type argument.
gen.data(br.2htm[["parameters"]][,1], 1, m.2htm, 
	n.per.item.type = c(240, 2160, 600, 1800, 1200, 1200, 1800, 600, 2160, 240))

gen.data(br.sdt[["parameters"]][,1], 1, m.sdt, 
	n.per.item.type = c(240, 2160, 600, 1800, 1200, 1200, 1800, 600, 2160, 240))

# sample one random dataset from the original data:
sample.data(colSums(d.broeder), 1, model.filename = m.2htm)
# above uses the model.filename argument

sample.data(colSums(d.broeder), 1, categories.per.type = rep(2,10))
# above uses the categories.per.type argument


# just get the predicted proportions:
predictions.mpt <- gen.predictions(br.2htm[["parameters"]][,1], m.2htm)
predictions.sdt <- gen.predictions(br.sdt[["parameters"]][,1], m.sdt)

# predicting using the proactive Inhibiton Model (Riefer & Batchelder, 1988, Figure 1)

model1 <- system.file("extdata", "rb.fig1.model", package = "MPTinR")

gen.predictions(c(r = 0.3, p = 1, q = 0.4944), model1)  
gen.predictions(c(r = 0.3, p = 1, q = 0.4944), model1, n.per.item.type = 180)

# the order of parameters is reordered (i.e., not alphabetically)
# but as the vector is named, it does not matter!
# Compare with:
data(rb.fig1.data, package = "MPTinR")
fit.mpt(rb.fig1.data[1,], model1, n.optim = 1)
}

