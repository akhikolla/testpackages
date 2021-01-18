\name{numero.prepare}
\alias{numero.prepare}

\title{Prepare datasets for analysis}

\description{Prepare training data by mitigating confounding factors and
standardizing values.}

\usage{
numero.prepare(data, variables = NULL, confounders = NULL, batch = NULL,
               method = "standard", pipeline = NULL)
}

\arguments{
  \item{data}{
A matrix or a data frame.
  }
  \item{variables}{
A character vector of column names.  
  }
  \item{confounders}{
Names of columns that contain confounder data.
  }
  \item{batch}{
The name of the column that contains batch labels.
  }
\item{method}{
Method to standardize values, see \code{\link{nroPreprocess}()}.
}
\item{pipeline}{
Processing parameters from a previous use of the function.
}
}

\details{
We recommend first applying \code{\link{numero.clean}()} to the full
dataset, then selecting a subset for training using the input argument
\code{variables}. This preserves any attributes that may be used in
Numero functions.

If a previous \code{pipeline} is available, it overrides all processing
parameters irrespective of other input arguments.

Due to safeguards against numerical instability, the standardized values
may deviate slightly from the expected range (<0.1 percent error is typical).
}

\value{
A matrix with the attributes 'pipeline' that contains the processing
parameters and 'subsets' that contains row names divided into batches if
batch correction was applied.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Set identities and manage missing data.
dataset <- numero.clean(dataset, identity = "INDEX")

# Prepare training variables using default standardization.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- numero.prepare(data = dataset, variables = trvars)
print(summary(trdata))

# Prepare training values adjusted for age and sex and
# standardized by rank-based method.
trdata <- numero.prepare(data = dataset, variables = trvars,
                         batch = "MALE", confounders = "AGE",
			 method = "tapered")
print(summary(trdata))
}