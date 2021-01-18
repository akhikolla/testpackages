\name{nroPostprocess}
\alias{nroPostprocess}

\title{Standardization using existing parameters}

\description{Process a new dataset using a standardization procedure that
was created for another dataset}

\usage{
nroPostprocess(data, mapping, reverse = FALSE, trim = FALSE)
}

\arguments{
  \item{data}{
A matrix or a data frame with column names.
  }

  \item{mapping}{
A list object or a matrix or a data frame. 
  }
  
  \item{reverse}{
If true, standardized data will be reverted back to original scale. 
  }
  
  \item{trim}{
If true, unusable rows and columns are removed. 
  }
}

\details{
The input argument can be a data frame with the attribute 'mapping' as
returned from \code{\link{nroPreprocess}()} or a list object with the
elements \code{input} and \code{output} that each contain a data frame or
a matrix of equal size.

The function projects the input data to the values in \code{mapping$input}
to determine the positions of the input values with respect to the rows
in the model. These positions are then used to interpolate corresponding
output values in \code{mapping$output}.

The mapping elements must have consistent row and column names.
}

\value{
A matrix or data frame of processed values.
}

\author{
Ville-Petteri Makinen
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Show original data characteristics.
print(summary(dataset))

# Preprocess a subset of data.
ds.pre <- nroPreprocess(dataset[1:100,])
print(summary(ds.pre))

# Repeat preprocessing for the whole dataset (approximation).
ds.post <- nroPostprocess(dataset, ds.pre)
print(summary(ds.post))
}