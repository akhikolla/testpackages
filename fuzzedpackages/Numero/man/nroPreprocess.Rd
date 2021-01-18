\name{nroPreprocess}
\alias{nroPreprocess}

\title{
Data cleaning and standardization
}

\description{
Convert to numerical values, remove unusable rows and columns, and
standardize scale of each variable.
}

\usage{
nroPreprocess(data, method = "standard", clip = 3.0,
    resolution = 100, trim = FALSE)
}

\arguments{
  \item{data}{
A matrix or a data frame. 
  }

  \item{method}{
Method for standardizing scale and location, see details below.
  }

  \item{clip}{
Range for clipping extreme values in multiples of standard deviations.
  }
  
  \item{resolution}{
Maximum number of sampling points to capture distribution shape.
  }
  
  \item{trim}{
if TRUE, empty rows and columns are removed.
  }
}

\details{
Standardization methods include empty string for no action, "standard"
for centering by mean and division by standard deviation, "uniform" for
normalized ranks between -1 and 1 and "tapered" for a version of the
rank-based method that puts more samples around zero.

Clipping is not applied if the method is rank-based.
}

\value{
A matrix of numerical values. A value mapping model is stored in the
attribute 'mapping'. The names of binary columns are stored in the
attribute 'binary'.
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

# Detect binary columns.
ds <- nroPreprocess(dataset, method = "")
print(attr(ds,"binary"))

# Centering and scaling cholesterol.
ds <- nroPreprocess(dataset$CHOL)
print(summary(ds))

# Centering and scaling.
ds <- nroPreprocess(dataset)
print(summary(ds))

# Tapered ranks.
ds <- nroPreprocess(dataset, method = "tapered")
print(summary(ds))
}