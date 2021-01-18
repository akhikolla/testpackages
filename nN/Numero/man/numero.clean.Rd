\name{numero.clean}
\alias{numero.clean}

\title{Clean datasets}

\description{Sets row names and removes unusable columns and rows.}

\usage{
numero.clean(..., identity = NULL, na.freq = 0.9,
             num.only = TRUE, select = "")
}

\arguments{
  \item{...}{
Matrices or a data frames. 
  }
  \item{identity}{
Name(s) of the column(s) that contain identification information.
  }
  \item{na.freq}{
The proportion of how many missing values are allowed in
each column and in each row.
  }
  \item{num.only}{
If true, only numeric columns are included.
  }
  \item{select}{
Indicate if only identities present in all datasets
or in exactly one of the datasets are included.
  }
}

\details{
If multiple identity columns are provided, composite identity keys
are constructed by concatenating elements from each column with
"_" added as a separator.

The frequency of missing values (against /code{na.freq}) is tested first
by column then by row.

Selection can take three values: "" for no selection, "union" for all
identities expanded to every dataset, "shared" for only those data points
present in all usable datasets or "distinct" for excluding any points that
can be found in more than one dataset. Note that the union may result in
rows with no usable values.
}

\value{
A data frame if only one input dataset, or a list of data frames if
multiple datasets.
}

\author{
Ville-Petteri Makinen
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Create new versions for testing.
dsA <- dataset[1:250, c("INDEX","AGE","MALE","uALB")]
dsB <- dataset[151:300, c("INDEX","AGE","MALE","uALB","CHOL")]
dsC <- dataset[201:500, c("INDEX","AGE","MALE","DIAB_RETINO")]

# Select all rows.
results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX")
cat("\n\nNo selection:\n")
print(nrow(results$a))
print(nrow(results$b))
print(nrow(results$c))

# Select all rows and expanded for all identities.
results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                        select = "union")
cat("\n\nUnion:\n")
print(nrow(results$a))
print(nrow(results$b))
print(nrow(results$c))

# Select only rows that are shared between all datasets.
results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                        select = "intersection")
cat("\n\nIntersection:\n")
print(nrow(results$a))
print(nrow(results$b))
print(nrow(results$c))

# Select only rows with a unique INDEX ('dsB' has none).
results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                        select = "exclusion")
cat("\n\nExclusion:\n")
print(nrow(results$a))
print(nrow(results$b))
print(nrow(results$c))

# Add extra identification information.
dsA$GROUP <- "A"
dsB$GROUP <- "B"
dsC$GROUP <- "C"

# Select rows with a unique identifier.
results <- numero.clean(a = dsA, b = dsB, c = dsC,
                        identity = c("GROUP","INDEX"),
                        select = "exclusion")
cat("\n\nMulti-identities:\n")
print(nrow(results$a))
print(nrow(results$b))
print(nrow(results$c))
}