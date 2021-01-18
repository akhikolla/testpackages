\name{calculate_p_value}
\alias{calculate_p_value}
\title{
Calculate asymptotic p-values
}
\description{
Calculate p-values given first order moment mu and second order moment std
}
\usage{
calculate_p_value(mu, std, val)
}
\arguments{
  \item{mu}{
  First order moment: mean for the GHD statistic
}
  \item{std}{
  Second order moment: std for the GHD statistic
}
  \item{val}{
  Value of GHD statistic
}
}
\value{
Returns the p-value for the GHD statistic (whose distribution follows normal distribution). 
}
\author{
Raghvendra Mall <rmall@hbku.edu.qa>
}