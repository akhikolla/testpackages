\name{auto}
\alias{auto}
\docType{data}
\title{Automobile data}
\description{
This data set is taken from UCI repository, see reference. Past usage includes price prediction of cars using 
all numeric and boolean attributes (Kibler et al., 1989).
}
\usage{data(auto)}
\format{
  A data frame with 205 observations on the following 26 variables, of which 15 are quantitative and 11 are categorical. The following description is extracted
from UCI repository (Frank and Asuncion, 2010):  

	\tabular{ll}{
		\code{Normalized-losses} \tab the relative average loss payment per insured vehicle year; ranged from 65 to 256 \cr
		\code{Make} \tab Vehicle's make\cr
		\code{Fuel-type} \tab diesel, gas\cr
		\code{Aspiration} \tab std, turbo\cr
		\code{Num-of-doors} \tab four, two\cr
		\code{Body-style} \tab hardtop, wagon, sedan, hatchback, convertible\cr
		\code{Drive-wheels} \tab 4wd, fwd, rwd\cr
		\code{Engine-location} \tab front, rear\cr
		\code{Wheel-base} \tab continuous from 86.6 120.9\cr
		\code{Length} \tab continuous from 141.1 to 208.1\cr
		\code{Width} \tab continuous from 60.3 to 72.3\cr
		\code{Height} \tab continuous from 47.8 to 59.8\cr
		\code{Curb-weight} \tab continuous from 1488 to 4066\cr
		\code{Engine-type} \tab dohc, dohcv, l, ohc, ohcf, ohcv, rotor\cr
		\code{Num-of-cylinders} \tab eight, five, four, six, three, twelve, two\cr
		\code{Engine-size} \tab continuous from 61 to 326\cr
		\code{Fuel-system} \tab 1bbl, 2bbl, 4bbl, idi, mfi, mpfi, spdi, spfi\cr
		\code{Bore} \tab continuous from 2.54 to 3.94\cr
		\code{Stroke} \tab continuous from 2.07 to 4.17\cr
		\code{Compression-ratio} \tab continuous from 7 to 23\cr
		\code{Horsepower} \tab continuous from 48 to 288\cr
		\code{Peak-rpm} \tab continuous from 4150 to 6600\cr
		\code{City-mpg} \tab continuous from 13 to 49\cr
		\code{Highway-mpg} \tab continuous from 16 to 54\cr
		\code{Price} \tab continuous from 5118 to 45400\cr
		\code{Symboling} \tab assigned insurance risk rating: -3, -2, -1, 0, 1, 2, 3\cr
  }
}
\references{
Frank, A. & Asuncion, A. (2010). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. 
Irvine, CA: University of California, School of Information and Computer Science.

Kibler, D., Aha, D.W., & Albert,M. (1989). Instance-based prediction of real-valued attributes. Computational Intelligence, Vol 5, 51--57. 
}
\source{
The original data have been taken from the UCI Repository Of Machine Learning
Databases at
  \itemize{
    \item \url{http://archive.ics.uci.edu/ml/datasets/Automobile}.
  }  
 }
\keyword{datasets}


