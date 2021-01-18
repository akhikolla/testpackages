\name{SimHistory}
\alias{SimHistory}

\title{ Simulation History creation. }

\description{
This is a "Type 1" stochastic simulation engine.  It creates a simulation history dataframe for a single operating line.
}

\usage{
SimHistory(Model,SimulationYears=2000, SimulationYearsPerPage=1000, ProgRpt=FALSE) 
}

\arguments{
\item{Model}{ A datafreme constructed of one or more EventElement objects as combined by rbind.  All elements must have the same OpLine integer value in the first column.  It is expected that this dataframe will be constructed in a spreadsheet then transferred to R using the "Put R dataframe" right-click menu selection using RExcel}
\item{SimulationYears}{ A value for total size of simulation.  If not a multiple of SimulationYearsPerPage the actual simulation will be rounded down since number of pages is determined by Pages<-as.integer(SimulationYears/SimulationYearsPerPage)}
\item{SimulationYearsPerPage}{ A value for sub-setting the overall simulation; must be less then 2000 for accuracy.}
\item{ProgRpt}{ A boolean value indicating whether a progress bar should be displayed, if sensible, during execution
 of the function.}
}

\value{
Returns a dataframe containing columns for Time and Duration of each simulated event.  Added fields record the OpLine and EventID integer values pertaining to the event.
}

\references{
  Robert, Christian P., G. Casella (2010) Introducing Monte Carlo Methods with R.
  Springer
  
    Taylor HM,  Karlin S (1998) An Introduction to Stochastic Modeling, 3rd Edition,
  Acadmic Press.
}

\examples{
plantA_DF <- EventElement("generic.pump",1,101,"E", 28260,0,0,"N",8,2,0,87)
## note simulation drastically reduced for example run
PlantA <- SimHistory(plantA_DF,100,100)
}

\keyword{ model engine }

