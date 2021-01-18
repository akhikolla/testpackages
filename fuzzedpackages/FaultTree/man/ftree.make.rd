\name{ftree.make}
\alias{ftree.make}

\title{ Fault Tree Creation }

\description{ftree.make prepares an initial, single row, dataframe, which will be built upon as the tree is developed.
}

\usage{
ftree.make(type, reversible_cond=FALSE, cond_first=TRUE, human_pbf=NULL,
		 start_id=1, system_mission_time=NULL, label="", name="", name2="",
		 description="")
}

\arguments{
\item{type}{ The logic type for the top gate. Possible gates are "or", "and", "inhibit", "priority", and "alarm". Gate types "vote" and "atleast" are disallowed in top event, place these under an OR gate. }
\item{reversible_cond}{A boolean value used only by the conditional gate type indicating whether repair of the input condition is viable to the model.}
\item{cond_first}{A boolean signifying whether the first child to an INHIBIT, ALARM, or PRIORITY logic gate should be taken as the condition.
 Alternatively, if FALSE, the second child will be taken as a condition.  This is primarily a graphic display issue.}
\item{human_pbf}{A probability of failure for a human to respond as needed to an alarm. This value is only used by the alarm gate.}
\item{start_id}{ An integer value for the starting unique ID, useful for transfer objects. }
\item{system_mission_time}{ An optional method for setting mission_time for the entire tree. It is also possible to set mission_time as a global variable. Use of system_mission_time in ftree.make will take precidence.}
\item{label}{An identifying string for the logic gate. Use of label defines the convention for rest of tree construction.}
\item{name}{ An identifying string for the logic gate.}
\item{name2}{ A second line, if needed for the identifying string label}
\item{description}{ An optional string providing more detail for the resultant event.}
}

\value{
Returns a dataframe holding data, results, and connection information.
}

\references{
  Nicholls, David [Editor] (2005) System Reliability Toolkit  Reliability information Analysis 
  Center
  
  O'Connor, Patrick D.T. (1991) Practical Reliability Engineering  John Wiley & Sons
  
  Vesely, W.E., Goldberg, F.F., Roberts, N.H., Haasl, D.F. (1981)  Fault Tree Handbook
  U.S.  Nuclear Regulatory Commission 
  
  Vesely, W.E., Stamatelato, M., Dugan, J., Fragola, J., Minarick, J., Railsback, J. (2002)
  Fault Tree Handbook with Aerospace Applications   NASA
  
  Doelp, L.C., Lee, G.K., Linney, R.E., Ormsby R.W. (1984) Quantitative fault tree analysis: Gate-by-gate method Plant/Operations Progress
  Volume 3, Issue 4 American Institute of Chemical Engineers
}

\examples{
mytree <- ftree.make(type="and", name="a specific undesired event")
}


