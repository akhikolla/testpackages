\name{labkey.makeRemotePath}
\alias{labkey.makeRemotePath}
\title{Build a file path to data on a remote machine}
\description{ Replaces a local root with a remote root given a full path }

\usage{
labkey.makeRemotePath(localRoot, remoteRoot, fullPath) 
}
\arguments{
  \item{localRoot}{local root part of the fullPath}
  \item{remoteRoot}{remote root that will replace the local root of the fullPath}
  \item{fullPath}{the full path to make remote}
}

\details{  
A helper function to translate a file path on a LabKey web server to a path accessible 
by a remote machine.  For example, if an R script is run on an R server that is a different 
machine than the LabKey server and that script references data files on the LabKey server, 
a remote path needs to be created to correctly reference these files.  The local and 
remote roots of the data pipeline are included by LabKey in the prolog of an R View report 
script.  Note that the data pipeline root references are only included if an administrator 
has enabled the Rserve Reports experimental feature on the LabKey server.  If the 
remoteRoot is empty or the fullPath does not contain the localRoot then the fullPath is 
returned without its root being changed.
}

\value{ A character array containing the full path. 
}
\author{Dax Hawkins}
\examples{

# library(Rlabkey)

labkey.pipeline.root <- "c:/data/fcs"
labkey.remote.pipeline.root <- "/volumes/fcs"
fcsFile <- "c:/data/fcs/runA/aaa.fcs"

# returns "/volumes/fcs/runA/aaa.fcs
labkey.makeRemotePath(
	localRoot=labkey.pipeline.root,
	remoteRoot=labkey.remote.pipeline.root,
	fullPath=fcsFile);

}
\keyword{file}
