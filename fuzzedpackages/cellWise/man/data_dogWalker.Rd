\name{data_dogWalker}
\alias{data_dogWalker}
\docType{data}
\title{
  Dog walker dataset
}
\description{
  A dataset containing the image sequence of a video. The sequence consists of 
  54 frames of 144 by 180 pixels pixels in Red/Geen/Blue (RGB) format.
}
\usage{data("data_dogWalker")}

\format{
  An array of dimensions \eqn{54 \times 144 \times 180 \times 3}.
}

\source{
 \url{http://www.wisdom.weizmann.ac.il/~vision/SpaceTimeActions.html}
}

\examples{
data("data_dogWalker")
# For more examples, we refer to the vignette:
vignette("Wrap_examples")
}
