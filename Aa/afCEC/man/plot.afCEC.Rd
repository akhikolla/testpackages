\name{plot.afCEC}
\alias{plot.afCEC}

\title{
plot.afCEC
}

\description{
Plots the chart showing the clustering of the data. Depending on the dimensionality of the data passed to the afCEC function, it draws either 2D or 3D chart. The function doesn't work with the other dimensions.
\cr
\cr
\strong{Remarks:}
\itemize{
\item{
The particular items are drawn cluster-by-cluster in the following order:
\itemize{
\item{points}
\item{means}
\item{curves / surfaces}
\item{ellipses / ellipsoids} of confidence
}
}
\item{
The functionality concerning drawing the means (see the \code{draw_means} argument description for more details), ellipses / ellipsoids of confidence (see the \code{draw_ellipsoids} argument description for more details) and curves / surfaces (see the \code{draw_surfaces} argument description for more details) doesn't work if object x was retrieved after clustering the data with the active function family defined by the matrix of the active function components values on the data set. In that case, the above-mentioned arguments as well as the other ones related to them (for example \code{meansSize2D}) are ignored and only data points are drawn on the chart. See \link{afCEC}.
}
}
\strong{Example:} (due to the long running time it is not placed in the "Examples" section)
\cr
\cr
\code{library(afCEC);}\cr
\code{data(airplane);}\cr
\code{result <- afCEC(airplane, 17);}\cr
\code{plot(result);}
}

\usage{
\method{plot}{afCEC} (
    x,
    draw_points=TRUE,
    draw_means=TRUE,
    draw_ellipsoids=TRUE,
    draw_surfaces=TRUE,
    confidence=0.95,
    grid_resolution=32,

    pointsSize2D=1, pointsColor2D="cluster",
    meansSize2D=2, meansColor2D="black",
    ellipsesHeight2D=1, ellipsesColor2D="black",
    surfacesHeight2D=1, surfacesColor2D="black",
    XLabel2D="X", YLabel2D="Y",

    pointsAlpha3D=1, pointsSize3D=0.01, pointsColor3D="cluster",
    meansAlpha3D=0.5, meansSize3D=0.01, meansColor3D="black",
    ellipsoidsAlpha3D=0.25, ellipsoidsColor3D="cluster",
    surfacesAlpha3D=0.5, surfacesColor3D="cluster",
    XLabel3D="X", YLabel3D="Y", ZLabel3D="Z",
	
	...
)
}

\arguments{

\item{x}{
Object of class afCEC returned from the afCEC function.
}

\item{draw_points}{
If the value is TRUE, the function draws the points belonging to the particular clusters. The default value is TRUE.
}

\item{draw_means}{
If the value is TRUE, the function draws the means of the clusters as the big black dots / spheres. The default value is TRUE.
}

\item{draw_ellipsoids}{
If the value is TRUE, the function draws the curved ellipses / ellipsoids of confidence of the particular clusters. The default value is TRUE.
}

\item{draw_surfaces}{
If the value is TRUE, the function draws the curves / surfaces corresponding to the subspace spanned by the inactive axes of the particular clusters with respect to the curvilinear coordinate system determined by the active function with 0's in the means of the clusters. The default value is TRUE.
}

\item{confidence}{
Determines the percentile of data belonging to the particular clusters, the corresponding ellipses / ellipsoids of confidence should contain. For example, if the value 0.5 is specified, then, for each cluster, it's ellipse / ellipsoid of confidence will contain 0.5 * 100 [\%] of assigned points. The default value is 0.95.
}

\item{grid_resolution}{
Determines the grid resolution using to approximate the curves / surfaces drawn in the plot function. The default value is 32.
}

\item{pointsSize2D}{
Determines the size of the points (in pixels) drawn in the plot function. The argument is used only in the case of the 2D data. The default value is 1.
}

\item{pointsColor2D}{
Determines the color of the points drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for points belonging to the different clusters.}
}
The argument is used only in the case of the 2D data. The default value is "cluster".
}

\item{meansSize2D}{
Determines the size of the means (in pixels) drawn in the plot function. The argument is used only in the case of the 2D data. The default value is 2.
}

\item{meansColor2D}{
Determines the color of the means drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for means of the different clusters.}
}
The argument is used only in the case of the 2D data. The default value is "black".
}

\item{ellipsesHeight2D}{
Determines the thickness of the ellipses (in pixels) drawn in the plot function. The argument is used only in the case of the 2D data. The default value is 1.
}

\item{ellipsesColor2D}{
Determines the color of the ellipses drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for ellipses of confidence of the different clusters.}
}
The argument is used only in the case of the 2D data. The default value is "black".
}

\item{surfacesHeight2D}{
Determines the thickness of the curves (in pixels) drawn in the plot function (see the \code{draw_surfaces} argument description for more details). The argument is used only in the case of the 2D data. The default value is 1.
}

\item{surfacesColor2D}{
Determines the color of the curves drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for curves (see the \code{draw_surfaces} argument description for more details) belonging the different clusters.}
}
The argument is used only in the case of the 2D data. The default value is "black".
}

\item{XLabel2D}{
Determines the label on the X axis of the chart drawn by the plot function. The argument is used only in the case of the 2D data. The default value is "X".
}

\item{YLabel2D}{
Determines the label on the Y axis of the chart drawn by the plot function. The argument is used only in the case of the 2D data. The default value is "Y".
}

\item{pointsAlpha3D}{
Determines the opacity of the points drawn in the plot function. The argument is used only in the case of the 3D data. The default value is 1.
}

\item{pointsSize3D}{
Determines the size of the points drawn in the plot function. The size is expressed as the fraction of the maximal extent of the bounding box containing the whole of the data points. The argument is used only in the case of the 3D data. The default value is 0.01.
}

\item{pointsColor3D}{
Determines the color of the points drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for points belonging to the different clusters.}
}
The argument is used only in the case of the 3D data. The default value is "cluster".
}

\item{meansAlpha3D}{
Determines the opacity of the means drawn in the plot function. The argument is used only in the case of the 3D data. The default value is 0.5.
}

\item{meansSize3D}{
Determines the size of the means drawn in the plot function. The size is expressed as the fraction of the maximal extent of the bounding box containing the whole of the data points. The argument is used only in the case of the 3D data. The default value is 0.01.
}

\item{meansColor3D}{
Determines the color of the means drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for points belonging to the different clusters.}
}
The argument is used only in the case of the 3D data. The default value is "black".
}

\item{ellipsoidsAlpha3D}{
Determines the opacity of the ellipsoids drawn in the plot function. The argument is used only in the case of the 3D data. The default value is 0.25.
}

\item{ellipsoidsColor3D}{
Determines the color of the ellipsoids drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for ellipses of confidence of the different clusters.}
}
The argument is used only in the case of the 3D data. The default value is "cluster".
}

\item{surfacesAlpha3D}{
Determines the opacity of the surfaces drawn in the plot function. The argument is used only in the case of the 3D data. The default value is 0.5.
}

\item{surfacesColor3D}{
Determines the color of the surfaces drawn in the plot function. It can hold the following values:
\itemize{
\item{Any color format acceptable by the standard R plot function. In that case the color is fixed across the particular clusters.}
\item{The "cluster" string indicating, that the plot function will use unique colors for ellipses of confidence of the different clusters.}
}
The argument is used only in the case of the 3D data. The default value is "cluster".
}


\item{XLabel3D}{
Determines the label on the X axis of the chart drawn by the plot function. The argument is used only in the case of the 3D data. The default value is "X".
}

\item{YLabel3D}{
Determines the label on the Y axis of the chart drawn by the plot function. The argument is used only in the case of the 3D data. The default value is "Y".
}

\item{ZLabel3D}{
Determines the label on the Y axis of the chart drawn by the plot function. The argument is used only in the case of the 3D data. The default value is "Z".
}

\item{...}{
The remaining unused arguments of the built-in R plot function.
}

}

\value{
None.
}

\seealso{
\code{\link{afCEC}}
}

\examples{
	library(afCEC);
	data(dog);
	res<-afCEC(dog,5,numberOfStarts=1);
	plot(res);
}