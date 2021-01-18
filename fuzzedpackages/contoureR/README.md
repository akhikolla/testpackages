[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/contoureR)](http://cran.r-project.org/web/packages/contoureR)
[![Downloads](http://cranlogs.r-pkg.org/badges/contoureR)](http://cran.rstudio.com/package=contoureR)

contoureR
======
Create contour lines for a non regular series of points, potentially from a non-regular canvas.

The contoureR package executes linear interpolation on a delaunay triangulated mesh strung between three-dimensional (3D) points supplied by the user. Contours are calculated across the surface constrained by the convex hull of the supplied data.

Usually, the well known functions such as contourLines from the grDevices package, expect (or rather, require) data to be regular, this means that a rectangular array or matrix of x and y coordinate pairs, each with a corresponding z value is to be modelled – that is to say the cartesian product of a numeric vector of x values of length n, with a numeric vector of y values having length m, used to produce a set of (m x n) unique points that have been concurrently provided with exactly (m x n) z values.

By restricting values to the above format, this in turn limits the region of analysis to square/rectangular canvasses (ie plane defined by geometric and orthogonal vectors parallel to the x and y axes and range bound by the [xmin,xmax] and [ymin,ymax] in the above x and y input numeric vectors, respectively). This restriction, from time-to-time, can be very inconvenient, and is a primary objective and purpose for the creation of this package.

As suggested in the previous paragraph, the contoureR package, on the other hand, has no such orthogonality / regularity requirement and can therefore be applied over obscurely shaped regions such as triangles, circles, polygons and the like.

For further examples and documentation, please proceed to the [contoureR website](http://contourer.com).