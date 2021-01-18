afCEC <- function(points, maxClusters, initialLabels="k-means++", cardMin=0.01, costThreshold=-0.000001, minIterations=10, maxIterations=100, numberOfStarts=10, method="Hartigan", values="quadratic", interactive=FALSE) {
    if (class(values) == "character") {
        if (length(values) == 1) {
            if (values == "quadratic") {
                # Quadratic function in R^(d-1) of the form: x1^2+...+x(d-1)^2+x1+...+x(d-1)+1
                pointsNum = dim(points)[1];
                dimension = dim(points)[2];
                valuesArray = matrix(rep(0,((2*(dimension-1))+1)*dimension*pointsNum),((2*(dimension-1))+1)*dimension,pointsNum);
                for (i in 1:pointsNum) {
                    tmp = points[i,2:dimension];
                    for (j in 1:dimension) {
                        for (k in 1:(dimension - 1)) {
                            valuesArray[((j-1)*((2*(dimension-1))+1))+k,i]=tmp[k]^2;
                            valuesArray[((j-1)*((2*(dimension-1))+1))+dimension-1+k,i]=tmp[k];
                            valuesArray[((j-1)*((2*(dimension-1))+1))+(2*(dimension-1))+1,i]=1;
                        }
                        if (j < dimension) tmp[j] = points[i, j];
                    }
                }
            } else {
                stop("Wrong type of the \"active function\" parameter");
            }
        } else {
            stop("Wrong type of the \"active function\" parameter");
        }
    } else {
        if (is.matrix(values) && is.numeric(values)) {
            valuesArray = values;
        } else {
            stop("Wrong type of the \"active function\" parameter");
        }
    }
    res=afCECCppRoutine(t(points),maxClusters,initialLabels,cardMin,costThreshold,minIterations,maxIterations,numberOfStarts,method,valuesArray,interactive);
    if (length(res) > 0) {
        if (class(values) == "character") {
            if (!interactive) {
                UpdateMeansForQuadraticFunction(res);
                res[["data"]] = points;
                res[["formula"]] = values;
            } else {
                for (i in 1:length(res)) {
                    UpdateMeansForQuadraticFunction(res[[i]]);
                    res[[i]][["data"]] = points;
                    res[[i]][["formula"]] = values;
                }
            }
        } else {
            if (!interactive) {
                res[["data"]] = points;
                res[["formula"]] = NULL;
            } else {
                for (i in 1:length(res)) {
                    res[[i]][["data"]] = points;
                    res[[i]][["formula"]] = NULL;
                }
            }
        }

        if (!interactive) {
            res=structure(res,class="afCEC");
        } else {
            for (i in 1:length(res)) res[[i]] = structure(res[[i]], class="afCEC");
        }
        return(res);
    } else {
        stop("The afCEC function terminated with errors");
    }
}

#   -*-   -*-   -*-

plot.afCEC <- function(
    x, draw_points=TRUE, draw_means=TRUE, draw_ellipsoids=TRUE, draw_surfaces=TRUE, confidence=0.95, grid_resolution=32,

    pointsSize2D=1,pointsColor2D="cluster",
    meansSize2D=2,meansColor2D="black",
    ellipsesHeight2D=1,ellipsesColor2D="black",
    surfacesHeight2D=1,surfacesColor2D="black",
    XLabel2D="X",YLabel2D="Y",

    pointsAlpha3D=1,pointsSize3D=0.01,pointsColor3D="cluster",
    meansAlpha3D=0.5,meansSize3D=0.01,meansColor3D="black",
    ellipsoidsAlpha3D=0.25,ellipsoidsColor3D="cluster",
    surfacesAlpha3D=0.5,surfacesColor3D="cluster",
    XLabel3D="X",YLabel3D="Y",ZLabel3D="Z",

    ...
) {
    if (dim(x$data)[2] == 2) {
        if (draw_points) {
            if (pointsColor2D=="cluster") {
                plot(x$data,col=x$labels+1,asp=1,pch=20,cex=pointsSize2D,xlab=XLabel2D,ylab=YLabel2D);
            } else {
                plot(x$data,col=pointsColor2D,asp=1,pch=20,cex=pointsSize2D,xlab=XLabel2D,ylab=YLabel2D);
            }
        }
        if ((!is.null(x$formula)) && (draw_means || draw_ellipsoids || draw_surfaces)) {
            if (draw_ellipsoids || draw_surfaces) {
                if (x$formula == "quadratic") {
                    # quadratic function
                    ellipses=CalculateEllipsesOfConfidenceForQuadraticFunction(x, confidence, grid_resolution);
                }
            }
            for (i in 1:length(x$means)) {
                if (!is.null(x$means[[i]])) {
                    if (draw_means) points(t(x$means[[i]]),pch=20,cex=meansSize2D,col=ifelse(meansColor2D=="cluster",i,meansColor2D));
                    if (draw_surfaces) lines(ellipses[[2]][[i]],lwd=surfacesHeight2D,col=ifelse(surfacesColor2D=="cluster",i,surfacesColor2D));
                    if (draw_ellipsoids) lines(ellipses[[1]][[i]],lwd=ellipsesHeight2D,col=ifelse(ellipsesColor2D=="cluster",i,ellipsesColor2D));
                }
            }
        }
    } else if (dim(x$data)[2] == 3) {
        dataSize=dim(x$data)[1];
        dimX=max(x$data[1:dataSize,1])-min(x$data[1:dataSize,1]);
        dimY=max(x$data[1:dataSize,2])-min(x$data[1:dataSize,2]);
        dimZ=max(x$data[1:dataSize,3])-min(x$data[1:dataSize,3]);
        maxDim=max(dimX,dimY,dimZ);

        if (draw_points) {
            if (pointsColor3D=="cluster") {
                plot3d(x$data,col=x$labels+1,aspect=FALSE,alpha=pointsAlpha3D,size=maxDim*pointsSize3D,xlab=XLabel3D,ylab=YLabel3D,zlab=ZLabel3D);
            } else {
                plot3d(x$data,col=pointsColor3D,aspect=FALSE,alpha=pointsAlpha3D,size=maxDim*pointsSize3D,xlab=XLabel3D,ylab=YLabel3D,zlab=ZLabel3D);
            }
        }
        if ((!is.null(x$formula)) && (draw_means || draw_ellipsoids || draw_surfaces)) {
            if (draw_ellipsoids || draw_surfaces) {
                if (x$formula == "quadratic") {
                    # quadratic function
                    ellipsoids=CalculateEllipsoidsOfConfidenceForQuadraticFunction(x, confidence, grid_resolution);
                }
            }
            for (i in 1:length(x$means)) {
                if (!is.null(x$means[[i]])) {
                    if (draw_means) spheres3d(t(x$means[[i]]),radius=maxDim*meansSize3D,alpha=meansAlpha3D,col=ifelse(meansColor3D=="cluster",i,meansColor3D));
                    if (draw_surfaces) triangles3d(ellipsoids[[3]][[i]],normals=ellipsoids[[4]][[i]],alpha=surfacesAlpha3D,col=ifelse(surfacesColor3D=="cluster",i,surfacesColor3D));
                    if (draw_ellipsoids) triangles3d(ellipsoids[[1]][[i]],normals=ellipsoids[[2]][[i]],alpha=ellipsoidsAlpha3D,col=ifelse(ellipsoidsColor3D=="cluster",i,ellipsoidsColor3D));
                }
            }
        }
    } else {
        stop("Cannot plot data of dimensionality other than 2 or 3");
    }
}
