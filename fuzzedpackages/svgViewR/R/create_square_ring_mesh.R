create_square_ring_mesh <- function(center, axes, dims, width, seg, parallel = TRUE, curve.seg=1:2, 
	mar = NULL, rev.normals = FALSE){

	# Add 0 to dims if only of length 2 (no curvature)
	if(length(dims) == 2) dims <- c(dims, 0)

	if(parallel){

		## One set of parallel planes
		n_center <- center + 0.5*(dims[1] - width)*axes[1,]
		n_dims <- c(width, dims[2])
		mesh1 <- create_plane_mesh(corners=set_plane_corners(center=n_center, dims=n_dims, vecs=axes), seg=seg, rev.normals=rev.normals)

		n_center <- center - 0.5*(dims[1] - width)*axes[1,]
		n_dims <- c(width, dims[2])
		mesh2 <- create_plane_mesh(corners=set_plane_corners(center=n_center, dims=n_dims, vecs=axes), seg=seg, rev.normals=rev.normals)

		## Other set of parallel planes
		n_center <- center + 0.5*(dims[2] - width)*axes[2,]
		n_dims <- c(dims[1] - 2*width, width)
		mesh3 <- create_plane_mesh(corners=set_plane_corners(center=n_center, dims=n_dims, vecs=axes), seg=seg, rev.normals=rev.normals)

		n_center <- center - 0.5*(dims[2] - width)*axes[2,]
		n_dims <- c(dims[1] - 2*width, width)
		mesh4 <- create_plane_mesh(corners=set_plane_corners(center=n_center, dims=n_dims, vecs=axes), seg=seg, rev.normals=rev.normals)

		return_mesh <- addMeshes(mesh1, mesh2, mesh3, mesh4)

	}else{

		## One set of parallel planes
		n_center <- center + 0.5*dims[1]*axes[1,]
		n_dims <- c(width, dims[2:3])
		mesh1 <- create_curved_plane_mesh(center=n_center, dims=n_dims, axes=axes[c(3,2,1),], 
			seg=seg, curve.seg=curve.seg, mar=mar, rev.normals=!rev.normals)

		n_center <- center - 0.5*dims[1]*axes[1,]
		n_dims <- c(width, dims[2:3])
		mesh2 <- create_curved_plane_mesh(center=n_center, dims=c(1,1,-1)*n_dims, axes=axes[c(3,2,1),], 
			seg=seg, curve.seg=curve.seg, mar=mar, rev.normals=rev.normals)

		## Other set of parallel planes
		n_center <- center + 0.5*dims[2]*axes[2,]
		n_dims <- c(dims[1], width, dims[3])
		mesh3 <- create_curved_plane_mesh(center=n_center, dims=n_dims, axes=axes[c(1,3,2),], 
			seg=seg, curve.seg=curve.seg, mar=rev(mar), rev.normals=!rev.normals)

		n_center <- center - 0.5*dims[2]*axes[2,]
		n_dims <- c(dims[1], width, dims[3])
		mesh4 <- create_curved_plane_mesh(center=n_center, dims=c(1,1,-1)*n_dims, axes=axes[c(1,3,2),], 
			seg=seg, curve.seg=curve.seg, mar=rev(mar), rev.normals=rev.normals)

		return_mesh <- addMeshes(mesh1, mesh2, mesh3, mesh4)
	}

	return_mesh
}