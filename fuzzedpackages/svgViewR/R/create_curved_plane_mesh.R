create_curved_plane_mesh <- function(center, dims, axes, seg, curve.seg, rev.normals, mar = NULL){

	# mar: fraction between 0 and 1 indicating the width of the margin on either side if "ramping up" to total curve height, can be two different values for each dimension

	# Create plane
	plane <- create_plane_mesh(corners=set_plane_corners(center=center, dims=dims[1:2], vecs=axes), 
		seg=seg, rev.normals=rev.normals)

	# Make sure of length 2
	if(length(seg) == 1) seg <- rep(seg, 2)
	
	# Get number of vertices along each dimension
	n_vertices <- seg+1
	
	# Get total number of vertices
	t_vertices <- prod(n_vertices)
	
	# Get normalized distance along each dimension for each vertex
	d1_norm_dist <- seq(0, 1, length=n_vertices[2])
	d2_norm_dist <- seq(0, 1, length=n_vertices[1])

	# Set number of indices to make margin
	if(!is.null(mar)){
	
		# Check that length is 2
		if(length(mar) == 1) mar <- rep(mar, 2)

		mar_n <- round(mar[1]*n_vertices[2])
		if(mar_n == 0) mar_n <- 1

		d1_norm_dist <- rep(0.5, n_vertices[2])
		d1_norm_dist[1:mar_n] <- seq(0, 0.5, length=mar_n)
		if(mar_n > 1){ d1_norm_dist[(n_vertices[2]-mar_n+1):n_vertices[2]] <- seq(0.5, 1, length=mar_n) }else{ d1_norm_dist[n_vertices[2]] <- 1 }

		mar_n <- round(mar[2]*n_vertices[1])
		if(mar_n == 0) mar_n <- 1

		d2_norm_dist <- rep(0.5, n_vertices[1])
		d2_norm_dist[1:mar_n] <- seq(0, 0.5, length=mar_n)
		if(mar_n > 1){ d2_norm_dist[(n_vertices[1]-mar_n+1):n_vertices[1]] <- seq(0.5, 1, length=mar_n) }else{ d2_norm_dist[n_vertices[1]] <- 1 }
	}

	# Ensure unit normal vector
	u_nvector <- uvector_svg(axes[3,])

	# Transform each vertex
	for(vn in 1:t_vertices){

		# Get row and column numbers
		d1 <- (vn-1) %% (seg[2]+1) + 1
		d2 <- (vn-1) %/% (seg[2]+1) + 1
		
		# Get normalized position along each dimension
		norm_xy <- c(d2_norm_dist[d2], d1_norm_dist[d1])

		## Get z-value (distance to raise from plane)
		# Curve along both axes
		if(length(curve.seg) == 2) z_val <- sign(dims[3])*prod(sqrt(abs(dims[3])) * (0.5 * sin(2*pi*norm_xy[1:2] - pi/2) + 0.5))
		if(length(curve.seg) == 1) z_val <- dims[3] * (0.5 * sin(2*pi*norm_xy[curve.seg] - pi/2) + 0.5)
		
		# Add vector to vertex position
		plane$vertices[vn, ] <- plane$vertices[vn, ] + z_val*u_nvector
	}

	plane	
}