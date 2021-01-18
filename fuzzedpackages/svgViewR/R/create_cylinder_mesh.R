create_cylinder_mesh <- function(center, axis, radius, height, rseg=10, hseg=2, 
	theta.start=0, theta.length=2*pi, rev.normals=FALSE){

	# Make sure vector is unit length
	axis <- uvector_svg(axis)
	
	# Make sure length 2
	if(length(radius) == 1) radius <- rep(radius, 2)

	# Set ends from center
	ends <- rbind(
		center + (height/2)*axis,
		center - (height/2)*axis
	)
	
	# Get vector orthogonal to axis
	o_axis <- vorthogonal_svg(axis)

	# Set thetas for end circle
	thetas <- seq(from=theta.start, to=theta.start+theta.length, length=rseg+1)[1:rseg]

	# Create matrix for end points
	ends_circle <- matrix(NA, nrow=rseg, ncol=3)

	# Get points on circumference of cylinder end circle
	for(i in 1:rseg) ends_circle[i, ] <- uvector_svg(o_axis %*% tMatrixEP_svg(axis, thetas[i]))

	# Draw vertices at lengths
	at_lengths <- seq(from=0, to=height, length=2+hseg-1)
	
	# Set radii
	radii <- radii <- seq(radius[1], radius[2], length=hseg+1)

	# Create vertices matrix
	n_vertices <- rseg*(hseg+1) + 2
	vertices <- matrix(NA, n_vertices, ncol=3)

	# Add vertices
	vertices[1, ] <- ends[1,]
	vertices[nrow(vertices), ] <- ends[2,]
	for(i in 1:length(at_lengths)){
		
		# Get rows to add
		at_rows <- 1 + (i-1)*rseg + (1:rseg)
		
		# Set vertices
		vertices[at_rows, ] <- radii[i]*ends_circle + matrix(ends[2,]+at_lengths[i]*axis, nrow=rseg, ncol=3, byrow=TRUE)
	}
	
	# Create side faces matrix
	faces <- matrix(NA, nrow=2*rseg*hseg, 3)

	# Fill side faces matrix
	for(i in 0:(hseg-1)){

		rows1 <- i*rseg*2 + (1:rseg)
		faces[rows1, 1] <- i*rseg + (1:rseg)
		faces[rows1, 2] <- faces[rows1, 1] + rseg
		faces[rows1, 3] <- faces[rows1, 2] + 1
		faces[tail(rows1,1), 3] <- faces[rows1[1], 1] + rseg

		rows2 <- rows1 + rseg
		faces[rows2, 1] <- i*rseg + (1:rseg)
		faces[rows2, 3] <- faces[rows2, 1] + 1
		faces[tail(rows2,1), 3] <- faces[rows2[1], 1]
		faces[rows2, 2] <- faces[rows1, 2] + 1
		faces[tail(rows2,1), 2] <- faces[rows1[1], 1] + rseg
	}
	
	# Remove ends
	vertices <- vertices[2:(nrow(vertices)-1),]

	# Shift down since end center is removed
	faces <- faces - 1
	
	#
	if(!rev.normals) faces <- faces[,c(1,3,2)]

	# Output vertices and faces
	list('vertices'=vertices, 'faces'=faces, 'normals'=NULL, 'vertices_norm'=NULL)
}