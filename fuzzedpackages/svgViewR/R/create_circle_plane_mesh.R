create_circle_plane_mesh <- function(center, nvector, radius, seg=10, rev.normals=FALSE){

	# Define circle
	circle <- defineCircle_svg(center=center, nvector=nvector, radius=radius)

	# Get vertices on edge of circle
	vertices_edge <- circlePoint_svg(circle, seq(0, 2*pi, length=seg+1)[1:seg])
	
	# Add center
	vertices <- rbind(center, vertices_edge)
	
	# Set normals
	normals <- matrix(nvector, nrow=nrow(vertices), ncol=3, byrow=TRUE)
	
	#
	faces <- matrix(NA, nrow=seg, ncol=3)
	faces[, 1] <- 1
	faces[, 2] <- 2:nrow(vertices)
	faces[1:(nrow(faces)-1),3] <- 3:nrow(vertices)
	faces[nrow(faces), 3] <- 2
	
	faces <- faces - 1
	
	if(rev.normals) faces <- faces[, c(1,3,2)]

	# Output vertices and faces
	list('vertices'=vertices, 'faces'=faces, 'normals'=normals, 'vertices_norm'=normals)
}