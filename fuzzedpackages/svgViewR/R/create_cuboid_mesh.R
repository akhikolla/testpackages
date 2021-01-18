create_cuboid_mesh <- function(ends, width, axes, rev.faces = FALSE){

	# Create vertices matrix
	vertices <- matrix(NA, 8, 3)

	# Add vertices
	vertices[1,] <- ends[1,] + (width[1]/2)*axes[2,] + (width[2]/2)*axes[3,]
	vertices[2,] <- ends[1,] + (width[1]/2)*axes[2,] - (width[2]/2)*axes[3,]
	vertices[3,] <- ends[1,] - (width[1]/2)*axes[2,] + (width[2]/2)*axes[3,]
	vertices[4,] <- ends[1,] - (width[1]/2)*axes[2,] - (width[2]/2)*axes[3,]
	vertices[5,] <- ends[2,] + (width[1]/2)*axes[2,] + (width[2]/2)*axes[3,]
	vertices[6,] <- ends[2,] + (width[1]/2)*axes[2,] - (width[2]/2)*axes[3,]
	vertices[7,] <- ends[2,] - (width[1]/2)*axes[2,] + (width[2]/2)*axes[3,]
	vertices[8,] <- ends[2,] - (width[1]/2)*axes[2,] - (width[2]/2)*axes[3,]

	# Create faces matrix
	faces <- matrix(NA, 12, 3)
	
	# Add faces
	faces[1,] <- c(2,1,0)
	faces[2,] <- c(2,3,1)
	faces[3,] <- c(4,5,6)
	faces[4,] <- c(5,7,6)
	faces[5,] <- c(6,3,2)
	faces[6,] <- c(6,7,3)
	faces[7,] <- c(1,3,7)
	faces[8,] <- c(1,7,5)
	faces[9,] <- c(6,2,4)
	faces[10,] <- c(4,2,0)
	faces[11,] <- c(1,5,4)
	faces[12,] <- c(1,4,0)
	
	if(rev.faces) faces <- faces[, c(3,2,1)]
	
	# Create normals
	normals <- matrix(NA, 8, 3)
	normals[1,] <- c(1,1,1)
	normals[2,] <- -axes[3,]
	normals[3,] <- axes[3,]
	normals[4,] <- axes[3,]
	normals[5,] <- -axes[1,]
	normals[6,] <- -axes[1,]
	normals[7,] <- -axes[1,]
	normals[8,] <- -axes[1,]
	
	# Normals
	center <- (ends[1,]+ends[2,]) / 2
	for(i in 1:8) normals[i,] <- c(vertices[i,] - center)
	normals <- uvector_svg(normals)

	list('vertices'=vertices, 'faces'=faces, 'normals'=normals)
}
