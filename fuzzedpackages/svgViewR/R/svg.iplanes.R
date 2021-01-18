svg.iplanes <- function(corners1, corners2 = NULL, center = FALSE, depth = NULL, col='blue', 
	emissive=rgb(0.03, 0.15, 0.21), opacity = 1, name = 'cplanes', seg = 30, ontop = FALSE, 
	create.uvs = FALSE, return.shape = FALSE, plot = TRUE){

	# Set number of faces
	n_faces <- 2*(nrow(corners1) - 2) + 2*nrow(corners1)
	n_vertices <- 2*nrow(corners1)
	n_vertices_h <- nrow(corners1)

	if(is.null(corners2)){

		# Create vertex normal matrix
		vnormals <- matrix(NA, n_vertices_h, 3)
	
		for(j in 1:(nrow(corners1)/2 - 1)){
	
			# Set vertices of plane
			v_c <- c(j, j+1, nrow(corners1)-j+1, nrow(corners1)-j)
			
			#
			vnormals[v_c[1],] <- uvector_svg(cprod_svg(corners1[v_c[2],]-corners1[v_c[1],], corners1[v_c[3],]-corners1[v_c[1],]))
			vnormals[v_c[3],] <- uvector_svg(cprod_svg(corners1[v_c[2],]-corners1[v_c[1],], corners1[v_c[3],]-corners1[v_c[1],]))
			
			if(j > 1){
				if(abs(sum(vnormals[v_c[1],])) == 0) vnormals[v_c[1],] <- prev_vnormal
				if(abs(sum(vnormals[v_c[3],])) == 0) vnormals[v_c[3],] <- prev_vnormal
			}
			
			if(j == nrow(corners1)/2 - 1){
				vnormals[v_c[2],] <- vnormals[v_c[1],]
				vnormals[v_c[4],] <- vnormals[v_c[3],]
			}
			
			prev_vnormal <- vnormals[v_c[1],]
		}
		
		# Standardize normals
		for(i in 2:nrow(vnormals)){
			if(abs(avec_svg(vnormals[i-1,], vnormals[i,])) > pi/2) vnormals[i,] <- -vnormals[i,]
		}
		
		# Check depth length
		if(length(depth) > nrow(vnormals)) stop(paste0("The length of 'depth' (", length(depth), ") cannot exceed the number of rows in 'vnormals' (", nrow(vnormals), ")"))

		# If center, project corners1 and corners2 half of depth using normals
		if(center){

			# Use vnormals to project corners2 if NULL
			corners2 <- corners1 + 0.5*depth*vnormals
			corners1 <- corners1 - 0.5*depth*vnormals
		}else{

			# Use vnormals to project corners2 if NULL
			corners2 <- corners1 + depth*vnormals
		}
	}

	# Create faces matrix
	faces <- matrix(NA, nrow=n_faces, ncol=3)
	
	# Create vertex list
	vertex_list <- list(corners1, corners2)
	vertex_mat <- rbind(corners1, corners2)

	## Create faces on "front" and "back"
	n <- 1
	for(i in 1:length(vertex_list)){
		for(j in 1:(nrow(vertex_list[[i]])/2 - 1)){
	
			# Set vertices of plane
			plane_vertices <- c(j, j+1, nrow(vertex_list[[i]])-j+1, nrow(vertex_list[[i]])-j)
		
			# Set faces
			faces[n,] <- plane_vertices[1:3] + (i-1)*n_vertices_h
			n <- n + 1
			faces[n,] <- plane_vertices[c(3,2,4)] + (i-1)*n_vertices_h
			n <- n + 1
		}
	}
	
	## Create faces "around the edges"
	for(j in 1:(nrow(vertex_list[[1]])-1)){

		# Set vertices of plane
		plane_vertices <- c(j, j+1, j + (i-1)*n_vertices_h, j + 1 + (i-1)*n_vertices_h)

		# Set faces
		faces[n,] <- plane_vertices[1:3]
		n <- n + 1
		faces[n,] <- plane_vertices[c(3,2,4)]
		n <- n + 1
	}

	# Add last face to close face
	plane_vertices <- c(1, n_vertices_h, n_vertices_h+1, 2*n_vertices_h)
	faces[nrow(faces)-1, ] <- plane_vertices[1:3]
	faces[nrow(faces), ] <- plane_vertices[c(3,2,4)]
	
	# Combine corners to set vertices
	vertices <- rbind(vertex_list[[1]], vertex_list[[2]])

	# Start face indices at 0
	faces <- faces - 1
	
	faces <- faces[!is.na(faces[,1]),]
	
	# If not plotting, return vertices and faces
	if(!plot) return(list('vertices'=vertices, 'faces'=faces))

	# Get viewer environment
	env <- as.environment(getOption("svgviewr_glo_env"))

	# Add to meshes
	add_at <- length(svgviewr_env$svg$mesh)+1

	# Add vertices
	svgviewr_env$svg$mesh[[add_at]] <- list()
	svgviewr_env$svg$mesh[[add_at]]$name <- name
	svgviewr_env$svg$mesh[[add_at]]$vertices <- t(vertices)
	svgviewr_env$svg$mesh[[add_at]]$faces <- t(faces)
	svgviewr_env$svg$mesh[[add_at]]$col <- setNames(webColor(col), NULL)
	svgviewr_env$svg$mesh[[add_at]]$opacity <- setNames(opacity, NULL)
	svgviewr_env$svg$mesh[[add_at]]$emissive <- setNames(webColor(emissive), NULL)
	svgviewr_env$svg$mesh[[add_at]]$computeVN <- TRUE
	svgviewr_env$svg$mesh[[add_at]]$parseModel <- FALSE
	svgviewr_env$svg$mesh[[add_at]]$depthTest <- !ontop
	svgviewr_env$svg$mesh[[add_at]]$material <- 'lambert'

	# Add object reference data
	svgviewr_env$ref$names <- c(svgviewr_env$ref$names, name)
	svgviewr_env$ref$num <- c(svgviewr_env$ref$num, add_at)
	svgviewr_env$ref$type <- c(svgviewr_env$ref$type, 'mesh')

	# Add limits
	obj_ranges <- apply(vertices, 2, 'range', na.rm=TRUE)
	
	# Set corners
	corners <- lim2corners(obj_ranges)
	
	# Add limits to object
	svgviewr_env$svg$mesh[[add_at]][['lim']] <- obj_ranges
	svgviewr_env$svg$mesh[[add_at]][['corners']] <- corners

	if(return.shape){
		return(list('vertices'=vertices, 'faces'=faces))
	}else{
		# Suppress return of value in console
		return(NULL)
	}
}
