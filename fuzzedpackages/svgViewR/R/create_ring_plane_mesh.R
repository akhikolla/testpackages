create_ring_plane_mesh <- function(center, axis, outer.radius, inner.radius, rseg=10, hseg=2, 
	theta.start=0, theta.length=2*pi, rev.normals=FALSE){

	mesh <- create_cylinder_mesh(center=center, axis=axis, radius=c(outer.radius, inner.radius), 
		height=0, rseg=rseg, hseg=hseg, theta.start=theta.start, theta.length=theta.length, 
		rev.normals=rev.normals)

	mesh
}