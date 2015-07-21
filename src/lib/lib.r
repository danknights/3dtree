# Given phylo object, finds n-D embedding
# that is faithful to branch lengths and minimizes
# stress from total tree-based distances
# topo.weight indicates how much to favor tree-based distance for non-connected nodes
# branch.weight indicates how much to favor branch lengths for connected nodes
# if coords.only, only returns embedded coords; else includes other data
"embed.phylo" <- function(tr, topo.weight=1, branch.weight=10, verbose=FALSE,coords.only=FALSE,
		optim.method=c('CG','BFGS','SANN')[1], maxit=100, reltol=1e-8){
	require('ape')
	ntips <- length(tr$tip.label)
	root.ix <- root.of(tr)
	ddt <- dist.nodes(tr)
	
	# function that returns stress of 3D coords vs. tree distance mat
	# takes a true distance matrix and an importance weighting matrix
	# stress is root mean squared error
	# note: coords are passed as a vector, as.numeric(coords)
	"tree.stress" <- function(coords, ddt, ddw){

		##### fill in here #######

	}
	
	# initial values based on pcoa
	pc <- cmdscale(ddt,k=3)
	if(verbose) cat('Optimizing...\n')

	# some parameters to control the optimization
	control.list <- list()
	if(verbose) control.list$trace <- 1
	control.list$maxit <- maxit
	control.list$reltol <- reltol

	# run the optimization on vector version of the coords
# 	res <- optim(as.numeric(pc),fn=tree.stress,gr=NULL,ddt=ddt,ddw=ddw,method=optim.method,control=control.list)

	# convert the final coords back to a matrix
# 	pc.hat <- matrix(res$par,nc=ncol(pc))
	pc.hat <- pc
	if(coords.only)	return(pc.hat)
	
	result <- list()
	result$coords <- pc.hat
	
	### uncomment when optimization is working
# 	result$tree.dist <- ddt
# 	result$weight.mat <- ddw
# 	result$stress <- res$value
# 	result$embedded.dist <- as.matrix(dist(pc.hat))
	return(result)
}

# plots a 3D interactive tree
# takes a tree and 3D coords
"plot.phylo.3D" <- function(tr,pc,color.by=NULL){
	require('rgl')
	require('RColorBrewer')
	cols <- c(brewer.pal(9,'Set1')[-6],brewer.pal(8,'Set2'),brewer.pal(9,'Set3'))
	cols <- sprintf('%saa',cols)
	if(length(unique(color.by)) > length(cols)) stop('Not enough colors available.')
	
	if(is.null(color.by)){
		cols <- rep('black',tr$Nnode + length(tr$tip.labels))
	} else {
		cols <- cols[1:length(unique(color.by))]
		cols <- cols[as.numeric(as.factor(color.by))]
		# if colors only provided for the tips, add black for internal nodes
		if(length(cols) < max(tr$edge[,1])) cols <- c(cols,rep('black',tr$Nnode))
	}
	
	sphere.sizes <- c(rep(.5,length(tr$tip.label)),rep(.25,tr$Nnode))
	sphere.sizes[root.of(tr)] <- 1

	open3d()
	par3d('FOV'=90,zoom=1,windowRect=c(0,0,600,600))
	decorate3d(range(pc[,1]), range(pc[,2]), range(pc[,3]),box=FALSE,axes=FALSE,aspect=TRUE,xlab='',ylab='',zlab='')
	plot3d(pc[,1], pc[,2], pc[,3],type='s',size=sphere.sizes,col=cols,add=TRUE)

	# draw lines
	for(i in 1:nrow(tr$edge)){
		a = tr$edge[i,1]
		b = tr$edge[i,2]
		plot3d(pc[c(a,b),1], pc[c(a,b),2], pc[c(a,b),3],lwd=.1,type='l',col='#00000055',add=TRUE)
	}
}

# returns index of root of tree
"root.of" <- function(tree) return(length(tree$tip.label) + 1)

# scales coords that are approaching the edge by shrinking their radial
# distance increasingly as they near the radius of the sphere
# center.type == 'root' means put the root node at origin
# center.type == 'range' means center the range of each axis on origin
# eps is the sharpness of the shrinkage, e.g. 10 is very sharp, 1 is no shrinkage
"spherify.coords" <- function(tr, pc, eps=4, radius=1, center.type=c('root','range')[1]){
	
	if(center.type == 'range'){
		# make each dim range from -1 to 1
		pc <- apply(pc,2, function(xx) (xx - min(xx))/abs(diff(range(xx))) * 2 - 1)
	} else {
		# or instead, make the root at (0,0,0)
		pc <- sweep(pc,2,pc[root.of(tr),],'-')
		# still shrink all coords so that max radius is 1
		pc <- pc / max(apply(pc,1,function(xx) sqrt(sum(xx^2))))
	}
	
	# shrink exponentially as radial distance approaches 1
	# doesn't make them "fan out"
	for(i in 1:nrow(pc)){
		di <- sqrt(sum(pc[i,]^2))
		if(di > 0){
			if(di > radius){
				dr <-  radius / di
			} else {
				dr <- 1-(1-di)^eps
				dr <- dr / di
			}
			pc[i,] <- pc[i,] * dr
		}
	}
	return(pc)
}
