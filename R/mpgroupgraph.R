mpgroupgraph <- function(mpcross, threshold=0.2)
{
	if (missing(mpcross)) 
	{
		stop("Must input mpcross object to create linkage groups")
	}

	if (is.null(mpcross$rf))
	{
		stop("Must calculate recombination fractions prior to grouping loci")
	}
	if(is.null(mpcross$lg))
	{
		stop("Must have an existing grouping structure to call mpsubgroup")
	}
	if(!(require(graph) && require(igraph) && require(RBGL)))
	{
		stop("Packages 'graph', 'igraph' and 'RBGL' are required to use mpgroupgraph")
	}

	n.groups <- length(mpcross$lg$all.groups)
	
	distances <- matrix(0, n.groups, n.groups)
	for(i in mpcross$lg$all.groups)
	{
		markerI <- names(which(mpcross$lg$groups==i))
		for(j in mpcross$lg$all.groups)
		{
			markerJ <- names(which(mpcross$lg$groups==j))
			distances[i, j] <- distances[j, i] <- mean(mpcross$rf$theta[markerI, markerJ], na.rm=TRUE)
		}
	}
	diag(distances) <- 0.5

	zeros <- (distances > threshold) | is.na(distances)
	distances[zeros] <- 0
		
	width <- (0.5 - distances)*10
		
	adjacency <- matrix(1, n.groups, n.groups)
	adjacency[as.vector(zeros)] <- 0
	dim(adjacency) <- c(n.groups, n.groups)

	graph <- new("graphAM", adjMat=adjacency, edgemode="undirected")
	
	nodes(graph) <- as.character(mpcross$lg$all.groups)
	graph <- as(graph, "graphNEL")
	ig <- igraph.from.graphNEL(graph)
	E(ig)$color <- "black"

	cliques <- maxClique(graph)$maxCliques
	return(list(igraph = ig, cliques = cliques))
}
