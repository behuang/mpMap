outputPedigreeGraph <- function(pedigree, fileName, omitFunnel = TRUE)
{
	library(igraph)
	library(graph)
	if(missing(fileName))
	{
		stop("Input fileName cannot be missing")
	}
	if(missing(pedigree))
	{
		stop("Input pedigree cannot be missing")
	}
	pedigree <- convertped(pedigree)
	nFounders <- sum(pedigree[,2] == 0 & pedigree[,3] == 0)
	if(!(nFounders  %in% c(2, 4, 8)))
	{
		stop("Number of founders must be 2, 4 or 8")
	}
	if("Design" %in% colnames(pedigree))
	{
		pedigreeDesign <- pedigree[,"Design"]
	}
	else
	{
		pedigreeDesign <- identifyDesign(pedigree)
	}
	
	#Remove large chunks of the pedigree before we form the graph. This includes the 8wayG3 / 4wayG2 individuals. They'll be present in the graph because they're listed as mother / father though. 
	if(omitFunnel)
	{
		remove <- which(pedigreeDesign %in% c("8wayG1", "8wayG2", "8wayG0", "8wayG3", "4wayG0", "4wayG1", "4wayG2"))
		initialNodes <- pedigree[pedigreeDesign == "8wayG3" | pedigreeDesign == "4wayG2", 1]
	}
	else
	{
		remove <- which(pedigreeDesign %in% c("8wayG0", "4wayG0"))
		initialNodes <- pedigree[pedigreeDesign == "8wayG0" | pedigreeDesign == "4wayG0", 1]
	}
	cutDownPedigree <- pedigree[-remove,]

	edges <- list()
	f <- function(row) 
	{
		as.character(unique(c(cutDownPedigree[row, 2], cutDownPedigree[row, 3])))
	}
	edges <- lapply(1:nrow(cutDownPedigree), f)
	names(edges) <- as.character(cutDownPedigree[,1])
	
	otherEdges <- replicate(length(edges), character(0))
	names(otherEdges) <- names(edges)
	for(vertex in names(edges))
	{
		for(otherVertex in edges[[vertex]])
		{
			otherEdges[[otherVertex]] <- c(otherEdges[[otherVertex]], vertex)
		}
	}
	
	nodes <- c(cutDownPedigree[,1], initialNodes)
	subsettedDesign <- pedigreeDesign[match(nodes, pedigree[,1])]
	subsettedObserved <- pedigree[match(nodes, pedigree[,1]), "Observed"]
	
	if(is.numeric(pedigree[,1]) && !is.null(rownames(pedigree)))
	{
		edges <- lapply(edges, function(x) rownames(pedigree)[match(x, pedigree[,1])])
		names(edges) <- rownames(pedigree)[match(as.numeric(names(edges)), pedigree[,1])]
		nodes <- rownames(pedigree)[match(nodes, pedigree[,1])]
	}
	graph <- new("graphNEL", nodes=nodes, edgeL=edges, edgemode="directed")
	
	iPedGraph <- igraph.from.graphNEL(graph)
	iPedGraph <- set.vertex.attribute(iPedGraph, "design", value=subsettedDesign)
	iPedGraph <- set.vertex.attribute(iPedGraph, "isObserved", value = subsettedObserved)
	write.graph(iPedGraph, file=fileName, format = "dot")
}