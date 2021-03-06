#' Construct linkage groups using 2-point recombination fraction estimates
#'
#' Use two-point recombination fraction estimates to group markers into the specified number of linkage group, The input \code{initial} can be used to select groups of markers that will be assigned to the same group. Grouping is performed using hierarchical clustering (\code{hclust}) using either average, complete or single linkage. The matrix used for clustering can be either theta (the recombination fraction matrix), lod (the log likelihood ratio) or a combination of the two. 
#' @importFrom stats na.omit
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom stats cutree
#' @export
#' @param mpcross Object of class \code{mpcross}
#' @param groups The number of groups to be formed
#' @param initial A list, with each entry containing markers which are to be assigned to the same group. Markers can be referenced by name or by index
#' @param clusterBy The type of data to cluster by. Can be one of "theta" for recombination fraction, "lod" for log likelihood ration, or "combined" for a combination of the two
#' @param method The clustering method to use. Must be one of "single", "complete" or "average"
#' @return A copy of the input mpcross object, with an additional "lg" entry containing the groupings of the markers. In addition the recombination fraction estimates and genetic data are reordered according to the created groupings of markers. 
#' \item{lg$groups}{ Numeric vector giving the group to which each marker belongs}
#' \item{lg$all.groups}{ Numeric vector giving the numbers for any groups that are present}
#' @examples
#' map <- qtl::sim.map(len=rep(100, 2), n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, 
#'		qtl=matrix(data=c(1, 50, .4, 0, 0, 0), 
#'		nrow=1, ncol=6, byrow=TRUE), seed=1)
#' dat.rf <- mpestrf(sim.dat)
#' grouped <- mpgroup(dat.rf, groups=2, clusterBy="combined", method="average")
#' grouped$lg
mpgroup <- function(mpcross, groups, initial = NULL, clusterBy="combined", method="average")
{
	if(!(clusterBy %in% c("combined", "theta", "lod")))
	{
		stop("Input clusterBy must be one of 'combined', 'theta' or 'lod'")
	}
	if(!(method %in% c("average", "complete", "single")))
	{
		stop("Input method must be one of 'average', 'complete' or 'single'")
	}
	if (missing(mpcross)) 
	{
		stop("Input mpcross cannot be missing")
	}

	if (is.null(mpcross$rf)&is.null(mpcross$map))
	{
		stop("Must calculate recombination fractions prior to grouping loci")
	}

	if(!is.null(mpcross$map)) {
	  cat("Using map groupings for groups. Remove map object if you want to regroup.\n")
	  initial <- lapply(mpcross$map, function(x) match(names(x), colnames(mpcross$finals)))
	  groups <- length(mpcross$map)
   	  output <- mpcross
    grpassignment <- vector(length=ncol(mpcross$finals))
    for (ii in 1:groups) grpassignment[initial[[ii]]] <- ii
    output$lg <- list(all.groups=1:groups, groups=grpassignment)
    names(output$lg$groups) <- unlist(lapply(mpcross$map, names))
	  return(output)
    ## don't rearrange order
   	}

	lod <- mpcross$rf$lod
	theta <- mpcross$rf$theta
	
	#Reverse lod so that small values indicate similarity
	lod[is.na(mpcross$rf$lod)] <- 0
	lod <- max(lod) - lod
	diag(lod) <- 0
	
	theta[is.na(mpcross$rf$theta)] <- 0.5

	if(method == "average")
	{
		linkFunction <- function(x) mean(x, na.rm=TRUE)
	}
	else if(method == "complete")
	{
		linkFunction <- function(x) max(x, na.rm=TRUE)
	}
	else
	{
		linkFunction <- function(x) min(x, na.rm=TRUE)
	}
	if(clusterBy == "combined")
	{
		distMat <- theta + lod / max(lod) * min(abs(diff(mpcross$rf$r)))
	}
	else if(clusterBy == "theta")
	{
		distMat <- theta
	}
	else
	{
		distMat <- lod
	}
	if(is.null(initial))
	{
		clustered <- hclust(as.dist(distMat), method=method)
	}
	else
	{
		#The number of values we're going to shrink the distance matrix by, as a result of the pre-grouping
		nLess <- sum(unlist(lapply(initial, function(x) length(x) - 1)))
		if(nLess == 0)
		{
			clustered <- hclust(as.dist(distMat), method=method)
		}
		else
		{
			#We actually have groups input. These groups go first, followed by the ungrouped values
			ungrouped <- setdiff(colnames(mpcross$rf$theta), unlist(initial))
			new.dist <- matrix(data=0, nrow=ncol(mpcross$rf$theta) - nLess, ncol = ncol(mpcross$rf$theta) - nLess)
			#Are there even any ungrouped values? If there aren't, skip this bit
			onlyGrouped <- length(initial) == ncol(mpcross$rf$theta) - nLess
			if(!onlyGrouped)
			{
				new.dist[(length(initial)+1):(ncol(mpcross$rf$theta) - nLess), (length(initial)+1):(ncol(mpcross$rf$theta) - nLess)] <- distMat[ungrouped, ungrouped]
			}
			for(i in 1:length(initial))
			{
				markersI <- initial[[i]]
				for(j in 1:length(initial))
				{
					markersJ <- initial[[j]]
					new.dist[i, j] <- new.dist[j, i] <- linkFunction(distMat[markersI, markersJ])
				}
				if(!onlyGrouped)
				{
					for(j in 1:length(ungrouped))
					{
						new.dist[i, j+length(initial)] <- linkFunction(distMat[markersI, ungrouped[j]])
					}
				}
			}
			clustered <- hclust(as.dist(new.dist), members=c(unlist(lapply(initial, length)), rep(1, length(ungrouped))), method=method)
		}
	}
	cut <- cutree(clustered, k=groups)
	if(!is.null(initial))
	{
		specifiedGroups <- cut[1:length(initial)]
		cut <- cut[-(1:length(initial))]
		names(cut) <- ungrouped
		for(i in 1:length(initial))
		{
			new <- rep(specifiedGroups[i], length(initial[[i]]))
			names(new) <- initial[[i]]
			cut <- c(cut, new)
		}
		cut <- cut[colnames(mpcross$rf$theta)]
	}
	else
	{
		names(cut) <- colnames(mpcross$rf$theta)
	}
	
	output <- mpcross
	output$lg <- list(all.groups=1:groups, groups=cut)
	return(subset(output, markers = colnames(output$founders)[order(output$lg$groups)]))
}
mpsplitgroup <- function(mpcross, toSplit, nSplits, clusterBy="combined", method="average")
{
	if(!(method %in% c("average", "complete", "single")))
	{
		stop("Input method must be one of 'average', 'complete' or 'single'")
	}
	if(!(clusterBy %in% c("combined", "theta", "lod")))
	{
		stop("Input clusterBy must be one of 'combined', 'theta' or 'lod'")
	}
  
	if (missing(mpcross)) 
	{
		stop("Input mpcross cannot be missing")
	}

	if (is.null(mpcross$rf))
	{
		stop("Must calculate recombination fractions prior to grouping loci")
	}
	if(is.null(mpcross$lg))
	{
		stop("Must have an existing grouping structure to call mpsubgroup")
	}
	lod <- mpcross$rf$lod
	theta <- mpcross$rf$theta
	
	#Reverse lod so that small values indicate similarity
	lod[is.na(mpcross$rf$lod)] <- 0
	lod <- max(lod) - lod
	diag(lod) <- 0
	
	theta[is.na(mpcross$rf$theta)] <- 0.5

	if(clusterBy == "combined")
	{
		distMat <- theta + lod / max(lod) * min(abs(diff(mpcross$rf$r)))
	}
	else if(clusterBy == "theta")
	{
		distMat <- theta
	}
	else
	{
		distMat <- lod
	}
	new.groups <- mpcross$lg$groups
	
	current.group <- which(mpcross$lg$groups == toSplit)
	subdist <- as.dist(distMat[current.group, current.group])
	clustered <- hclust(subdist, method=method)
	cut <- cutree(clustered, k=nSplits)
	new.groups[new.groups > toSplit] <- new.groups[new.groups > toSplit] + nSplits-1
	new.groups[current.group] <- cut + toSplit - 1
	
	output <- mpcross
	all.groups <- mpcross$lg$all.groups
	all.groups[all.groups > toSplit] <- all.groups[all.groups > toSplit] + nSplits - 1
	all.groups <- unique(c(all.groups, toSplit:(toSplit + nSplits-1)))
	output$lg <- list(all.groups=all.groups, groups=new.groups)
	return(subset(output, markers = colnames(output$founders)[order(output$lg$groups)]))
}
mpsubgroup <- function(mpcross, subgroups, clusterBy="combined", method="average")
{
	if(!(method %in% c("average", "complete", "single")))
	{
		stop("Input method must be one of 'average', 'complete' or 'single'")
	}
	if(!(clusterBy %in% c("combined", "theta", "lod")))
	{
		stop("Input clusterBy must be one of 'combined', 'theta' or 'lod'")
	}
  	if (missing(mpcross)) 
	{
		stop("Input mpcross cannot be missing")
	}
	if (is.null(mpcross$rf))
	{
		stop("Must calculate recombination fractions prior to grouping loci")
	}
	if(is.null(mpcross$lg))
	{
		stop("Must have an existing grouping structure to call mpsubgroup")
	}
	
	lod <- mpcross$rf$lod
	theta <- mpcross$rf$theta
	
	#Reverse lod so that small values indicate similarity
	lod[is.na(mpcross$rf$lod)] <- 0
	lod <- max(lod) - lod
	diag(lod) <- 0
	
	theta[is.na(mpcross$rf$theta)] <- 0.5

	if(clusterBy == "combined")
	{
		distMat <- theta + lod / max(lod) * min(abs(diff(mpcross$rf$r)))
	}
	else if(clusterBy == "theta")
	{
		distMat <- theta
	}
	else
	{
		distMat <- lod
	}
	new.groups <- vector(mode="integer", length=length(mpcross$lg$groups))
	names(new.groups) <- names(mpcross$lg$groups)
	
	for(index in 1:length(mpcross$lg$all.groups))
	{
		group <- mpcross$lg$all.groups[index]
		current.group <- which(mpcross$lg$groups == group)
		subdist <- as.dist(distMat[current.group, current.group])
		clustered <- hclust(subdist, method=method)
		cut <- cutree(clustered, k=subgroups)
		new.groups[current.group] <- cut + ((index -1)*subgroups)
	}
	
	output <- mpcross
	output$lg <- list(n.groups=length(mpcross$lg$all.groups) * subgroups, groups=new.groups)
	return(subset(output, markers = colnames(output$founders)[order(output$lg$groups)]))
}
fineTuneGroups <- function(grouped, excludeGroups=c())
{
	if(is.null(names(grouped$lg$groups))) stop("Invalid mpcross object input")
	originalChromosome <- newChromosome <- markerName <- newAverage <- c()
	for(i in 1:length(grouped$lg$all.groups))
	{
		indicesI <- which(grouped$lg$groups == grouped$lg$all.groups[i])
		diagonal <- grouped$rf$theta[indicesI, indicesI, drop=FALSE]
		originalAverages <- apply(diagonal, 1, function(x) mean(x, na.rm=TRUE))
		for(j in grouped$lg$all.groups[-i])
		{
			indicesJ <- which(grouped$lg$groups == j)
			offDiagonal <- grouped$rf$theta[indicesI, indicesJ, drop=FALSE]
			averages <- apply(offDiagonal, 1, function(x) mean(x, na.rm=TRUE))
			
			#Ok, for these markers names chromosome j is better than chromosome i. 
			betterMarkerNames <- names(which(averages < originalAverages))
			
			#For these ones we have no previously better chromosome, so add them straight in
			firstSelection <- betterMarkerNames[!(betterMarkerNames %in% markerName)]
			#For these we've already made a selection. But is j also better than other any other previously looked choices. If it's mostly unlinked to its current chromosome then a lot of other bad chromosomes will probably be a little bit better too
			secondSelection <- betterMarkerNames[betterMarkerNames %in% markerName]
			
			markerName <- c(markerName, firstSelection)
			originalChromosome <- c(originalChromosome, rep(grouped$lg$all.groups[i], length(firstSelection)))
			newChromosome <- c(newChromosome, rep(j, length(firstSelection)))
			newAverage <- c(newAverage, averages[firstSelection])
			
			#Which of these ones to keep 
			keep <- newAverage[match(secondSelection, markerName)] > averages[secondSelection]
			#And the corresponding bits of the previous relocation data that has to be removed
			overwrite <- match(secondSelection[keep], markerName)
			secondSelection <- secondSelection[keep]
			
			newChromosome[overwrite] <- rep(j, length(secondSelection))
			newAverage[overwrite] <- averages[secondSelection]
		}
	}
	#Actually do the re-arranging
	for(index in 1:length(markerName))
	{
		currentNewChromosome <- newChromosome[index]
		if(!(currentNewChromosome %in% excludeGroups))
		{
			currentMarkerName <- markerName[index]
			original <- originalChromosome[index]
			grouped$lg$groups[currentMarkerName] <- currentNewChromosome
		}
	}
	remove <- c()
	for(i in grouped$lg$all.groups)
	{
		indicesI <- which(grouped$lg$groups == i)
		diagonal <- grouped$rf$theta[indicesI, indicesI, drop=FALSE]
		additionalRemove <- which(apply(diagonal, 1, function(x) mean(x, na.rm=TRUE)) > 0.41)
		remove <- c(remove, additionalRemove)
	}
	
	markers <- colnames(grouped$founders)[order(grouped$lg$groups)]
	markers <- markers[-match(names(remove), markers)]	

	newGrouped <- subset(grouped, markers = markers)
	return(newGrouped)
}
joinGroups <- function(mpcross, join)
{
	#The value contained at index i is the new group ID for that group
	new.group.id <- vector(mode="numeric", length=max(mpcross$lg$all.groups))
	new.group.id[] <- NA
	next.group.id <- 1
	for(i in 1:length(join))
	{
		group <- as.integer(join[[i]])
		#None of these have existing groups, so we're creating a new group
		if(all(is.na(new.group.id[group])))
		{
			new.group.id[group] <- next.group.id
			next.group.id <- next.group.id + 1
		}
		else
		{
			#At least one of these is already part of a group, in fact there might be two groups being joined here. 
			existing <- unique(na.omit(new.group.id[group]))
			#We're joining different groups
			if(length(existing) > 1)
			{
				new.group.id[unique(c(group, which(new.group.id %in% existing)))] <- min(existing)
			}
			else
			{
				new.group.id[group] <- existing
			}
		}
	}
	for(i in mpcross$lg$all.groups)
	{
		if(is.na(new.group.id[i]))
		{
			new.group.id[i] <- next.group.id
			next.group.id <- next.group.id + 1
		}
	}
	previous.names <- names(mpcross$lg$groups) 
	mpcross$lg$groups <- new.group.id[mpcross$lg$groups]
	names(mpcross$lg$groups) <- previous.names
	
	mpcross$lg$all.groups <- unique(mpcross$lg$groups)
	return(mpcross)
}
findJoinPoint <- function(mpcross, marker1, marker2, joins)
{
	group1 <- mpcross$lg$groups[marker1]
	group2 <- mpcross$lg$groups[marker2]
	for(i in 1:length(joins))
	{
		join <- joins[[i]]
		if(join[[1]] == "join")
		{
			joinCommand <- c(join[[2]], join[[3]])
			if(all(c(group1, group2) %in% joinCommand))
			{
				return(i)
			}
			else if(group1 %in% joinCommand)
			{
				group1 <- min(joinCommand)
			}
			else if(group2 %in% joinCommand)
			{
				group2 <- min(joinCommand)
			}
		}
	}
	return(-1)
}
