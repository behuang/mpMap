#' Order markers within linkage groups
#' 
#' Orders markers within linkage groups using two-point or multipoint probabilities. Two-point ordering is based on estimated recombination fractions; multi-point ordering is based on R/qtl ripple function. 
#' @export
#' @param object Object of class \code{mpcross}
#' @param chr Selected chromosomes or linkage groups to order
#' @param type Which type of ordering to perform - two-point or multipoint
#' @param window Window size for multipoint ordering
#' @param repeats Number of times to repeat multipoint ordering
#' @param mapfx Map function to use to compute final cM positions
#' @param criterion Criterion used in 2-pt ordering to determine best order
#' @param use.identity Options to improve 2-pt ordering via seriation 
#' @param seriate.control Options to improve 2-pt ordering via seriation
#' @param ... Additional arguments
#' @details \emph{Two-point ordering}\cr
#' To use the two-point ordering, the recombination fractions between all pairs of markers must first be estimated. If there are missing values in this matrix, the markers with the largest number of missing values will be removed until there are no missing values left. These markers will not be used in the ordering and are recommended to be inserted into the resulting framework map using \code{\link[mpMap]{add3pt}} later. 
#'
#' Multiple methods are used to investigate optimal two-point orderings. These are taken from the package \code{seriation} and include simulated annealing, hierarchical clustering, and traveling salesman solver. The orders are compared on the basis of the argument \code{criterion}. Thus the total path length, or sum of adjacent recombination fractions can be minimized; or the number of Anti-Robinson events/deviations; or the number of crossovers; or the sum of the adjacent two-point LOD scores. 
#'
#' \emph{Multi-point ordering} \cr
#' The multi-point ordering assumes that there is a pre-existing map, and then repeatedly applies the ripple function in R/qtl to investigate local permutations of the order. These orderings are constrained by the arguments \code{window} and \code{repeats}, which determine how large the perturbations are and how many are considered. Large values of \code{window} are very time consuming; recommended values are 5 or less, due to the number of permutations which must be considered. Large values of \code{repeats} will eventually converge to an ordering in which all local rearrangements of size \code{window} have been optimized with respect to the number of crossovers. 
#' @return The original object with a new map component. 
#' @seealso \code{\link[mpMap]{mpestrf}}, \code{\link[mpMap]{mpgroup}}, \code{\link[mpMap]{add3pt}}, \code{\link[seriation]{seriate}}, \code{\link[qtl]{ripple}}


mporder <-
function(object, chr, type=c("2", "m"), mapfx=c("haldane", "kosambi"), window=3, repeats=1, criterion=c("Path_length", "AR_events", "AR_deviations", "Gradient_raw", "Inertia", "Least_squares", "minXO"), use.identity = TRUE, seriate.control=NULL, ...)
{
	if (missing(object)) 
	{
		stop("Input object cannot be missing")
	}
	if (!inherits(object, "mpcross")) 
	{
		stop("Input object must have class 'mpcross'")
	}
	if (is.null(object$rf))
	{
		stop("Must calculate recombination fractions prior to ordering")
	}
	if(type != "2" && type != "m")
	{
		stop("Input type must have value '2' or 'm'")
	}

	if (missing(criterion)) criterion <- "Path_length"
	if (missing(mapfx)) 	mapfx <- "haldane" 
	
	decreasing <- FALSE
	if (criterion %in% c("Path_length", "AR_events", "AR_deviation", "Least_squares", "minXO")) decreasing <- TRUE

	if (mapfx=="haldane") mf <- haldaneR2X 
	else mf <- kosambiR2X

	#If we have a map use that, otherwise create dummy map using groups
	if (is.null(object$map) & is.null(object$lg))
	{
		stop("No grouping of markers input")
	}
	if (is.null(object$map) & !is.null(object$lg))
	{
		object$map <- vector(mode="list",length=length(object$lg$all.groups))
		for (i in 1:length(object$lg$all.groups))
		{
			object$map[[i]] <- rep(0, sum(object$lg$groups==object$lg$all.groups[i], na.rm=TRUE))
			names(object$map[[i]]) <- colnames(object$finals)[which(object$lg$groups==object$lg$all.groups[i])]
		}
	}
	output <- object
	
	if(any(is.na(object$rf$theta)))
	{
		object <- mpimputerf(object)
	}

	if (missing(chr))	chr <- c(1:length(object$map))  
	
	#If chr is a list of character names, convert them to the index of that chromosome within the map object
	if (is.character(chr)) chr <- match(chr, names(object$map))
	# do 2-pt ordering
	if (type=="2")
	{
		order <- list()
		if (criterion=="minXO") 
		{
			write2cross(object, "tmp", chr=chr)
			cr <- qtl:::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(object, "type"))
		}
		
		#Create a copy of the old map, purely to preserve the chromosome names
		newmap <- list()
		for (i in chr)
		{
			cat(paste("Ordering chromosome ", i, "...\n", sep=""))
			#Get the indices of the columns of theta for the markers that are in group / chromosome i
			nam <- match(names(object$map[[i]]), colnames(object$rf$theta))
			mat <- originalmat <- object$rf$theta[nam, nam, drop=FALSE] 
			diag(mat) <- diag(originalmat) <- 0
			#For the subgroups, replace that chunk with a single averaged column, and accept that afterwards we may get something that's flipped relative to what we want - fix it up later, in a later stage. 
			if("subgroups" %in% names(object$lg))
			{
				for(subgroup in 1:10)
				{
					markers <- intersect(names(which(object$lg$subgroups == subgroup)), names(object$map[[i]]))
					if(length(markers) > 0)
					{
						indices <- match(markers, colnames(mat))
						relevantColumns <- mat[,indices][-indices,]
						mat <- mat[-indices, -indices]
						replacementColumn <- apply(relevantColumns, 1, function(x) mean(x, na.rm=TRUE))
						
						mat <- rbind(cbind(mat, replacementColumn), c(replacementColumn, 0))
						newColName <- paste("subgroup", subgroup, sep="")
						if((newColName %in% rownames(mat)) || (newColName %in% colnames(mat))) stop("Markers with names beginning with 'subgroup' are reserved for internal use")
						rownames(mat)[nrow(mat)] <- colnames(mat)[ncol(mat)] <- newColName
					}
				}
			}
			dmat <- as.dist(mat)
		
			#If there are more than two markers, actually look for orderings
			if (length(object$map[[i]])>2)
			{
				order <- mporderchunk(object, dmat, criterion, cr, decreasing, use.identity, chromosome = i, seriate.control=seriate.control, ...)
				if(!is.null(dim(order))) order <- as.vector(order)
			} 
			#....otherwise there's only one ordering
			else 
			{
				order <- 1:length(object$map[[i]])
			}
			markerOrder <- rownames(mat)[order]
			#OK, now for the subgroups we need to order these separately
			if("subgroups" %in% names(object$lg))
			{
				for(subgroup in 1:10)
				{
					markers <- intersect(names(which(object$lg$subgroups == subgroup)), names(object$map[[i]]))
					if(length(markers) > 0)
					{
						submat <- originalmat[markers, markers]
						dsubmat <- as.dist(submat)
						suborder <- mporderchunk(object, dsubmat, criterion, cr, decreasing, use.identity, seriate.control=seriate.control, ...)
						if(!is.null(dim(suborder))) suborder <- as.vector(suborder)
					
						subgroupIndex <- match(paste("subgroup", subgroup, sep=""), markerOrder)
						if(subgroupIndex == 1)
						{
							markerOrder <- c(colnames(submat)[suborder], markerOrder[2:length(markerOrder)])
						}
						else if(subgroupIndex == length(markerOrder))
						{
							markerOrder <- c(markerOrder[1:(length(markerOrder)-1)], colnames(submat)[suborder])
						}
						else
						{
							markerOrder <- c(markerOrder[1:(subgroupIndex-1)], colnames(submat)[suborder], markerOrder[(subgroupIndex+1):length(markerOrder)])
						}
					}
				}
			}
			newmap[[i]] <- rep(0, length(markerOrder))
			names(newmap[[i]]) <- markerOrder
		}
		class(newmap) <- "map"
	}
	# do multi-point ordering
	else if (type=="m")
	{
		write2cross(object, "tmp", chr=chr)
		cr <- qtl:::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(object, "type"))
		newmap <- list()
		order <- list()
		chr <- match(names(object$map)[chr], names(cr$geno))
		for (i in chr)
		{
		rip <- ripple(cr, window=window, chr=names(cr$geno)[i])
		nmrk <- nmar(cr)[i]
		cat("Minimum XO for starting order: ", rip[1,nmrk+1], " for best order: ", rip[2,nmrk+1], "\n")
		cr2 <- cr
		order[[i]] <- rip[2, 1:nmrk]
		repeats2 <- repeats

			while(repeats2 >0 & (rip[1,nmrk+1]!=rip[2,nmrk+1]))
			{
		  ## actually have to go in and reorder the markers each time
		  cr2$geno[[i]]$data <- cr2$geno[[i]]$data[, order[[i]]]
		  rip <- ripple(cr2, window=window, chr=names(cr$geno)[i])	
		  cat("Minimum XO for starting order: ", rip[1,nmrk+1], " for best order: ", rip[2,nmrk+1],"\n")
		  order[[i]] <- rip[2, 1:nmrk]
		  repeats2 <- repeats2-1
			}
		## construct new map from ordering
		nam <- match(colnames(cr2$geno[[i]]$data), colnames(object$rf$theta))
		mat <- object$rf$theta[nam,nam]
		mat[mat==.5] <- .49
		newmap[[i]] <- cumsum(mf(c(0, mat[row(mat)==(col(mat)+1)][1:(length(nam)-1)])))
		names(newmap[[i]]) <- colnames(cr2$geno[[i]]$data)
		}
		names(newmap) <- names(cr$geno)[chr]
		class(newmap) <- "map"
	}
	output <- subset(output, markers = unlist(lapply(newmap, names)))
	#There's no guarantee any previously generated map will be correct any more, so remove it. 
	output$map <- newmap
	return(output)
}

mporderchunk <- function(object, dmat, criterion, cr, decreasing, use.identity, chromosome, seriate.control=NULL, ...)
{
	#Temporary workaround for a seriation package bug
	if(all(dmat == 0))
	{
		size <- attr(dmat, "Size")
		return(1:size)
	}
	## test out each of the different ordering techniques, pick the 
	## one with the shortest path length, according to our given loss function
	methods <- c("TSP", "OLO", "ARSA", "MDS", "GW", "HC")
	ser <- lapply(methods, function(x) return(seriate(dmat, method=x, control=seriate.control, ...))) 
	o2 <- do.call("rbind", lapply(ser, get_order))
	if(use.identity)
	{
		o2 <- rbind(1:ncol(o2), o2)	
	}
	if (criterion!="minXO") 
	{
		crit <- lapply(ser,function(x) return(criterion(dmat,x,criterion)))
		if(use.identity) crit <- c(criterion(dmat, method=criterion), crit)
		minx <- which.min(unlist(crit))
		if (!decreasing) minx <- which.max(unlist(crit))
	} 
	else 
	{
		crit <- compare_orders(cr, chr=chromosome, orders=o2, method="countxo")
		minx <- which.min(crit[,ncol(o2)+1]) 
	}
	if(use.identity && minx == 1) warning("Existing order was assessed as optimal, no reordering performed")
	return(o2[minx,1:ncol(o2), drop=FALSE])
}
