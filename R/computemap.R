#' Computes map distances
#' 
#' Given an mpcross object with a map order and matrix of recombination fractions, this function will estimate map positions, using the subdiagonal entries of the recombination fraction matrix and (optionally) other non-diagonal entries.
#' @export
#' @param object Object of class \code{mpcross}
#' @param mapfx Map function to convert recombination fractions to cM
#' @param maxOffset Controls which regions of the recombination fraction matrix are used to estimate distances. Specify 1 for super-diagonal elements, 2 to use one additional set of diagonal entries, etc
#'
#' @details
#' The distance between marker i and i+1 can be estimated as \code{mapfunction(object$rf$theta[i, i+1])}. Specifying maxOffset = 1 will use this method to estimate map distances. Specifying maxOffset = 2 will use \code{mapfunction(object$rf$theta[i, i+2])}, 
#' \code{mapfunction(object$rf$theta[i, i+1])} and \code{mapfunction(object$rf$theta[i+1, i+2])} to estimate the distances between markers i, i+1 and i+2 jointly. Similarly, larger values of maxOffset use recombination fractions between more distant markers in the estimation of map disatnce.
#' All estimation is done using non-negative linear least squares
#' @return An mpcross object is returned whose map component has been estimated based on the map order and matrix of recombination fractions. Missing recombination fractions are imputed either by filling in the closest non-missing value (missfx=1) or by averaging the distance between other nearby markers (missfx=2). 
#' @seealso \code{\link[mpMap]{mpcross}} \code{\link[nnls]{nnls}}


computemap <- function(object, mapfx=c("haldane", "kosambi"), maxOffset = 1)
{
	if (missing(mapfx)) mapfx <- "haldane"
	if (mapfx=="haldane") mf <- haldaneR2X else mf <- kosambiR2X

	if (is.null(object$rf))
	{
		warning("Recombination fractions have not been calculated. Recalculating....")
		object <- mpestrf(object)
	}

	if (is.null(object$map)) 
	{
		if("lg" %in% names(object))
		{
			map <- list()
			for(group in object$lg$all.groups)
			{
				markerNames <- names(which(object$lg$groups == group))
				chromosome <- rep(0, length(markerNames))
				names(chromosome) <- markerNames
				map[[length(map)+1]] <- chromosome
			}
			class(map) <- "map"
			object$map <- map
			names(object$map) <- paste("Chr", 1:length(map), sep="")
		}
		else
		{
			cat("Warning: no map order listed. Assuming all markers are on the same chromosome\n")
			newMap <- rep(0, length=nrow(object$rf$theta))
			names(newMap) <- colnames(object$rf$theta)
			object$map <- list(newMap)
			class(object$map) <- "map"
			names(object$map) <- "Chr1"
		}
	}
	#Remove lg component because we now have a map object, and mpimpute will warn if there are both.
	object$lg <- NULL
	object <- mpimputerf(object)
	for (chr in 1:length(object$map))
	{
		#Now rearrange recombination fractions into the same order as the design matrix
		m <- match(names(object$map[[chr]]), colnames(object$finals))
		rf <- object$rf$theta[m, m]
		maxOffset <- min(maxOffset,(nrow(rf)-1))
		#maxOffset <- (nrow(rf)-1)
		
		#Construct design matrix
		#d <- designMat(length(object$map[[chr]])-1, maxOffset)
		offset <- maxOffset
		if(offset > length(object$map[[chr]])-1) offset <- length(object$map[[chr]])-1
		d <- .Call("generateDesignMatrix", length(object$map[[chr]])-1, offset, package="mpMap")
		indices <- matrix(nrow=maxOffset *nrow(rf) - maxOffset * (maxOffset + 1) / 2, ncol = 2)
		counter <- 1
		for(offset in 1:maxOffset)
		{
			for(i in 1:(nrow(rf)-offset))
			{
				indices[counter,] <- c(i+offset, i)
				counter <- counter + 1
			}
		}
		#B vector for nnls
		b <- rf[indices]
		b[b == 0.5] <- 0.49
		if (!requireNamespace("nnls", quietly=TRUE))
			stop("nnls needed for computemap to work. Please install it.\n", call.=FALSE) 
		result <- nnls::nnls(d, mf(b)) 
		object$map[[chr]] <- c(0, cumsum(result$x[which(indices[,1] == indices[,2]+1)]))
		names(object$map[[chr]]) <- colnames(object$finals)[m]
		#m <- match(names(object$map[[chr]]), colnames(object$finals))
		#rf <- fill(fill(object$rf$theta[m,m], missfx), 1)
		#rf[rf==.5] <- .49
		#object$map[[chr]] <- c(0, cumsum(mf(rf[row(rf)==(col(rf)+1)])))
		#names(object$map[[chr]]) <- colnames(object$finals)[m]
	}
	return(object)
}
designMat <- function(n, maxOffset)
{
	resultMat <- matrix(0, nrow=n*maxOffset - maxOffset*(maxOffset - 1)/2, ncol=n)
	#column of matrix
	for(i in 1:n)
	{
		offset <- 1
		#j is the section going by rows
		#for(j in 1:n)
		for(j in 1:maxOffset)
		{
			resultMat[offset + max(0, i-j):min(n-j, i-1) ,i] <- 1
			offset <- offset + (n-j+1)
		}
	}
	return(resultMat)
}
