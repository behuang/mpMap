#' Impute unknown recombination fraction values
#' 
#' Construct a new mpcross object, with missing recombination fraction values imputed.  
#' @export
#' @param mpcross Object of class \code{mpcross}
#' @details 
#' Imputation is performed by replacing missing values with a value from the `nearest' column. Here the distance between two columns of the recombination fraction matrix is determined by treating both columns as points in euclidean space and taking the euclidean distance between these points. In theory using rows instead of columns may result in a different imputed value, but in practice the difference between the two approaches is probably negligible. 
#' 
#' @return The original object with any NA values of rf$theta replaced by imputed values

mpimputerf <- function(mpcross)
{
	if(missing(mpcross))
	{
		stop("Input mpcross cannot be missing")
	}
	if(!("rf" %in% names(mpcross)))
	{
		stop("Input mpcross must contain an entry named 'rf'")
	}
	if(!all(c("theta", "lod", "lkhd") %in% names(mpcross$rf)))
	{
		stop("Input mpcross must contain entries rf$theta, rf$lod and rf$lkhd")
	}
	if(all(c("lg", "map") %in% names(mpcross)))
	{
		warning("Both lg and map entries present in input mpcross object. Using map entry, discarding lg entry")
	}
	#If we're using map entry, construct lg entry and delete it later. This only issues a warinng if there was an existing lg entry which is getting wipped.
	discard.lg <- FALSE
	if("map" %in% names(mpcross))
	{
		mpcross$lg <- list()
		mpcross$lg$groups <- vector(mode="numeric", length=sum(unlist(lapply(mpcross$map, length))))
		mpcross$lg$all.groups <- 1:length(mpcross$map)
		names(mpcross$lg$groups) <- unlist(lapply(mpcross$map, names))
		for(i in 1:length(mpcross$map))
		{
			mpcross$lg$groups[names(mpcross$map[[i]])] <- i
		}
		discard.lg <- TRUE
	}
	if(!("lg" %in% names(mpcross)))
	{
		stop("Input mpcross must contain an entry named 'lg'")
	}
	if(!all(c("groups", "all.groups") %in% names(mpcross$lg)))
	{
		stop("Input mpcross must contain entries lg$groups and lg$all.groups")
	}
	
  if (sum(is.na(mpcross$rf$theta))>0)
	  ret <- .Call("impute", mpcross, PACKAGE="mpMap")
  else ret <- mpcross
	#Imputation .Call function only works on columns and won't give a symmetric result. So symmetrise..
	ret$rf$theta <- (ret$rf$theta + t(ret$rf$theta))/2
	ret$rf$lod <- (ret$rf$lod + t(ret$rf$lod))/2
	ret$rf$lkhd <- (ret$rf$lkhd + t(ret$rf$lkhd))/2
	for(group in mpcross$lg$all.groups)
	{
		indices <- which(mpcross$lg$groups == group)
		if(any(is.na(ret$rf$theta[indices, indices]))) stop("Internal error")
	}
	if(discard.lg)
	{
		ret$lg <- NULL
	}
	return(ret)
}
