#' Function to combine two data sets which have the same pedigree and genetic lines, but different (potentially overlapping) marker data.
#' 
#' @export mpadd
#' @param mpcross1 First of the two objects to be combined
#' @param mpcross2 Second of the two objects to be combined
#' @param gpu Should a GPU be used to do any required recombination fraction calculations, if available. Defaults to \code{FALSE}.
#' @return An object of class \code{mpcross}
#' 
#' If recombination fractions have been calculated for either of the input objects, the recombination fractions will also be calculated, for all the markers in the final object. Similarly, if LD has been calculated for either of the input objects, LD will be calculated for all markers in the final object.
#' Previous calculations are re-used in the process, so this is much more efficient than discarding the recombination fraction / LD data from both objects and recomputing again. 
#'
mpadd <- function(mpcross1, mpcross2, gpu)
{
	if(missing(mpcross1) || missing(mpcross2))
	{
		stop("Inputs mpcross1 and mpcross2 are required")
	}
	if(missing(gpu)) gpu <- FALSE

	nmrk1 <- ncol(mpcross1$founders)
	nmrk2 <- ncol(mpcross2$founders)
	#Drop common markers from object 2
	markers.intersect <- intersect(colnames(mpcross1$founders), colnames(mpcross2$founders))
	remainingMarkers <- setdiff(colnames(mpcross2$finals), markers.intersect)
	if(length(remainingMarkers) == 0) 
	{
		warning("Input object 2 did not contain any new markers. Returning object 1 unchanged")
		return(mpcross1)
	}
	mpcross2 <- subset(mpcross2, markers=remainingMarkers)

	#Check for inconsistent ID names
	if(length(mpcross1$id) != length(mpcross2$id) || any(sort(mpcross1$id) != sort(mpcross2$id))) stop("Different lines in mpcross1 and mpcross2")

	#Do we have to re-order the lines though?
	if(any(mpcross1$id != mpcross2$id))
	{
		mpcross2$finals <- mpcross2$finals[match(mpcross2$id, mpcross1$id),]
		mpcross2$id <- mpcross1$id
	}

	#Check for inconsistent founder names
	if(length(mpcross1$fid) != length(mpcross2$fid) || any(sort(mpcross1$fid) != sort(mpcross2$fid))) stop("Different founder lines in mpcross1 and mpcross2")
	#Reorder founders if necessary
	if(any(mpcross1$fid != mpcross2$fid))
	{
		mpcross2$founders <- mpcross2$founders[match(mpcross2$fid, mpcross1$fid),]
		mpcross2$fid <- mpcross1$fid
	}

	#Check for same pedigree
	if(any(dim(mpcross1$pedigree) != dim(mpcross2$pedigree)) || any(mpcross1$pedigree != mpcross2$pedigree))
	{
		stop("Input objects have differing pedigrees")
	}
	output <- mpcross1
	output$finals <- cbind(mpcross1$finals, mpcross2$finals)
	output$founders <- cbind(mpcross1$founders, mpcross2$founders)
	firstHasRF <- "rf" %in% names(mpcross1)
	secondHasRF <- "rf" %in% names(mpcross2)

	firstHasLD <- "ld" %in% names(mpcross1)
	secondHasLD <- "ld" %in% names(mpcross2)
	if(firstHasRF || secondHasRF)
	{
		if(firstHasRF && secondHasRF)
		{
			if(!all(mpcross1$rf$r == mpcross2$rf$r)) stop("Set of tested recombination fraction values was inconsistent across input objects. These objects cannot be combined. ")
			r <- mpcross1$rf$r
		}
		if(!firstHasRF)
		{
			mpcross1 <- mpestrf(mpcross1, mpcross2$rf$r, gpu=gpu)
			r <- mpcross2$rf$r
		}
		if(!secondHasRF)
		{
			mpcross2 <- mpestrf(mpcross2, mpcross1$rf$r, gpu=gpu)
			r <- mpcross1$rf$r
		}
		#Convert inputs to .Call to standardised forms / types
		output$pedigree <- as.matrix(output$pedigree)
		if(class(output$founders) != "matrix") output$founders <- as.matrix(output$founders)
		if(class(output$finals) != "matrix") output$finals <- as.matrix(output$finals)
		if(class(output$pedigree) != "matrix") output$pedigree <- as.matrix(output$pedigree)
		
		if(mode(output$founders) != "integer") mode(output$founders) <- "integer"
		if(mode(output$finals) != "integer") mode(output$finals) <- "integer"
		if(mode(output$fid) != "integer") mode(output$fid) <- "integer"
		if(mode(output$pedigree) != "integer") mode(output$pedigree) <- "integer"
		
		if(class(output$id) != "integer") output$id <- as.integer(output$id)

		offDiagonal <- .Call("rfhaps", list(output), r, c(1,ncol(mpcross1$finals)+1), c(ncol(mpcross1$finals)+1, ncol(output$finals)+1), list(rep(1, nrow(output$finals))), gpu, -2, PACKAGE="mpMap")
		output$rf$lkhd <- rbind(cbind(mpcross1$rf$lkhd, offDiagonal$lkhd), cbind(t(offDiagonal$lkhd), mpcross2$rf$lkhd))
		output$rf$lod <- rbind(cbind(mpcross1$rf$lod, offDiagonal$lod), cbind(t(offDiagonal$lod), mpcross2$rf$lod))
		output$rf$theta <- rbind(cbind(mpcross1$rf$theta, offDiagonal$theta), cbind(t(offDiagonal$theta), mpcross2$rf$theta))
		
		if (firstHasLD || secondHasLD)
		{
			if(!secondHasLD) mpcross2 <- mpcalcld(mpcross2)
			if(!firstHasLD) mpcross1 <- mpcalcld(mpcross1)
			output$ld <- combine_ld(mpcross1, mpcross2, output$rf$theta[1:nmrk1, nmrk1+1:nmrk2])
		}
	}
	
	return(output)
}

