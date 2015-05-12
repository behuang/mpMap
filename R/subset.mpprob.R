#' Subset mpprob object
#'
#' Reduces an mpprob object down to a specified set of chromosomes, markers and/or lines
#' @export subset mpprob
#' @method subset mpprob
#' @param x Object of class \code{mpprob}
#' @param chr Selected chromosomes TO KEEP. Default is all
#' @param markers Selected markers TO KEEP. Default is all
#' @param lines Selected lines TO KEEP. Default is all
#' @param ... Additional arguments
#' @note Chromosomes can be input either as the character names of chromosomes or the index of the chromosomes in the map. Markers can be input as either character values matching the colnames of x$finals, or indices of columns in that matrix. Note that if markers are removed, the founder probabilities will be recomputed for the new map with previous settings for mpprob. Lines can be input as either character values (matching the rownames of x$finals) or indices of rows in that matrix. 
#' @return The original object with chromosomes/lines removed which are not listed in the arguments.
#' @seealso \code{\link[mpMap]{mpprob}}

subset.mpprob <-
function(x, chr=NULL, markers=NULL, lines=NULL, ...)
{
   n.founders <- nrow(x$founders)

   if (all(sapply(c(chr, markers, lines), length)==0)) return(x)
   if(!is.null(markers) && !is.null(chr)) stop("Inputs chr and markers cannot be used simultaneously")

   output <- subset.mpcross(x, chr=chr, markers=markers, lines=lines)

	if (!is.null(chr)) 
	{
		output$estfnd <- output$estfnd[chr]
		#subsetting drops attributes. So use subassignment instead
		#output$prob <- output$prob[chr]
		if(is.character(chr)) chr <- match(chr, names(output$prob))
		if(any(is.na(chr))) stop("Internal error")
		output$prob[setdiff(1:length(chr), chr)] <- NULL
		attr(output$prob, "map") <- attr(output$prob, "map")[chr]
		
		#Copy all other attributes
		attr(output$prob, "step") <- attr(x$prob, "step")
		attr(output$prob, "program") <- attr(x$prob, "program")
		attr(output$estfnd, "threshold") <- attr(x$estfnd, "threshold")
		attr(output$prob, "mapfx") <- attr(x$prob, "mapfx")
	}
 
	if (!is.null(lines)) 
	{
		if (is.character(lines)) linnum <- match(lines, rownames(x$finals)) else linnum <- lines
		output$estfnd <- lapply(output$estfnd, function(x) return(x[linnum,]))
		output$prob <- lapply(output$prob, function(x) return(x[linnum,]))
		
		#Copy all attributes
		attr(output$prob, "step") <- attr(x$prob, "step")
		attr(output$prob, "program") <- attr(x$prob, "program")
		attr(output$estfnd, "threshold") <- attr(x$estfnd, "threshold")
		attr(output$prob, "mapfx") <- attr(x$prob, "mapfx")
		attr(output$prob, "map") <- attr(x$prob, "map")
	}

	if (!is.null(markers))
	{
		if(attr(x$prob, "step") == 0)
		{
      for (i in names(x$map)) {
        removedMarkers <- setdiff(names(x$map[[i]]), names(output$map[[i]]))
        removedMarkers <- match(removedMarkers, names(x$map[[i]]))
        if (length(removedMarkers)>0) {
    			output$estfnd[[i]] <- output$estfnd[[i]][, -removedMarkers]
          output$prob[[i]] <- output$prob[[i]][, -(rep((removedMarkers-1)*n.founders, each=n.founders)+rep(1:4, length(removedMarkers)))]
        }
        if (length(removedMarkers)==length(x$map[[i]])) {
          output$estfnd[[i]] <- NULL
          output$prob[[i]] <- NULL
        }
      }
			#Copy attributes. The map was subsetted by the subset.mpcross function
			attr(output$prob, "map") <- output$map
			
			attr(output$prob, "step") <- attr(x$prob, "step")
			attr(output$prob, "program") <- attr(x$prob, "program")
			attr(output$estfnd, "threshold") <- attr(x$estfnd, "threshold")
			attr(output$prob, "mapfx") <- attr(x$prob, "mapfx")
		}
		else
		{	
			output <- mpprob(output, step=attr(x$prob, "step"), program=attr(x$prob, "program"), threshold=attr(x$estfnd, "threshold"), mapfx=attr(x$prob, "mapfx"))
		}
	}
	return(output)  
}

