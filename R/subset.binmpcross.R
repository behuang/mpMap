#' Subset binned mpcross object
#'
#' Reduces a binned mpcross object down to a specified set of chromosomes, markers and/or lines
#' @export subset.binmpcross 
#' @method subset binmpcross
#' @param x Object of class \code{binmpcross}
#' @param chr Selected chromosomes TO KEEP. Default is all
#' @param groups selected groups to KEEP. Default is all
#' @param markers Selected markers TO KEEP. Default is all
#' @param lines Selected lines TO KEEP. Default is all
#' @param ... Additional arguments
#' @note This is essentially the same function as subset.mpcross with the addition of dealing with the bins component. 
#' @return The original object with chromosomes/lines/markers removed which are not listed in the arguments.
#' @seealso \code{\link[mpMap]{mpcollapse}}, \code{\link[mpMap]{mpexpand}}

subset.binmpcross <-
function(x, groups=NULL, chr=NULL, markers=NULL, lines=NULL, ...)
{
  if (all(sapply(c(chr, markers, lines, groups), length)==0)) return(x)
  if (sum(sapply(c(chr, markers, groups), length) > 1))
  {
	stop("Only one of chr, markers and groups can be input")
  }
 
  output <- NextMethod(x)
  output$bins <- x$bins

  ## Now do something special to handle the bins - needs to be a matrix of
  ## MarkerName bin group
  binmarkers <- output$bins$binMarkerName
  keepbm <- which(binmarkers %in% colnames(output$finals))
  output$bins <- output$bins[keepbm,]
  binmarkers <- binmarkers[keepbm]

  ## Now just need to deal with any reordering
  output$bins <- output$bins[order(match(binmarkers, colnames(output$finals))),]
  
  ## need to recode the groups to match what exists in the object
  output$bins$group <- match(output$bins$group, unique(output$bins$group))
  output
}

