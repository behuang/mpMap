#' Subset mpcross object
#'
#' Reduces an mpcross object down to a specified set of chromosomes, markers and/or lines
#' @export 
#' @method subset mpcross
#' @param x Object of class \code{mpcross}
#' @param chr Selected chromosomes TO KEEP. Default is all
#' @param groups selected groups to KEEP. Default is all
#' @param markers Selected markers TO KEEP. Default is all
#' @param lines Selected lines TO KEEP. Default is all
#' @param ... Additional arguments
#' @note Chromosomes can be input either as the character names of chromosomes or the index of the chromosomes in the map. Markers can be input as character names or the index in the matrix x$finals. Lines can be input as either character values (matching the rownames of x$finals) or indices of rows in that matrix. Groups must be input as numbers corresponding to the specific groups to keep (the groups present in the object may not form a consecutive set of numbers). Note that only one of chr, groups and markers can be selected. 
#' @return The original object with chromosomes/lines/markers removed which are not listed in the arguments.
#' @seealso \code{\link[mpMap]{mpcross.object}}
#' @examples
#' map <- qtl::sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=ped, 
#'		qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), 
#'		nrow=2, ncol=6, byrow=TRUE), seed=1)
#' sim.dat
#' red.dat <- subset(sim.dat, chr=1, lines=1:50)
#' red.dat

subset.mpcross <-
function(x, groups=NULL, chr=NULL, markers=NULL, lines=NULL, ...)
{
  if (all(sapply(c(chr, markers, lines, groups), length)==0)) return(x)
  if (sum(sapply(c(chr, markers, groups), length) > 1))
  {
	stop("Only one of chr, markers and groups can be input")
  }
 
  output <- x

  if (!is.null(chr)) {
    if (is.numeric(chr)) chr <- names(x$map)[chr]
    output$map <- as.list(output$map[chr])
	class(output$map) <- "map"

    if (is.null(x$lg)) x <- mpgroup(x)
    groups <- match(chr, names(x$map))

    for (ii in chr)
      markers <- unique(c(markers, names(x$map[[ii]])))
  }

  if(!is.null(groups))
  {
	if(!("lg" %in% names(x)))
		stop("If groups is specified, input mpcross object must contain an entry named lg created by mpgroup")
	
	bool <- x$lg$groups %in% groups
	names(bool) <- names(x$lg$groups)
	markers <- which(bool)
  }  

  if (!is.null(markers)) {
    	if (is.character(markers)) 
      	  mrknum <- match(markers, colnames(x$finals)) else mrknum <- markers

	if (is.numeric(markers)) {
	    mrknum <- markers
	    markers <- colnames(x$finals)[mrknum]
    	}

    	output$founders <- as.matrix(output$founders[,mrknum])
    	output$finals <- as.matrix(output$finals[,mrknum])
    	colnames(output$founders) <- colnames(output$finals) <- markers

    	if (!is.null(x$rf)) {
	  m2 <- match(markers, colnames(output$rf$theta))
 	  output$rf[1:3] <- lapply(output$rf[1:3], function(x) x[m2,m2])
    	}
    	if (!is.null(x$ld)) {
	  m2 <- match(markers, colnames(output$ld$W))
	  output$ld <- lapply(output$ld, function(x) x[m2,m2])
    	}

  #If there's a map, we MIGHT keep it, unless things get reordered. In which case we drop the map and put in an lg structure
 	if (!is.null(output$map))
	{   
		#If there's a map and we're just subsetting stuff, then keep it. If we're reordering markers though, strip it out. 
		isReordering <- FALSE
		for(chr in names(output$map))
		{
			currentChrMarkers <- names(output$map[[chr]])
			
			#Get out the markers on this chromosome, in the final order
			currentKeptMarkers <- markers[markers %in% currentChrMarkers]
			#If when reordered there are any decreases in position, we must have reordered stuff. In that case remove make
			if(any(diff(output$map[[chr]][currentKeptMarkers]) < 0))
			{
				isReordering <- TRUE
				break
			}
		}
		if(!isReordering)
		{
		 	output$map <- lapply(x$map, function(y) return(y[which(names(y) %in% markers)]))
			class(output$map) <- "map"
			output$map[which(lapply(output$map, length)==0)] <- NULL
		}
		else
		{
			#If we can't keep the map, put in an LG structure
			output$lg <- list(groups = rep(1:length(output$map), unlist(lapply(output$map, length))), all.groups = 1:length(output$map))
			names(output$lg$groups) <- markers
			output$map <- NULL
		}
	}
	#If we kept the map, strip out any lg structure
#	if(!is.null(output$map)) output$lg <- NULL
	
	#If there's still an LG structure, reorder it
        if (!is.null(output$lg)) 
	{
		m2 <- match(markers, names(output$lg$groups))
		output$lg$groups <- x$lg$groups[m2]
		output$lg$all.groups <- unique(output$lg$groups)
        } 
  }

  if (!is.null(lines)) {
    if (is.character(lines)) linnum <- match(lines, rownames(x$finals)) else linnum <- lines
	  output$id <- output$id[linnum]
	  output$pedigree[,4] <- rep(0, nrow(output$pedigree))
	  output$pedigree[output$id,4] <- 1
      output$finals <- output$finals[linnum,, drop=FALSE]
    if (!is.null(output$pheno)){
      	output$pheno <- output$pheno[linnum,, drop=FALSE]
 	rownames(output$pheno) <- rownames(output$finals)
    }
  }
  output
}

