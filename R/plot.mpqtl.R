#' Plot output from interval mapping with detected QTL
#'
#' Plot -log10(p-value) or test statistic against cM position for (composite) interval mapping in multi-parent crosses. QTL support intervals are indicated with rectangles surrounding peaks. 
#' @importFrom graphics plot
#' @importFrom stats na.exclude
#' @importFrom graphics rect
#' @export 
#' @method plot mpqtl
#' @param x Object of class \code{mpqtl}
#' @param wald Flag for whether to plot the Wald statistic or -log10(p)
#' @param chr Set of chromosomes to plot
#' @param lodsupport x-LOD support interval plotted in green
#' @param ... Additional arguments to plotting function
#' @return Plots the -log10(p) or Wald statistic for all chromosomes against the total genome in cM. QTL support intervals are indicated with shaded rectangles surrounding peaks 
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- qtl::sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, 
#'		qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), 
#' 		nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' plot(mpq.dat)

plot.mpqtl <- function(x, wald=FALSE, chr, lodsupport=1, ...)
{
  output <- list()
  map <- attr(x$prob, "map")
  if (missing(chr)) chr <- 1:length(map)
  if (is.factor(chr)) chr <- as.character(chr)
  if (is.numeric(chr)) chr <- names(map)[chr]

  ## set up as a scanone output
  psc <- data.frame(chr=rep(names(x$QTLresults$wald), unlist(lapply(x$QTLresults$wald, length))), pos=unlist(map), lod=-log10(unlist(x$QTLresults$pval)))
  waldsc <- data.frame(chr=rep(names(x$QTLresults$wald), unlist(lapply(x$QTLresults$wald, length))), pos=unlist(map), lod=unlist(x$QTLresults$wald))

  waldsc[,1] <- factor(waldsc[,1], levels=unique(waldsc[,1]))
  psc[,1] <- factor(psc[,1], levels=unique(psc[,1]))
  class(psc) <- class(waldsc) <- c("scanone", "data.frame")

  if (any(psc[,3]==Inf)) { 
	wald <- TRUE
	cat("Some p-values=0; plotting Wald score\n")
    } 

  if (wald==TRUE) plot(waldsc, chr=chr, ylab="Wald", ...) else plot(psc, chr=chr, ylab="-log10(p)", ...)
  output$waldsc <- waldsc
  output$psc <- psc
  
  map2 <- map[chr]
  qtlpos <- vector()
  chrpos <- c(0,cumsum(unlist(lapply(map2, function(x) diff(range(x))))+25))
  for (i in na.exclude(match(names(x$QTLresults$qtl), names(map2))))
	qtlpos <- c(qtlpos, x$QTLresults$qtl[[names(map2)[i]]][,1]+chrpos[i])

#  si <- supportinterval(x, lodsupport=lodsupport)$support[, which(names(x$QTLresults$qtl) %in% names(map2)),drop=FALSE]
   si <- supportinterval(x, lodsupport=lodsupport)
   si <- si$support[, which(names(si$qtlpos)%in% names(map2)), drop=FALSE]

  
  xqtl <- x$QTLresults$qtl
  if (length(intersect(chr, names(xqtl)))>0) {
  qtlmat <- do.call("rbind", xqtl[na.exclude(match(chr, names(xqtl)))])
  rownames(qtlmat) <- rep(names(xqtl)[na.exclude(match(chr, names(xqtl)))], unlist(lapply(xqtl[na.exclude(match(chr, names(xqtl)))], nrow)))
  for (i in 1:length(qtlpos)) {
  	chrnam <- vector()
	 m <- match(rownames(qtlmat)[i], names(map2))
	 chrstart <- nrow(waldsc) - max((nrow(waldsc)-(1:nrow(waldsc)))*(waldsc[,1]==rownames(qtlmat)[i]))
	 chrend <- max((1:nrow(waldsc))*(waldsc[,1]==rownames(qtlmat)[i]))

	rect(xleft=chrpos[m]+si[1,i], xright=si[2,i]+chrpos[m], ybottom=-1000, ytop=1000, col="#00990022")
  }

  output$support <- si
  }

  invisible(output)
}
