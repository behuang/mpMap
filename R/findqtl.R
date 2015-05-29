#' Detect a QTL peaks in a QTL profile from (composite) interval mapping
#'
#' Given the output from a scan of a chromosome, locates local maxima exceeding a significance threshold in a QTL profile. 
#' @export
#' @param mpqtl Object of class \code{mpqtl}
#' @param dwindow Window over which to smooth p-values - default is five markers
#' @param threshold Threshold peaks must exceed to be detected (-log10(p))
#' @return The original input object with additional entries for newly detected QTL. 
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{plot.mpqtl}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' mpq2 <- findqtl(mpq.dat, dwindow=5, threshold=3)
#' plot(mpq2)
#' summary(mpq2)

findqtl <- function(mpqtl, dwindow=5, threshold)
{
  output <- mpqtl 
  output$QTLresults$qtl <- NULL
  map <- attr(mpqtl$prob, "map")

  ## Need to form scanone object from mpqtl. 
  sc <- data.frame(chr=rep(names(mpqtl$QTLresults$wald), unlist(lapply(mpqtl$QTLresults$wald, length))), pos=unlist(map), lod=-log10(unlist(mpqtl$QTLresults$pval)))
  sc[,1] <- factor(sc[,1], levels=unique(sc[,1]))
  class(sc) <- c("scanone", "data.frame")

  require(VPdtw)
  index <- list()
  for (j in unique(sc$chr)) {
    startind <- min(which(sc$chr==j))
    y <- sc$lod[sc$chr==j]
    x <- sc$pos[sc$chr==j]
    ap <- approx(x, y, xout=unique(sort(c(x, seq(min(x), max(x), .1)))))
    ts <- ap$y
    tmp1 <- dilation(ts, dwindow)
    tmp2 <- which(ts==tmp1)
    tmp3 <- tmp2[which(ts[tmp2]>threshold)]
  
    if (length(tmp3)>0) {
    locmax <- vector(length=length(tmp3))
    for (i in 1:length(tmp3)) {
      lower <- max(tmp3[i]-dwindow, 1)
      upper <- min(tmp3[i]+dwindow, length(tmp1))
      locmax[i] <- (tmp1[lower] <= tmp1[tmp3[i]] & tmp1[upper] <= tmp1[tmp3[i]])
    }
  
    dlm <- c(diff(tmp3), dwindow)
    locmax[dlm<dwindow] <- FALSE
    ## need to rematch the tmp3 indices up to the original ones. 
    oldindex <- vector(length=length(tmp3[locmax]))
    for (k in 1:length(oldindex)) oldindex[k] <- which.min(abs(ap$x[tmp3[locmax][k]]-x))

    index[[j]] <- oldindex+startind-1
    nqtl <- length(index[[j]])
    output$QTLresults$qtl[[j]] <- cbind(sc[index[[j]], 2], matrix(mpqtl$QTLresults$fndrfx[[j]][, oldindex], nrow=nqtl, byrow=T), matrix(mpqtl$QTLresults$se[[j]][, oldindex], nrow=nqtl, byrow=T))
    } else index[[j]] <- NULL
  }
  attr(output$QTLresults$qtl, "index") <- index
  attr(output$QTLresults$qtl, "nqtl") <- length(unlist(index))
  

  output
}
