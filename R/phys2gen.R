#' Estimate a genetic map based on physical map ordering in multi-parent cross
#'
#' Use an input physical map to group and order markers; estimate positions 
#' based on recombination fractions between adjacent markers
#'
#' @export
#' @importFrom stats approx
#' @param physmap Physical map to be used for groups and orders of markers
#' @param mpcross Mpcross objects from which to estimate recombination fractions
#' @param mapfx Map function to convert recombination fractions to cM 
#' @param maxrf Maximum rf to allow between markers before removing from map
#' 
#' @details
#' For each pair of adjacent markers recombination fractions are estimated, 
#' and the map positions computed as the sum of adjacent cM values. The
#' returned object is of class "pgmap", which contains both physical and 
#' genetic map information for ease of later comparison. 
#' 
#' @return A pgmap object is returned containing both genetic and physical map positions. 

phys2gen <- function(physmap, mpcross, mapfx, maxrf=.25)
{
  output <- list()

  if (missing(mapfx)) mapfx <- "haldane"
  if (mapfx=="haldane") mf <- haldaneR2X else mf <- kosambiR2X

  ## Check that there are the same markers in both cases
  gmrk <- colnames(mpcross$finals)
  pmrk <- unlist(lapply(physmap, names))
  if (length(setdiff(pmrk, gmrk))>0) 
	cat("Markers in physical map with no genetic data, will interpolate genetic map positions for those without data \n")

  if (length(setdiff(gmrk, pmrk))>0) {
	cat("Markers in genetic data but not physical map, will not be included in genetic map \n")
  }
  mpcross <- subset(mpcross, markers=intersect(gmrk, pmrk))
  gmrk <- colnames(mpcross$finals)
 
  ### Estimate recombination fractions between all pairs of adjacent markers
  pmap2 <- physmap
  gmap <- physmap
  badmrk <- list()
  ## remove markers which are not in the genetic data
  for (i in 1:length(pmap2)) {
    gmap[[i]][1:length(gmap[[i]])] <- NA
    mrk <- intersect(names(pmap2[[i]]), gmrk)
    pmap2[[i]] <- pmap2[[i]][mrk]
    rf <- vector(length=length(pmap2[[i]])-1)
    for (j in 1:length(rf)) {
	sub <- subset(mpcross, markers=names(pmap2[[i]])[j+0:1])
	rf[j] <- mpestrf(sub)$rf$theta[1,2]
     }
    ## need to deal with some values being unidentifiable
    ## remove values of 0.5 and missing - flag values of 0.5 as potentially bad
    index <- names(pmap2[[i]])
    ### once we drop some markers out need to estimate gaps
    while(any(is.na(rf)) | any(rf>maxrf)) {
      tokeep <- which(!(is.na(rf)) & rf<maxrf)
      rf <- rf[tokeep]
      tokeep <- c(tokeep, length(index))
      for (j in which(diff(tokeep)>1)) {
	sub <- subset(mpcross, markers=index[tokeep[j+0:1]])
        rf[j] <- mpestrf(sub)$rf$theta[1,2]
      }
      index <- index[tokeep]
    }
    cm <- mf(rf)
    gmap[[i]][match(index, names(gmap[[i]]))] <- c(0, cumsum(cm))
    ## now interpolate for rest of the positions
    ap <- approx(x=physmap[[i]], y=gmap[[i]], xout=physmap[[i]])
    gmap[[i]] <- ap$y
    names(gmap[[i]]) <- names(physmap[[i]])
  }
  output$Gmap <- gmap
  output$Pmap <- physmap
  class(output) <- c("pgmap", "map")
  output
}


