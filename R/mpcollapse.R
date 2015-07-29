#' Collapse an mpcross object into a binned object
#' 
#' Collapses an mpcross object into bins. Within bins, markers are grouped by founder distribution pattern and missing genotypes imputed. Haplotypes are formed for markers with unique FDPs and determine the bin alleles. Information on which markers belong to each bin is retained within the output object. 
#' @export
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom stats cutree
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @param method Choice of whether to bin based on recombination fraction or correlation 
#' @param adj Flag for whether require bins to be based on adjacent markers only 
#' @param cutoff Max RF allowed before starting a new bin
#' @param consensusProbability Proportion of genotypes within FDP groups which must match for imputation to occur. 
#' @param missingCutoff Proportion of missing data allowed within markers before not using to create haplotypes
#' @return A binned object of class \code{binmpcross}
#' @seealso \code{link{mpexpand}}
#' @note Markers are assigned to a bin if there is zero recombination between any of the markers within the bin (or correlation of 1). Next, a three-step process is performed to create the new set of bin markers. First, markers within a bin are grouped based on founder distribution pattern (FDP). Second, for each FDP, a single representative marker is chosen, with missing values imputed from all markers sharing that FDP. Third, haplotypes formed from these representative markers are used to replace genotypes for the bin marker. 

mpcollapse<-function(object, method=c("rf", "cor"), adj=FALSE, cutoff=0, consensusProbability=.75, missingCutoff=.25) {

  enlist <- function(X) lapply(1:ncol(X), function(j) X[,j])
  
  Matches <- function(x, X) {
    stopifnot(ncol(X) == length(x))
    
    mv <- which(is.na(X), arr.ind = TRUE)
    X[mv] <- x[mv[,2]]
    match(do.call(paste0, enlist(X)), paste(x, collapse=""), 0L)
  }
    
  ## Collapse a matrix of similar values down to one column with missing
  ## values imputed from all the rest
  formConsensus <- function(X, consensusProbability) {
    vec <- apply(X, 1, function(x) if (max(table(x))/sum(table(x))>consensusProbability) names(table(x))[which.max(table(x))] else NA)  
    vec
  }

  if (missing(method)) method <- "rf"

  if (is.null(object$rf$theta) & method=="rf") 
	stop("Need to estimate recombination fractions first\n")

  if (is.null(object$lg) & is.null(object$map)) {
	cat("Assuming all markers are on the same linkage group\n") 
	object <- mpgroup(object, groups=1) 
  }

  nFounders <- nrow(object$founders)
  nFinals <- nrow(object$finals)
  
  if (is.null(object$lg)) object <- mpgroup(object, groups=length(object$map))

  ## Note that correlation is not ideal as it does not take into
  ## account founder genotypes
  if (method=="cor")
	  dist <- (cor(object$finals, use="pairwise.complete.obs"))^2
  else
	  dist <- mpimputerf(object)$rf$theta
  
  n.chr <- length(unique(object$lg$all.groups))
  nmrkperchr <- table(object$lg$groups)

## Note that we remove recombination fractions - these will need to be 
## re-estimated based on the binned markers. Linkage groups we will 
## just collapse down since we're retaining that information with the bins
## Maps will be removed as the distances may change based on the binned markers. 
  output <- object

  if (!is.null(object$rf)) cat("RF need to be re-estimated based on binned markers\n")
  if (!is.null(object$map)) cat("Map need to be re-estimated based on binned markers\n") 

  binFounders <- binFinals <- list()
  binInfo <- list()
  
  for (i in 1:n.chr) {
    indexGroup <- which(object$lg$groups==i)
    
    ## create bins
    mat <- dist[indexGroup, indexGroup]
    
    if (!adj) {
      cl.obj <- hclust(as.dist(mat), method="complete")
      cl.obj$height <- round(cl.obj$height, 5)
      binInfo[[i]] <- cutree(cl.obj, h=cutoff)
      nbins <- max(binInfo[[i]])
    } else {
      ## how are we going to construct bins in this case?
      binInfo[[i]] <- 1:ncol(mat)
      for (k in 2:ncol(mat)) 
        if (mat[k-1, k]==0) binInfo[[i]][k] <- binInfo[[i]][k-1]
      binInfo[[i]] <- match(binInfo[[i]], names(table(binInfo[[i]])))
      nbins <- length(table(binInfo[[i]]))
    }
    ## Note that after this there will be a substantial amount of recoding. Do not return original object
    ## within each bin recode genotypes to haplotypes
    binFounders[[i]] <- matrix(nrow=nFounders, ncol=nbins)
    binFinals[[i]] <- matrix(nrow=nFinals, ncol=nbins)
    colnames(binFounders[[i]]) <- colnames(binFinals[[i]]) <- paste("C", i, "B", 1:nbins, sep="")
    rownames(binFounders[[i]]) <- rownames(object$founders)
    rownames(binFinals[[i]]) <- rownames(object$finals)
    
    for (j in 1:nbins) {
      indexBin <- indexGroup[which(binInfo[[i]]==j)]    
      if(length(indexBin)>1) {
##### Add in extra step here to pull out unique FDPs and impute missingness
	finstart <- object$finals[, indexBin]
	foustart <- object$founders[, indexBin]
	fdp <- unique(foustart, MARGIN=2)
	fin <- matrix(nrow=nFinals, ncol=ncol(fdp))
	if(ncol(fdp)<ncol(foustart))
	  for (k in 1:ncol(fdp)) {
	     cols <- which(apply(foustart, 2, function(x) all(x==fdp[,k]))==T)
	     if (length(cols)>1)
	      	vec <- formConsensus(finstart[, cols], consensusProbability) else vec <- finstart[,cols]
	     fin[,k] <- as.numeric(vec)
	  } else fin <- finstart

	## If any of the markers have too much missing data, don't use
	## them to form haplotypes (even if you end up with less informative 
	## markers) - better off having more values
	nm <- apply(fin, 2, function(x) sum(is.na(x))/length(x))
	if (all(nm>missingCutoff)) touse <- 1:ncol(fin) else touse <- which(nm<missingCutoff) 
###### next go back to previous routine, taking out missingness imputation
        uniqueHaps <- unique(fdp[, touse, drop=F], MARGIN=1)
        missFinals <- which(apply(fin[, touse, drop=F], 1, function(x) any(is.na(x))))
    #    uniqueMissFinals <- unique(object$finals[missFinals, indexBin], MARGIN=1)
        matches <- matrix(nrow=length(missFinals), ncol=nrow(uniqueHaps))
        for (k in 1:nrow(uniqueHaps)) 
            matches[,k] <- Matches(uniqueHaps[k,], fin[missFinals, touse,drop=F])
        for (k in 1:nrow(uniqueHaps)) {
	    toreplace <- which(matches[,k]==1 & rowSums(matches)==1) 
           fin[missFinals[toreplace], touse] <- matrix(rep(uniqueHaps[k,], length(toreplace)), ncol=ncol(uniqueHaps), byrow=TRUE)
        }
#        matchHap <- vector(length=nrow(uniqueMissFinals))
#        for (k in 1:nrow(uniqueMissFinals)) {
#            tmp <- apply(uniqueHaps, 1, function(x) all(x[which(!is.na(uniqueMissFinals[k,]))]==uniqueMissFinals[k,which(!is.na(uniqueMissFinals[k,]))]))
#            if (sum(tmp)==1) matchHap[k] <- which(tmp) else matchHap[k] <- NA
#            if (!is.na(matchHap[k])) {
#        matchUnique <- apply(object$finals[missFinals, indexBin], 1, 
#                  function(x) all(x[which(!is.na(uniqueMissFinals[k,]))]==uniqueMissFinals[k,which(!is.na(uniqueMissFinals[k,]))]))
#              object$finals[missFinals[which(matchUnique)],indexBin] <- uniqueHaps[matchHap[k],]
#            }  
#        }
        ## Now that we've filled in missing values, should better be able to fill out haps

        binFounders[[i]][, j] <- apply(fdp[, touse, drop=F], 1, function(x) match(paste(x, collapse=""), do.call(paste0, enlist(uniqueHaps))))
        binFinals[[i]][,j] <- apply(fin[, touse, drop=F], 1, function(x) match(paste(x, collapse=""), do.call(paste0, enlist(uniqueHaps))))
      } else {
        binFounders[[i]][,j] <- object$founders[, indexBin]
        binFinals[[i]][,j] <- object$finals[,indexBin]
      }
    }
    binInfo[[i]] <- data.frame(MarkerName=colnames(object$founders)[object$lg$groups==i], bin=binInfo[[i]], group=i, binMarkerName=paste("C", i, "B", binInfo[[i]], sep=""))
  }
  output$founders <- do.call("cbind", binFounders)
  output$finals <- do.call("cbind", binFinals)
  output$bins <- do.call("rbind", binInfo)
   
  ## Collapse down the linkage groups
  output$lg$groups <- rep(1:n.chr, unlist(lapply(binInfo, function(x) max(x$bin))))
  names(output$lg$groups) <- colnames(output$finals)

  class(output) <- c("binmpcross", "mpcross")
  return(output)
}  
