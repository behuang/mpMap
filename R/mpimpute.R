#' Impute missing values for an mpcross object
#' 
#' Use multi-point founder probabilities to fill in the most likely value for a missing genotype
#' @export
#' @param object Object of class \code{mpcross}
#' @param what Whether to impute founders, finals, or both
#' @param threshold Threshold for probability to call an allele
#' @param calls Form of the finals output - discrete (genotype calls if above the threshold), or continuous (expectation of allele). For multiallelic markers no value will be imputed if continuous output is selected.
#' @param \dots Additional parameters to be passed into mpprob function
#' @return An mpcross object with a new component $missfinals and $missfounders for the original set of genotypes. The components $finals and $founders replaced by the imputed values.
#' @seealso \code{\link[mpMap]{mpprob}}, \code{\link[mpMap]{mpcross}}
#' @references Huang BE, Raghavan C, Mauleon R, Broman KW, Leung H (under review) Imputation of low-coverage genotyping-by-sequencing in multi-parental crosses.  

mpimpute <- function(object, what=c("both", "founders", "finals"), threshold=.5, calls=c("discrete", "continuous"), ...)
{
  output <- object
  output$missfinals <- output$finals
  output$missfounders <- output$founders
  if (missing(calls)) calls <- "discrete"
  nmrk <- ncol(object$finals)
  if (missing(what)) what <- "both"
  ## First impute founders if necessary
  if (what!="finals") {
    nmissfou <- apply(object$founders, 2, function(x) sum(is.na(x)))
    mpp <- mpprob(object, program="qtl", threshold=threshold, ...)
    object <- mpp
    est <- do.call("cbind", mpp$est)
    
    for (i in 1:nmrk) {
      missfou <- which(is.na(object$founders[,i]))
      for (j in missfou)
      {
        tab <- table(object$finals[which(est[,i]==j),i])
        if (length(tab)>0)
          object$founders[j,i] <- as.numeric(names(tab)[which.max(tab)])
      }
    }
  }
  
  ## Then do the finals
  if (what!="founders") {
    
    if (!inherits(object, "mpprob")) stop("Must compute founder probabilities first")
    
    chr <- names(object$map)
    
    n.founders <- nrow(object$founders)
    output <- object
    output$prob <- NULL
    output$estfnd <- NULL
    
    probs <- do.call("cbind", object$prob)
    ifmat <- ipmat <- matrix(nrow=nrow(probs), ncol=ncol(probs)/n.founders)
    founder <- as.vector(object$founders)
    
    nall <- apply(object$founders, 2, function(x) length(table(x)))
    biall <- which(nall==2)
    multiall <- which(nall>2) 
    uniall <- which(nall==1)
    
    #  if (length(biall)+length(multiall)!=ncol(object$founders)) stop("Monomorphic markers included, please remove before imputation\n")
    
    missfx <- function(x) {
    tmp <- by(x, fd, sum)
    if (max(tmp)>threshold) return(as.numeric(names(which.max(tmp)))) else return(NA)}
    
    for (jj in seq(1, ncol(probs), n.founders))
    {
      index <- (jj-1)/n.founders+1
      fd <- founder[jj:(jj+n.founders-1)]
      
      if (index%in%multiall){
        alltab <- apply(probs[,jj:(jj+n.founders-1)], 1, missfx)
        ifmat[,index] <- alltab  
      } 
      else if (index %in% uniall)
        ifmat[,index] <- ipmat[, index] <- founder[jj]
      else {
        fvec2 <- fvec <- founder[jj:(jj+n.founders-1)]
        fvec[fvec==min(fvec)] <- 0
        fvec[fvec==max(fvec)] <- 1
        tmp <- probs[, jj:(jj + n.founders - 1)] %*% fvec
        yn <- as.numeric(tmp>threshold)
        ifmat[, index] <- yn*max(fvec2)+(1-yn)*min(fvec2)
        ipmat[, index] <- tmp
      }
    }
    rownames(ifmat) <- rownames(ipmat) <- rownames(object$finals)
    colnames(ifmat) <- colnames(ipmat) <- colnames(object$finals)
    output$missfinals <- output$finals
    
    if (calls=="discrete")
      output$finals <- ifmat
    else
      output$finals <- ipmat
    attr(output, "imputethreshold") <- threshold
  }
    
  class(output) <- "mpcross"
  return(output)
}
