### Note - this assumes all markers are generated as 0/1
### Does not generate errors in founders, just finals

generate_error2 <- function(geno, error.prob)
{
  obsgeno <- geno
  n.founders <- nrow(geno$founders)
  n.mrk <- ncol(geno$founders)/2
  n.finals <- nrow(geno$finals) 

#  fdr.err <- matrix(data=sample(c(TRUE,FALSE), n.founders*n.mrk, replace=TRUE, 
#	prob=c(error.prob, 1-error.prob)), nrow=n.founders, ncol=n.mrk)
  fin.err <- matrix(data=sample(c(TRUE,FALSE), n.finals*n.mrk, replace=TRUE,
	prob=c(error.prob, 1-error.prob)), nrow=n.finals, ncol=n.mrk)
  fin.err <- cbind(fin.err, fin.err)
#  fdr.err <- cbind(fdr.err, fdr.err)

#  obsgeno$founders[fdr.err==1] <- 1-obsgeno$founders[fdr.err==1]
  obsgeno$finals[fin.err==1] <- 1-obsgeno$finals[fin.err==1]

  # replacing IBD genotypes with something different
#  for (i in 1:n.founders) {
#     prob <- rep(error.prob/(n.founders-1), n.founders) 
#     prob[i] <- 0

#     obsgeno$founders[fdr.err==1 & obsgeno$founders==i] <- sample(1:n.founders, sum(fdr.err==1 & obsgeno$founders==i), replace=TRUE, prob=prob)
#     obsgeno$finals[fin.err==1 & obsgeno$finals==i] <- sample(1:n.founders, sum(fin.err==1 & obsgeno$finals==i), replace=TRUE, prob=prob)
#  }

  return(obsgeno)
}

