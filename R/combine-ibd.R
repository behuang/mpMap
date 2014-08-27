combine_ibd <- function(obsgeno)
{
	geno <- list()
	nmrk <- ncol(obsgeno$founders)/2
	geno$founders <- obsgeno$founders[,1:nmrk, drop=FALSE]
	geno$finals <- obsgeno$finals[,1:nmrk, drop=FALSE]
	geno$finals[obsgeno$finals[,1:nmrk]!=obsgeno$finals[,nmrk+1:nmrk]] <- NA
	return(geno)
}
