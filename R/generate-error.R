generate_error <- function(geno, error.prob, founderErrors = TRUE)
{
	#Copy input data to return value. 
	obsgeno <- geno
	
	n.founders <- nrow(geno$founders)
	n.mrk <- ncol(geno$founders)
	n.finals <- nrow(geno$finals) 

	founderGenotypes <- lapply(1:ncol(geno$founders), function(x) unique(geno$founders[,x]))
	
	if(founderErrors)
	{
		#Matrix telling us where the founder errors are
		founderErrors <- matrix(data=sample(c(TRUE,FALSE), n.founders*n.mrk, replace=TRUE, prob=c(error.prob, 1-error.prob)), nrow=n.founders, ncol=n.mrk)
		for(markerIndex in 1:n.mrk)
		{
			#The founders which are going to have the wrong value, for this marker
			founderErrorsThisMarker <- which(founderErrors[,markerIndex])
			for(founderIndex in founderErrorsThisMarker)
			{
				#Put in one of the other alleles (chosen uniformly)
				obsgeno$founders[founderIndex,markerIndex] <- sample(setdiff(founderGenotypes[[markerIndex]], obsgeno$founders[founderIndex, markerIndex]), 1)
			}
		}
	}
	
	#Matrix telling us where the final errors are.
	finalErrors <- matrix(data=sample(c(TRUE,FALSE), n.finals*n.mrk, replace=TRUE, prob=c(error.prob, 1-error.prob)), nrow=n.finals, ncol=n.mrk)
	for(markerIndex in 1:n.mrk)
	{
		#Possible observed values at this marker
		for (finalValue in founderGenotypes[[markerIndex]]) 
		{
			#The other possible observed values (the ones errors will be changed to)
			otherValues <- setdiff(founderGenotypes[[markerIndex]], finalValue)
			#For the values which are going to be errors and previously had observation finalValue, choose uniformly from the alternatives
			obsgeno$finals[finalErrors[,markerIndex] & obsgeno$finals[,markerIndex] == finalValue, markerIndex] <- sample(otherValues, sum(finalErrors[,markerIndex] & obsgeno$finals[,markerIndex] == finalValue), replace=TRUE)
		}
	}
	return(obsgeno)
}

