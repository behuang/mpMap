#' @export
#' @importFrom stats pchisq
mpsegrat <- function(object)
{
	if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")

	finals <- object$finals
	founders <- object$founders
	nmrk <- ncol(founders)
	chisq <- vector(length=nmrk)
	pval <- vector(length=nmrk)
	badmrk <- vector()
	for (i in 1:nmrk)
	{
		obs <- table(finals[,i])
		exp <- table(founders[,i])/nrow(founders)
		chisq[i] <- NA
		pval[i] <- NA
		possibleValues <- names(exp)
		if (!all(names(obs) %in% possibleValues))
		{
			badmrk <- c(badmrk, i) 
		}
		else 
		{
			#Cut out any hets from the observed values
			obs <- obs[intersect(names(obs), names(exp))]
			#If there are any actual genotypes (not hets) missing, put them in. 
			if (!all(names(exp) %in% names(obs))) 
			{
				#append zeros
				obs2 <- c(obs, rep(0, length(exp) - length(obs)))
				#put in right names
				names(obs2) <- c(names(obs), setdiff(names(exp), names(obs)))
			}
			else obs2 <- obs
			#Correct order
			obs2 <- obs2[names(exp)]
			
			exp <- exp * sum(obs)
			chisq[i] <- sum((obs2-exp)^2/exp)
			pval[i] <- 1-pchisq(chisq[i], length(obs)-1)
		}
	}
#	if (length(badmrk)>0)  cat("Markers ", badmrk, " had values appear in finals which are not in founders.\n They probably have genotyping errors.\n") 
# Move this into clean.mpcross
	df <- data.frame(MarkerName=colnames(founders), chisq, pval)
	df
}

