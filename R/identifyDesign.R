identifyDesignPair <- function(pair, nFounders)
{
	char1 <- nchar(pair[1])
	char2 <- nchar(pair[2])
	motherName <- pair[3]
	fatherName <- pair[4]
	maxG <- log(nFounders, base=2)
	headStr <- paste(nFounders, "wayG", maxG, sep="")
	#Are we begining a selfing or AIC line?
	if(pair[1] == headStr && pair[2] == headStr)
	{
		if(motherName == fatherName)
		{
			return(paste(nFounders, "wayG", maxG, "self1", sep=""))
		}
		return(paste(nFounders, "wayG", maxG, "aic1", sep=""))
	}
	#Are we doing part of the initial cross?
	headStr <- paste(nFounders, "wayG", sep="")
	if(char1 == 6 && char2 == 6 && substr(pair[1], 1, 5) == headStr && substr(pair[2], 1, 5) == headStr)
	{
		n <- substr(pair[1], 6, 6)
		if(n >= maxG)
		{
			stop("Too many generations of initial crossing")
		}
		return(paste(nFounders, "wayG", as.integer(n)+1, sep=""))
	}
	#Straight selfing line?
	headStr <- paste0(nFounders, "wayG", maxG, "self")
	if(motherName == fatherName && pair[1] == pair[2] && substr(pair[1], 1, 10) == headStr && char1 == 11)
	{
		generation <- as.integer(substr(pair[1], 11, 11))
		return(paste0(nFounders, "wayG", maxG, "self", generation+1))
	}
	headStr <- paste0(nFounders, "wayG", maxG, "aic")
	if(substr(pair[1], 1, 9) == headStr && pair[1] == pair[2])
	{
		#AIC line?
		if(motherName != fatherName && char1 == 10)
		{
			generation <- as.integer(substr(pair[1], 10, 10))
			return(paste0(nFounders, "wayG", maxG, "aic", generation+1))
		}
		#Initial selfing of an AIC line?
		if(motherName == fatherName && char1 == 10)
		{
			return(paste0(substr(pair[1], 1, 10), "self1"))
		}
		if(motherName == fatherName && substr(pair[1], 11, 14) == "self" && char1 == 15)
		{
			generation <- as.integer(substr(pair[1], 15, 15))
			return(paste0(substr(pair[1], 1, 14), generation+1))
		}
	}
	#isG1 <- regexec(paste(nFounders, "wayG[0-9]", sep=""), pair[1])
	#isG2 <- regexec(paste(nFounders, "wayG[0-9]", sep=""), pair[2])
	return(NA)
}

#' Identify subtype of design in pedigree
#'
#' For a pedigree containing multiple variations on MAGIC designs, identify
#' the subtype of design, e.g., MAGIC+
#' @export
#' @param pedigree Pedigree to categorize
#' @return An additional vector identifying designs of individuals within the pedigree

identifyDesign <- function(pedigree)
{
	if(inherits(pedigree, "mpcross")) pedigree <- pedigree$pedigree
	designCol <- rep(NA, nrow(pedigree))
	isInitial <- pedigree[,2] == 0 & pedigree[,3] == 0
	nFounders <- sum(isInitial)
	if(nFounders != 4 && nFounders != 8) stop("nFounders must be 4 or 8")
	designCol[which(isInitial)] <- paste(nFounders, "wayG0", sep="")
	names(designCol) <- rownames(pedigree)
	addedDesigns <- which(isInitial)
	#Standardise the parents founders as NA (Otherwise, 0 and "0" have a different behaviour when the occur in pedigree[,2] of designCol[pedigree[,2]] - 0 shortens the result, but "0" returns an NA)
	pedigree[1:nFounders, 2:3] <- NA
	while(length(addedDesigns) > 0)
	{
		canAddDesigns <- is.na(designCol) & !is.na(designCol[pedigree[,2]]) & !is.na(designCol[pedigree[,3]])
		#Now we're dealing with only the subset where canAddDesigns is true
		motherDesign <- designCol[pedigree[which(canAddDesigns),2]]
		fatherDesign <- designCol[pedigree[which(canAddDesigns), 3]]
		motherName <- pedigree[pedigree[which(canAddDesigns),2],1]
		fatherName <- pedigree[pedigree[which(canAddDesigns),3],1]
		designCombinations <- cbind(motherDesign, fatherDesign, motherName, fatherName)
		newDesign <- apply(designCombinations, 1, identifyDesignPair, nFounders = nFounders)
		#And we want to further deal with the subset for which newDesign is not NA
		identifiedDesigns <- which(canAddDesigns)[which(!is.na(newDesign))]
		designCol[identifiedDesigns] <- newDesign
		addedDesigns <- identifiedDesigns
	}
	return(designCol)
}
