#Convert a (possibly character-based) pedigere across to a numeric-based pedigree
convertped <- function(pedigree, nFounders)
{
# 	if(any(colnames(pedigree)[1:4] != c("id", "Male", "Female", "Observed"))) 
#	{
#		warning("Column names of pedigree should be 'id', 'Male', 'Female' and 'Observed'. Renaming...")
#		colnames(pedigree)[1:4] <- c("id", "Male", "Female", "Observed")
#	}
	if(ncol(pedigree) != 4)
	{
		if(ncol(pedigree) != 5 || colnames(pedigree)[5] != "Design")
		{
			warning("If a pedigree has five columns the fifth must be named 'Design'. Renaming...")
			colnames(pedigree)[5] <- "Design"
		}
	}
	if(missing(nFounders)) nFounders <- sum(pedigree[,2]==0 & pedigree[,3]==0)
	if (nFounders != 2 && nFounders !=4 && nFounders !=8)
	{
		stop("Cannot process a number of founders which is not 4 or 8")
	}
	# convert pedigree to numeric, if it's not already
	ped <- apply(pedigree, 2, as.character)
	maleID <- match(ped[,2], ped[,1])
	femaleID <- match(ped[,3], ped[,1])

	# Any parents not in the pedigree are unknown
	maleID[is.na(maleID)] <- 0
	femaleID[is.na(femaleID)] <- 0

	newPedigree <- data.frame(id = 1:nrow(ped), Male = maleID, Female = femaleID, Observed = pedigree[,4])
	if(ncol(pedigree) == 5) newPedigree <- cbind(newPedigree, Design = pedigree$Design, stringsAsFactors = FALSE)

	sorted <- (newPedigree[,1] > newPedigree[,3] & newPedigree[,1] > newPedigree[,2])
	if(!any(sorted))
	{
		stop("Pedigree is not sorted, parents should appear before children")
	}
  	return(newPedigree)
}

