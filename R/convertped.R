#Convert a (possibly character-based) pedigree across to a numeric-based pedigree
convertped <- function(pedigree)
{
	#Rename columns
 	if(any(colnames(pedigree)[1:4] != c("id", "Male", "Female", "Observed"))) 
	{
		warning("Column names of pedigree should be 'id', 'Male', 'Female' and 'Observed'. Renaming...")
		colnames(pedigree)[1:4] <- c("id", "Male", "Female", "Observed")
	}
	#If there are not four columns, then there must be five and the fifth one must be named Design
	if(ncol(pedigree) != 4)
	{
		if(ncol(pedigree) != 5 || colnames(pedigree)[5] != "Design")
		{
			warning("If a pedigree has five columns the fifth must be named 'Design'. Renaming...")
			colnames(pedigree)[5] <- "Design"
		}
	}
	lineNames <- NULL
	#If pedigree has rownames which are not just 1:nrow(pedigree), keep them
	if(!is.null(rownames(pedigree)) && any(rownames(pedigree) != as.character(1:nrow(pedigree))))
	{
		lineNames <- rownames(pedigree)
	}
	#Otherwise if the first column is a character and not just 1:nrow(pedigree), keep that
	else if(is.character(pedigree[,"id"]) && any(pedigree[,"id"] != as.character(1:nrow(pedigree))))
	{
		lineNames <- pedigree[,"id"]
	}
	#convert pedigree to character, and then numeric
	ped <- apply(pedigree, 2, as.character)
	maleID <- match(ped[,"Male"], ped[,"id"])
	femaleID <- match(ped[,"Female"], ped[,"id"])

	#Any parents not in the pedigree are unknown
	maleID[is.na(maleID)] <- 0
	femaleID[is.na(femaleID)] <- 0

	newPedigree <- data.frame(id = 1:nrow(ped), Male = maleID, Female = femaleID, Observed = pedigree[,"Observed"])
	if(ncol(pedigree) == 5) newPedigree <- cbind(newPedigree, Design = pedigree$Design, stringsAsFactors = FALSE)

	#Check that pedigree is sorted (parents before children)
	sorted <- (newPedigree[,"id"] > newPedigree[,"Female"] & newPedigree[,"id"] > newPedigree[,"Male"])
	if(!any(sorted))
	{
		stop("Pedigree is not sorted, parents should appear before children")
	}
	#Check number of founders
	nFounders <- sum(pedigree[,"Male"]==0 & pedigree[,"Female"]==0)
	if (nFounders != 2 && nFounders !=4 && nFounders !=8)
	{
		stop("Cannot process a number of founders which is not 2, 4 or 8")
	}
	#Keep line names, if there were any
	if(!is.null(lineNames))
	{
		rownames(newPedigree) <- lineNames
	}
	else rownames(newPedigree) <- paste0("L", 1:nrow(newPedigree))
  	return(newPedigree)
}

