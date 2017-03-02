isInterestingFunnel <- function(row, founder1Index, founder2Index)
{
	#not interested in anything without a full funnel
	if(length(unique(row)) != length(row)) return(FALSE)
	whichFounder1 <- which(row == founder1Index)
	whichFounder2 <- which(row == founder2Index)
	if(whichFounder1 %% 2 == 0 && whichFounder2 + 1 == whichFounder1) return(TRUE)
	if(whichFounder2 %% 2 == 0 && whichFounder1 + 1 == whichFounder2) return(TRUE)
	return(FALSE)
}
testFunnelEffects <- function(prob)
{
	if(!inherits(prob, "mpprob"))
	{
		stop("Input object must have class \"mpprob\"")
	}
	if(!("prob" %in% names(prob)))
	{
		stop("Haplotype probabilities must be calculated before testing funnels effects")
	}
		n.founders <- nrow(prob$founders)
	
	map <- attr(prob$prob, "map")
	locations <- unlist(lapply(map, names))
	chromosome <- rep(names(map), times = unlist(lapply(map, length)))
	
	results <- vector(mode="list", length = n.founders)
	names(results) <- rownames(prob$founders)
	
	funnels <- getAllFunnels(prob)
	
	#Set up storage for lines for each set of funnels
	relevantLines <- vector(mode="list", length=n.founders)
	for(founder1Index in 1:n.founders)
	{
		relevantLines[[founder1Index]] <- vector(mode="list", length=n.founders)
	}
	
	for(founder1Index in 1:n.founders)
	{
		for(founder2Index in 1:n.founders)
		{
			if(founder1Index != founder2Index)
			{
				relevantLines[[founder1Index]][[founder2Index]] <- relevantLines[[founder2Index]][[founder1Index]] <- which(apply(funnels, 1, isInterestingFunnel, founder1Index = founder1Index, founder2Index = founder2Index))
			}
		}
	}
	w <- options("warn")
	options(warn=-1)
	for(founder1Index in 1:n.founders)
	{
		founder1Name <- rownames(prob$founders)[founder1Index]
		resultsThisFounder <- data.frame(Location = locations, value = rep(0, length(locations)))
		
		for(locationIndex in 1:length(locations))
		{
			column <- paste(locations[locationIndex], ", Founder ", founder1Index, sep="")
			otherColumns <- paste(locations[locationIndex], ", Founder ", (1:n.founders)[-founder1Index], sep="")
			founder1ProbThisLocation <- notFounder1ProbThisLocation <- vector(mode="numeric", length=n.founders)
			for(founder2Index in 1:n.founders)
			{
				if(founder1Index != founder2Index)
				{
					founder1ProbThisLocation[founder2Index] <- sum(prob$prob[[chromosome[locationIndex]]][relevantLines[[founder1Index]][[founder2Index]],column])
					notFounder1ProbThisLocation[founder2Index] <- sum(prob$prob[[chromosome[locationIndex]]][relevantLines[[founder1Index]][[founder2Index]],otherColumns])
				}
			}
			founder1ProbThisLocation <- founder1ProbThisLocation[-founder1Index]
			notFounder1ProbThisLocation <- notFounder1ProbThisLocation[-founder1Index]
			names(founder1ProbThisLocation) <- rownames(prob$founders)[-founder1Index]
			names(notFounder1ProbThisLocation) <- rownames(prob$founders)[-founder1Index]
			data <- cbind(founder1ProbThisLocation, notFounder1ProbThisLocation)
			founders <- as.factor(rownames(data))
			total <- stats::glm(data ~ founders, family=stats::binomial)
			probability <- drop1(total, ~founders, test="LRT")[[5]][2]
			resultsThisFounder$value[locationIndex] <- probability
		}
		results[[founder1Name]] <- resultsThisFounder
	}
	options(warn=w$warn)
	return(results)
}
