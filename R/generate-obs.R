generate_obs <- function(geno, map, full.prob, fg, transpos, transval, founderld)
{
	if(founderld && !is.null(fg))
	{
		stop("Input fg cannot be specified if founderld == TRUE")
	}
	if(is.null(fg))
	{
		return(generate_obs_without_fg(geno, map, full.prob, transpos, transval, founderld))
	}
	else
	{
		return(generate_obs_with_fg(geno, map, full.prob, fg, transpos, transval))
	}
}
generate_obs_without_fg <- function(geno, map, full.prob, transpos, transval, founderld)
{
	#Retain everything from geno into obsgeno (the return value)
	obsgeno <- geno
	n.founders <- nrow(geno$founders)
	nMarkers <- ncol(geno$founders)/2
	#which markers are goign to be biallelic?
	biallelic <- sample(c(FALSE, TRUE), nMarkers, replace=TRUE, prob=c(full.prob, 1-full.prob))
	whichBiallelic <- which(biallelic)

	#Probability distribution for the number of '1' alleles at any given marker. Obviously they can't all be 1, so there's a 0 at the end. 
	domprob <- c(rep(1/(n.founders-1), n.founders-1), 0)
	dom <- matrix(data=0, nrow=n.founders, ncol=nMarkers)

	##founderld == TRUE means founder genotypes generated wrt rf
	if (founderld) 
	{
		#create genotypes for founders at first biallelic marker
		ndom <- sample(1:n.founders, size=1, prob=domprob)
		dom[sample(1:n.founders, ndom), whichBiallelic[1]] <- 1

		#Subset marker data to just the biallelic markers
		allmrk <- unlist(lapply(map, names))
		chr <- rep(1:length(map), unlist(lapply(map, length)))
		pos <- unlist(map)
		
		dommrk <- allmrk[whichBiallelic]
		domchr <- chr[whichBiallelic]
		dompos <- pos[whichBiallelic]

		#Recombination fractions between consecutive dominant markers
		recfr <- haldaneX2R(abs(diff(dompos)))
		recfr[which(diff(domchr)!=0)] <- 0.5

		#For each founder
		for (i in 1:n.founders) 
		{
			#now at each dominant marker position generate a uniform according to recombination fraction with "previous" dominant marker
			runi <- runif(length(recfr))
			#Rec tells us whether we should switch the 0/1 values. This happens with probability proportional to recombination fraction
			rec <- (runi < recfr)
			#we want to respect the 0/1 value already simulated for the first biallelic marker
			rec <- c(dom[i, whichBiallelic[1]], rec)
			dom[i, whichBiallelic] <- cumsum(rec) %% 2
		}
		biallelic <- c(biallelic,biallelic)
		dom <- cbind(dom, dom)
		obsgeno$founders[,biallelic] <- dom[,biallelic]
	} 
	else 
	{
		#In this case the number of 1 values is independent for different markers
		ndom <- sample(1:n.founders, sum(biallelic), replace=TRUE, prob=domprob)
		#Generate 0/1 values for each dominant marker
		for (i in 1:sum(biallelic)) 
		{
			dom[sample(1:n.founders, ndom[i]),whichBiallelic[i]] <- 1
		}
		#Duplicate dominance values
		dom <- cbind(dom, dom)

		biallelic <- c(biallelic, biallelic)
		#Copy dominant marker data into obsgeno (the return value)
		obsgeno$founders[,which(biallelic)] <- dom[,which(biallelic)]
	 
		#set special founders for translocation region if markers are biallelic - They should all have the same biallelic pattern
		if (transval>0) 
		{
			for (i in transpos) 
			{ 
				if (biallelic[i])
				{
					obsgeno$founders[,i]  <- obsgeno$founders[,i+nMarkers] <- c(rep(0, transval-1), 1, rep(0, n.founders - transval))
				}
			}
		}
	} 

	#Put in zero values where appropriate
	for (i in 1:n.founders)
	{
		markersWithZeroAlleleThisFounder <- (obsgeno$founders[i,]==0)
		obsgeno$finals[,markersWithZeroAlleleThisFounder][geno$finals[,markersWithZeroAlleleThisFounder]==i] <- 0
	}
	#Everything that is biallelic and has allele > 0 gets a value of 1
	obsgeno$finals[,biallelic][(obsgeno$finals[,biallelic]>0)] <- 1
	return(obsgeno)
}
generate_obs_with_fg <- function(geno, map, full.prob, fg, transpos, transval)
{
	#Retain everything from geno into the return value (obsgeno), except for founders and finals. 
	obsgeno <- geno
	obsgeno$founders <- NULL
	
	n.founders <- nrow(geno$founders)
	#At this point the genetic map INCLUDES the QTL. But the fg input DOESN'T. A copy of fg is made, which includes extra columns for the QTL markers. The fg values for these markers are copied from the closest linked marker. 
		
	allmrk <- unlist(lapply(map, names))
	chr <- rep(1:length(map), unlist(lapply(map, length)))
	pos <- unlist(map)

	#Recombination fractions between consectutive markers
	recfr <- haldaneX2R(abs(diff(pos)))
	recfr[which(diff(chr)!=0)] <- 0.5

	qtlpos <- grep("QTL", allmrk)

	#Copy of fg which contains extra markers
	fg2 <- matrix(nrow=nrow(fg), ncol=length(allmrk))
	colnames(fg2) <- allmrk
	rownames(fg2) <- rownames(fg)
	indices <- match(colnames(fg), colnames(fg2))
	#Copy across everything in fg to fg2
	fg2[, na.omit(indices)] <- as.matrix(fg)[,which(!is.na(indices))]

	#Put in fg for the most tightly linked marker, for each QTL marker
	for (i in 1:length(qtlpos))
	{
		if(qtlpos[i] == 1)
		{
			warning("QTL found as first marker")
			fg2[,1]  <- fg2[,2]
		}
		else if(qtlpos[i] == length(allmrk))
		{
			warning("QTL found as last marker")
			fg2[,length(allmrk)]  <- fg2[,length(allmrk)-1]
		}
		else if(recfr[qtlpos[i]] == 0.5 && recfr[qtlpos[i]-1] == 0.5)
		{
			warning("A QTL marker was the only one on a chromosome. Setting founders = 1:nFounders")
			fg2[,qtlpos[i]] <- 1:n.founders
		}
		else if(recfr[qtlpos[i]] <= recfr[qtlpos[i]-1])
		{
			fg2[,qtlpos[i]] <- fg2[,qtlpos[i]+1]
		}
		else
		{
			fg2[,qtlpos[i]] <- fg2[,qtlpos[i]-1]
		}
	}
	#fg2 is the founder genotypes
	obsgeno$founders <- cbind(fg2, fg2)
	#Overlay founder genotypes on finals. 
	for (i in 1:ncol(obsgeno$finals))
	{
		obsgeno$finals[,i] <- obsgeno$founders[,i][obsgeno$finals[,i]]
	}
	return(obsgeno)
}