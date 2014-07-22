#' Plot summary of founder probabilities and haplotype blocks
#' 
#' Plot the percentage of each chromosome inherited from each founder
#' @S3method plot mpprob
#' @method plot mpprob
#' @param x Object of class \code{mpprob}
#' @param chr Chromosomes to plot. Default is all
#' @param compositionPercent Flag for whether to plot the percent alleles inherited from each founder
#' @param nlines Number of most recombinant lines to plot
#' @param lines Flag for whether to plot haplotype reconstructions
#' @param compositionTrace Flag for ??
#' @param compositionTraceArgs List of ??
#' @param compositionTracePlotArgs List of ??
#' @param ... Additional arguments to plot function
#' @return Barplot of the percentage of each founder on each chromosome; individual heatmaps of which chunks of each chromosome are inherited from each founder.
#' @seealso \code{\link[mpMap]{mpprob}}, \code{\link[mpMap]{summary.mpprob}}, \code{\link[Heatplus]{heatmap_2}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' plot(mpp.dat)

plot.mpprob <-
function(x, chr, locations, compositionPercent = TRUE, lines=TRUE, nlines, compositionTrace = TRUE, compositionTraceArgs = list(), compositionTracePlotArgs = list(), linesPlotArgs = list(), ...)
{
  require(RColorBrewer)
  require(graphics)
  require(grDevices)
	if(!is.list(compositionTraceArgs))
	{
		stop("Input compositionTraceArgs must be a list")
	}
  if (missing(chr)) chr <- 1:length(x$map)
  n.founders <- nrow(x$founders)

  cts1 <- lapply(x$estfnd, function(y) {
	z <- factor(as.vector(y), levels=1:n.founders)
	return(round(table(z)/prod(dim(y))*100, 2)) })
  cts <- do.call("cbind", cts1)

  ## first plot is the stacked barplot of founder probabilities
  cts1 <- cbind(cts, rep(NA, nrow(cts))) 
  cts1 <- rbind(cts1, 100 - colSums(cts1))
    
  colours <- brewer.pal(n.founders, "Spectral") 
  #colours <- rainbow(n.founders)

	if(compositionPercent)
	{
		barplot(cts1, col=c(colours, "white"), main="Founder %age by Chromosome", xlab="Chromosome", legend.text = c(rownames(cts1)[1:n.founders], "Unknown"))
	}

	sum <- vector(length=nrow(x$finals))
  if (lines) {
  	for (i in chr)
		{
			nrec <- apply(x$estfnd[[i]], 1, function(x) return(sum(diff(x[!is.na(x)])!=0)))
      if(missing(nlines)) nlines <- nrow(x$finals)
			## if lines is not missing, need to select out the ones with the most
			## recombination events
			ord <- cbind(1:length(nrec), nrec)
			ord <- ord[order(ord[,2], decreasing=TRUE), ]    
			sel <- ord[1:nlines, 1]

			if(missing(locations)) mat <- x$estfnd[[i]][sort(sel),]
			else mat <- x$estfnd[[i]][sort(sel),locations]

			args <- list(t(mat), col=colours, main=names(x$map)[[i]], Rowv=NA, Colv=NA, scale="none", xlab="Lines", ylab="Markers")
			#If we have overrides for xlab, ylab or main, don't use the default values. So remove them. 
			removeIndices <- na.omit(match(intersect(c("xlab", "ylab", "main"), names(linesPlotArgs)), names(args)))
			if(length(removeIndices) > 0) args <- args[-removeIndices]
			do.call(heatmap, c(args, linesPlotArgs))
		}
  }
  
	if(compositionTrace)
	{
		relevantFounders <- 1:n.founders
		if("relevantFounders" %in% names(compositionTraceArgs))
		{
			relevantFounders <- compositionTraceArgs[["relevantFounders"]]
			compositionTraceArgs["relevantFounders"] <- NULL
			if(length(relevantFounders) < 1) stop("Invalid value for parameter relevantFounders")
		}
		percentMissing <- apply(x$finals, 1, function(x) sum(is.na(x))/length(x))
		probabilitiesMap <- attr(x$prob, "map")
		
		#Toss lines with more than 50% missing data
		for (i in 1:length(x$map)) x$estfnd[[i]] <- x$estfnd[[i]][which(percentMissing<.5),]
		x$finals <- x$finals[which(percentMissing<.5),]
		
		probabilities <- vector(mode="list", length=n.founders)
		#Get out a seperate data set for each founder
		for(founder in relevantFounders)
		{
			for (i in 1:length(probabilitiesMap)) 
			{
				newBit <- x$prob[[i]][,seq(founder, ncol(x$prob[[i]]), n.founders)]
				colnames(newBit) <- names(probabilitiesMap[[i]])
				probabilities[[founder]] <- c(probabilities[[founder]], apply(newBit, 2, function(x) mean(x, na.rm=TRUE)))
			}
		}
		
		#Concat up the map locations
		cm <- 0
		for (i in 1:length(probabilitiesMap)) cm <- c(cm, probabilitiesMap[[i]]+max(cm))
		cm <- cm[2:(length(cm))]

		joinedProbabilities <- do.call(c, probabilities[relevantFounders])
		joinedProbabilities <- data.frame(prob = joinedProbabilities)
		joinedProbabilities$founders <- rep(rownames(x$founders)[relevantFounders], each = length(unlist(probabilitiesMap)))
		joinedProbabilities$cm <- rep(cm, length(relevantFounders))
		joinedProbabilities$founders <- as.factor(joinedProbabilities$founders)
		
		chrEndPoints <- unlist(lapply(probabilitiesMap, max))
		chrEndPoints <- cumsum(chrEndPoints)

		library(lattice)

		colors <- brewer.pal(length(relevantFounders), "Spectral")
		xyargs <- list(prob~cm|founders, data=joinedProbabilities, ylab="Founder Probability", xlab="Chromosome", as.table=T, layout=c(1,length(relevantFounders)), panel=
		function(x,y) {panel.xyplot(x,y, col=colors[panel.number()], type="l", lwd=2) 
			 if(length(probabilitiesMap) > 1) for (i in 1:(length(probabilitiesMap)-1)) panel.abline(v=chrEndPoints[i], lwd=2, col="tomato")})
		tmp <- do.call(xyplot, c(xyargs, compositionTraceArgs))
		tmp$x.scales$at <- c(0, cumsum(unlist(lapply(probabilitiesMap, max))[1:(length(probabilitiesMap)-1)]))+unlist(lapply(probabilitiesMap, max))/2
		tmp$x.scales$labels <- names(probabilitiesMap)
		do.call(plot, c(list(tmp), compositionTracePlotArgs))
	}
}

