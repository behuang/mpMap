\name{plot.mpprob}
\alias{plot.mpprob}
\title{Plot summary of founder probabilities and haplotype blocks}
\usage{
  \method{plot}{mpprob} (x, chr, locations,
    compositionPercent = TRUE, lines = TRUE, nlines,
    compositionTrace = TRUE, compositionTraceArgs = list(),
    compositionTracePlotArgs = list(),
    linesPlotArgs = list(), ...)
}
\arguments{
  \item{x}{Object of class \code{mpprob}}

  \item{chr}{Chromosomes to plot. Default is all}

  \item{compositionPercent}{Flag for whether to plot the
  percent alleles inherited from each founder}

  \item{nlines}{Number of most recombinant lines to plot}

  \item{lines}{Flag for whether to plot haplotype
  reconstructions}

  \item{compositionTrace}{Flag for ??}

  \item{compositionTraceArgs}{List of ??}

  \item{compositionTracePlotArgs}{List of ??}

  \item{...}{Additional arguments to plot function}
}
\value{
  Barplot of the percentage of each founder on each
  chromosome; individual heatmaps of which chunks of each
  chromosome are inherited from each founder.
}
\description{
  Plot the percentage of each chromosome inherited from
  each founder
}
\examples{
sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl")
plot(mpp.dat)
}
\seealso{
  \code{\link[mpMap]{mpprob}},
  \code{\link[mpMap]{summary.mpprob}},
  \code{\link[Heatplus]{heatmap_2}}
}

