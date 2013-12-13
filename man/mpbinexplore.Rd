\name{mpbinexplore}
\alias{mpbinexplore}
\title{Compute Binning Summary Statistics}
\usage{
  mpbinexplore(object, chr,
    methods = c("single", "complete", "average", "ward"),
    parameters = seq(0, 0.025, 0.005),
    breakpoints = c(0, 0.01, 0.05, 0.1, 0.5))
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{Chromosomes to bin}

  \item{methods}{Methods to be used for clustering. Options
  are "single", "complete", "average" and "ward".}

  \item{parameters}{Vector of parameter values to use}

  \item{breakpoints}{Intervals to use in summarising the
  distribution of maximium bin recombination fractions}
}
\value{
  Returned object is a list with the following four
  components: \item{chr}{The chromosome the analysis was
  carried out on} \item{n.markers}{The number of markers in
  that linkage group} \item{summary}{A matrix with a row
  for each element of \code{parameters}, and a column for
  each element of \code{methods}. The values in the methods
  column give the number of bins that result from grouping
  based on the parameter value in each row.}
  \item{distribution}{A list with a component for each
  parameter value specified. Each component in the list is
  a table summarising the distirbution of maximum intra bin
  recombination fractions. The numbers in the cells give
  the number of bins which have a maximum internal bin
  recombination fraction falling in the intervals dictated
  by breakpoints.}
}
\description{
  Summarises the effect of grouping markers into bins under
  different parameters. Markers are clustered using the
  hclust package. Estimates for all pairwise recombination
  fractions must be present within linkage groups.
}
\examples{
map <- sim.map(len=rep(10,2), n.mar=30, eq.spacing=FALSE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
fgmat<-matrix(c(1,1,0,0), ncol=30*2, nrow=4)
colnames(fgmat)<-unlist(lapply(map, names))
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1, fg=fgmat)
dat.rf <- mpestrf(sim.dat)
dat.lg<-mpgroup(dat.rf, 2)
exp1<- mpbinexplore(dat.lg,chr=1, parameters=c(0, 0.005, 0.01))
exp2<- mpbinexplore(dat.lg,chr=2, methods=c("single", "complete"), breakpoints=c(0, 0.05, 0.2, 0.5))
}

