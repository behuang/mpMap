\name{mpbin}
\alias{mpbin}
\title{Group Clustered Markers into Bins}
\usage{
  mpbin(object, chr, parameter,
    method = c("single", "complete", "average", "ward"))
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{Linkage group to bin markers on. If omitted,
  the default is all chromosomes}

  \item{parameter}{Threshold value used for clustering
  algorithm}

  \item{method}{Clustering algorithm to use.}
}
\value{
  Returned object is of the class \code{mpcross} with the
  additional component \code{bin}. If n.mrk is the number
  of markers genotyped and n.chr is the number of linkage
  groups, this is a list with the following components
  \item{bins}{n.mrk x 2 matrix labelling each marker by
  assigned linkage group and bin number in that group}
  \item{methods}{A data frame with n.chr columns and 2
  rows, giving information on the method and parameter used
  to group each linkage group} \item{summary}{A matrix with
  n.chr columns and 2 rows. The first row gives the number
  of markers in the linkage group , the second row gives
  the number of bins}
}
\description{
  Places tightly linked markers into bins, grouping is
  based on recombination fractions.
}
\examples{
map <- sim.map(len=rep(10,2), n.mar=30, eq.spacing=FALSE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
fgmat<-matrix(c(1,1,0,0), ncol=30*2, nrow=4)
colnames(fgmat)<-unlist(lapply(map, names))
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1, fg=fgmat)
dat.rf <- mpestrf(sim.dat)
dat.lg<-mpgroup(dat.rf, 2)
dat.bin<-mpbin(dat.lg, method="complete", parameter=0)
}

