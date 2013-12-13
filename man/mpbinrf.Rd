\name{mpbinrf}
\alias{mpbinrf}
\title{Return Bin Recombination Fractions}
\usage{
  mpbinrf(object, chr, bin)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{The chromosome where the bin is located}

  \item{bin}{The bin of interest}
}
\value{
  A matrix of recombination fraction estimates between
  markers in the specified bin
}
\description{
  Returns the matrix of recombination fraction estimates
  between markers in a specified linkage group and bin
}
\examples{
map <- sim.map(len=rep(10,2), n.mar=30, eq.spacing=FALSE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
fgmat<-matrix(c(1,1,0,0), ncol=30*2, nrow=4)
colnames(fgmat)<-unlist(lapply(map, names))
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1, fg=fgmat)
dat.rf <- mpestrf(sim.dat)
dat.lg<-mpgroup(dat.rf, 2)
dat.bin<-mpbin(dat.lg, method="complete", parameter=0.05)
bin.rf <- mpbinrf(dat.bin, chr=1, bin=2)
}

