\name{mpskeleton}
\alias{mpskeleton}
\title{Create a skeleton mpcross object}
\usage{
  mpskeleton(object)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}
}
\value{
  Returned object is of the class \code{mpcross} which is a
  subset of the original object.
}
\description{
  This function is to reduce a binned mpcross object to its
  minimal form. The function returns a subset of the
  original object, keeping one marker per bin.
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
dat.skeleton<-mpskeleton(dat.bin)
summary(dat.bin)
summary(dat.skeleton)
}

