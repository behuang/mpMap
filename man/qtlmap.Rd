\name{qtlmap}
\alias{qtlmap}
\title{Select markers in a region around QTL}
\usage{
qtlmap(qtlpos, qtlchr, map, window = 10, qtlnam)
}
\arguments{
  \item{qtlpos}{Vector of QTL positions}

  \item{qtlchr}{Vector of QTL chromosomes (one for each
  position)}

  \item{map}{Linkage map to determine where QTL are
  relative to other markers}

  \item{window}{Number of cM to each side of QTL in which
  to find nearby markers}

  \item{qtlnam}{Optional vector of names for the QTL}
}
\value{
Returns a map selecting out regions +- window around the
QTL positions.
}
\description{
Outputs a list of markers in regions around QTL. Region is
defined by the window parameter of x cM to either side of
the QTL positions.
}
\examples{
sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl")
mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
qmap <- qtlmap(summary(mpq.dat)[,2], summary(mpq.dat)[,1], mpq.dat$map)
plotlink.map(qmap)
}
\seealso{
\code{\link[mpMap]{plotlink.map}}
}

