\name{clean.mpcross}
\alias{clean.mpcross}
\alias{mpsegrat}
\title{Check data format and compute summary statistics for genotypes}
\usage{
\method{clean}{mpcross}(object, ...)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{...}{Additional arguments}
}
\value{
\item{drop1}{List of markers which are monomorphic in
founders} \item{drop2}{List of markers which have missing
values in founders} \item{drop3}{List of markers with
alleles appearing in finals but not founders}
\item{alleles}{Number of alleles at each marker}
\item{missing}{Percent missing data at each marker}
\item{seg}{Matrix with one row for each marker and columns
for the marker name, the chisquare test for segregation
distortion, and the p-value of the test}
}
\description{
Given an object of class 'mpcross', the function checks
that the data is in the correct format, containing founder
and final genotypes, ids, and a pedigree. The number of
markers genotyped for both founders and finals should
coincide. The pedigree should be completely numeric.
Markers which are not polymorphic across the founders are
removed, as are markers which have missing values in the
founders.
}
\details{
Summary statistics for the genotypes are printed, included
the number of markers with varying levels of missing data,
with varying levels of segregation distortion, and with
different numbers of alleles.
}
\examples{
map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, qtl=matrix(data=c(1, 45, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE),seed=1)
dat.chk <- clean(sim.dat)
}
\seealso{
\code{\link[mpMap]{mpcross}}
}

