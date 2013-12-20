\name{subset.mpcross}
\alias{subset.mpcross}
\title{Subset mpcross object}
\usage{
\method{subset}{mpcross}(x, groups = NULL, chr = NULL, markers = NULL,
  lines = NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{mpcross}}

  \item{chr}{Selected chromosomes TO KEEP. Default is all}

  \item{groups}{selected groups to KEEP. Default is all}

  \item{markers}{Selected markers TO KEEP. Default is all}

  \item{lines}{Selected lines TO KEEP. Default is all}

  \item{...}{Additional arguments}
}
\value{
The original object with chromosomes/lines/markers removed
which are not listed in the arguments.
}
\description{
Reduces an mpcross object down to a specified set of
chromosomes, markers and/or lines
}
\note{
Chromosomes can be input either as the character names of
chromosomes or the index of the chromosomes in the map.
Markers can be input as character names or the index in the
matrix x$finals. Lines can be input as either character
values (matching the rownames of x$finals) or indices of
rows in that matrix. Groups must be input as numbers
corresponding to the specific groups to keep (the groups
present in the object may not form a consecutive set of
numbers). Note that only one of chr, groups and markers can
be selected.
}
\examples{
map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=map, pedigree=ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
sim.dat
red.dat <- subset(sim.dat, chr=1, lines=1:50)
red.dat
}
\seealso{
\code{\link[mpMap]{mpcross.object}}
}

