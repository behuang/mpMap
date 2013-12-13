\name{mpbinhap}
\alias{mpbinhap}
\title{Get Bin founder Haplotypes}
\usage{
  mpbinhap(object, chr, bin)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{The chromosome where the bin is located}

  \item{bin}{The bin to get hfounder haplotypes for}
}
\value{
  If n.founders is the number of founders and n.mrk is the
  number of markers in the bin,the outpuut is a n.founders
  x n.mrk matrix where each row is a founder haplotype.
}
\description{
  Returns the founder haplotypes for a given bin
}

