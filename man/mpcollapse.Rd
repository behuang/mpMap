\name{mpcollapse}
\alias{mpcollapse}
\title{Collapse an mpcross object into a binned object}
\usage{
  mpcollapse(object, method = c("rf", "cor"))
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{method}{Choice of whether to bin based on
  recombination fraction or correlation}
}
\value{
  A binned object of class \code{binmpcross}
}
\description{
  Collapses an mpcross object into bins. Markers are
  recoded as binned markers (BM) where the genotype at a BM
  is determined from the haplotypes observed in the
  founders within the bin. Information on which markers
  belong to each bin is retained within the output object.
}
\note{
  Markers are assigned to a bin if there is zero
  recombination between any of the markers within the bin
  (or correlation of 1). Missing values at markers are
  imputed based on the assumption that the true haplotype
  within a bin matches one of the founder haplotypes.
  Patterns matching more than one or none of the founder
  haplotypes will be coded as NA in the returned object.
  Bins are collapsed to "bin markers" after missing values
  are imputed, with each unique founder haplotype being
  assigned an allele and progeny matching those founder
  haplotypes assigned the same allele.
}
\seealso{
  \code{link{mpexpand}}
}

