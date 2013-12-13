\name{mpheat}
\alias{mpheat}
\title{Recombination Fraction Heatmap}
\usage{
  mpheat(obj, chr, bin, cols = c("white", "steelblue"),
    limits = c(0, 0.05))
}
\arguments{
  \item{obj}{Object of class \code{mpcross}}

  \item{chr}{The chromosome where the bin is located}

  \item{bin}{The bin of interest}

  \item{cols}{The colours to use in the heatmap}

  \item{limits}{The numerical range to use}
}
\value{
  None
}
\description{
  Generate a heatmap of recombination fraction for a given
  bin and chromosome
}

