\name{mpview}
\alias{mpview}
\title{Review possible genotyping errors and missing values}
\usage{
  mpview(object, chr, bin, method = "partials")
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{chromosome where errors are suspected.}

  \item{bin}{bin where errors are suspected}

  \item{method}{Type of classification of interest,
  "errors" referes to genotyping errors, "partials" refers
  to missing values}
}
\value{
  A list with components \code{chr}, \code{bin},
  \code{lines}, \code{markers},\code{haplotypes} and
  \code{linetypes} \item{chr}{The chromosome where the
  suspected errors are located (\code{chr})} \item{bin}{The
  bin where the the suspected errors are located
  (\code{bin}) } \item{lines}{The line indexes where
  suspected errors have been found} \item{markers}{The
  indexes of markers which belong to bin number \code{bin}
  on chromosome \code{chr}} \item{haplotypes}{A n.founder x
  n.mrk matrix, where n.founders is the number of founders,
  and n.mrk is the number of markers in \code{bin}. The
  matrix gives the genotype of each founder at each marker.
  Each row is treated as a haplotype for each founder.}
  \item{linetypes}{A n.line x n.mrk matrix, where n.line is
  the number of lines in \code{lines} amd n.mrk is number
  of markers in \code{bin}. The matrix gives the genotype
  of each line at each marker. Row i gives the genotyping
  information for the line with index given by the ith
  entry in \code{lines}. Column j reprents the marker with
  index given by the jth entry in \code{markers}}
}
\description{
  This function gives more comprehensive information about
  possible genotyping errors and missing values flagged by
  \code{mperrorcheck}. The entries in the \code{errors}
  matrix of the output give possible locations to review in
  more detail. The function returns the bin founder
  haplotypes, lines where missing values/ errors are
  suspected and their phenotypes on the specified
  chromosome and bin number.
}

