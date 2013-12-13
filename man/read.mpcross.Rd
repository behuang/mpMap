\name{read.mpcross}
\alias{read.mpcross}
\title{Construct mpcross objects from datafiles}
\usage{
  read.mpcross(founderfile, finalfile, pedfile, mapfile,
    phenofile)
}
\arguments{
  \item{founderfile}{File containing founder genotypes -
  should have 1+(number of founders) rows. The first row
  contains the marker names - first space left blank. The
  first column contains the founder name. Additional rows
  contain observed marker data for all founders.}

  \item{finalfile}{File containing final genotypes - should
  have 1+(number of lines) rows. The first row contains the
  marker names - first space left blank. The first column
  contains line names. Additional rows contain observed
  marker data for all lines.}

  \item{pedfile}{File containing pedigree data - should
  have three or four columns - first three columns indicate
  id, mother and father; fourth column is a flag for
  whether the lines was observed.}

  \item{mapfile}{File containing linkage map - should
  contain three columns - one for marker names, one for
  chromosome assignment and one for map position in cM}

  \item{phenofile}{File containing phenotypic data - should
  contain one column for each phenotype, with rows indexing
  lines.}
}
\description{
  Generate an mpcross object by reading in components from
  files - requires founders marker data, finals marker
  data, pedigree, and IDs for all lines. Marker map and
  phenotypic data are optional.
}
\seealso{
  \code{\link[mpMap]{sim.mpcross}},
  \code{\link[mpMap]{mpcross}}
}

