\name{runCommands}
\alias{runCommands}
\title{Apply interactive changes}
\usage{
  runCommands(mpcross, commands)
}
\arguments{
  \item{mpcross}{The object to which the changes in
  \code{commands} were made}

  \item{commands}{The set of commands to be run, as
  returned by plotGTK in the commands variable of the
  output environment}
}
\value{
  Returns \code{NULL}.
}
\description{
  Apply changes made interactively to a heatmap using
  \code{\link[mpMap]{plotGTK}} to an \code{mpcross} object
}
\details{
  This experimental function is part of a new workflow for
  grouping and ordering markers. It allows the user to take
  commands recorded by the \code{\link[mpMap]{plotGTK}}
  function and apply them to the corresponding
  \code{mpcross} object. Allowable commands include
  Deleting markers Deleting linkage groups Merging linkage
  groups Constructing linkage subgroups which will be
  ordered consecutively Constructing reordering groups
  which will be reordered independent of the rest of the
  linkage group

  See \code{\link[mpMap]{plotGTK}} for an example .
}
\seealso{
  \code{\link[mpMap]{plotGTK}}
}

