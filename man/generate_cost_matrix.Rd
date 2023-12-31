% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{generate_cost_matrix}
\alias{generate_cost_matrix}
\title{Generate a simple cost matrix}
\usage{
generate_cost_matrix(
  charset,
  match = 0L,
  mismatch = 1L,
  gap = NULL,
  gap_open = NULL
)
}
\arguments{
\item{charset}{A string representing all possible characters in both query and target sequences (e.g. "ACGT")}

\item{match}{The cost of a match}

\item{mismatch}{The cost of a mismatch}

\item{gap}{The cost of a gap or NULL if this parameter will be set later.}

\item{gap_open}{The cost of a gap opening or NULL. If this parameter is set, gap must also be set.}
}
\value{
A cost matrix
}
\description{
Generate a cost matrix for use with the \code{search} method
}
\examples{
generate_cost_matrix("ACGT", match = 0, mismatch = 1)
}
