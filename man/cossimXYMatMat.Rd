% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosinedistance.R
\name{cossimXYMatMat}
\alias{cossimXYMatMat}
\title{Cosine similarity between a pair of matrices}
\usage{
cossimXYMatMat(x, y)
}
\arguments{
\item{x}{Matrix}

\item{y}{Matrix}
}
\value{
A similarity value for each row.
}
\description{
Matrices \code{x} and \code{y} must be of equal dimension. The i-th row in
\code{x} is compared to the i-th row in \code{y}.
}
