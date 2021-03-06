% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processes.R
\name{wiener.process}
\alias{wiener.process}
\title{Create a Wiener Process}
\usage{
wiener.process(dispersion = 1)
}
\arguments{
\item{dispersion}{the dispersion parameter of the Wiener process.}
}
\value{
a function hanlde in the form of \code{X(tObs,n)} which generates \code{n} independent trajectories observed at \code{tObs}.
}
\description{
Create a Wiener Process
}
\examples{
X <- wiener.process()
X(regular.grid(50),25)
}
