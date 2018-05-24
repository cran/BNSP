\name{print.mvrm}
\alias{print.mvrm}

\title{Prints an mvrm fit}

\description{
Provides basic information from an mvrm fit.
}

\usage{
\method{print}{mvrm}(x, digits = 5, ...)
}

\arguments{
\item{x}{an object of class "mvrm", usually a result of a call to \code{mvrm}.}
\item{digits}{the number of significant digits to use when printing.}
\item{...}{other arguments.}
}

\details{The function prints information about mvrm fits.}

\value{The function provides a matched call, the number of posterior samples obtained and marginal inclusion probabilities of the terms 
in the mean and variance models.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}
