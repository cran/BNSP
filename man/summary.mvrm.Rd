\name{summary.mvrm}
\alias{summary.mvrm}

\title{Summary of an mvrm fit}

\description{
Provides basic information from an mvrm fit.
}

\usage{
\method{summary}{mvrm}(object, nModels = 5, digits = 5, ...)
}

\arguments{
\item{object}{an object of class "mvrm", usually a result of a call to \code{mvrm}.}
\item{nModels}{integer number of models with the highest posterior probability to be displayed.}
\item{digits}{the number of significant digits to use when printing.}
\item{...}{other arguments.}
}

\details{Use this function to summarize mvrm fits.}

\value{The functions provides a detailed description of the specified model and priors. 
In addition, the function provides information about the Markov chain ran (length, burn-in, thinning) and the
folder where the files with posterior samples are stored. Lastly, the function provides the mean
posterior and null deviance and the mean/variance models visited most often during posterior sampling.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}
