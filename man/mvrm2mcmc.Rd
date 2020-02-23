\name{mvrm2mcmc}
\alias{mvrm2mcmc}

\title{Convert posterior samples from function \code{mvrm} into an object of class `mcmc'}

\description{Reads in files where the posterior samples were written and creates an object of class 
`mcmc' so that functions like \code{summary} and \code{plot} from package \code{coda} can be used}

\usage{mvrm2mcmc(mvrmObj, labels)}

\arguments{
\item{mvrmObj}{An object of class `mvrm' as created by a call to function \code{mvrm}.}

\item{labels}{The labels of the files to be read in. These can be one or more of: "alpha", "beta", "gamma", "delta", "sigma2", "cbeta", "calpha", "R", "muR", "sigma2R", "nmembers", "nmembersV", "compAlloc", "compAllocV", and "DPconc" and they correspond to the parameters of the model that a call to functions \code{mvrm} fits. In addition, "deviance" can be read in. If left unspecified, all files are read in.}
}

\value{An object of class `mcmc' that holds the samples from the posterior of the 
selected parameter.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}
