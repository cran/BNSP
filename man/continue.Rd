\name{continue}
\alias{continue}

\title{Continues the sampler from where it stopped}

\description{
Allows the user to continue the sampler from the state it stopped in the previous call to \code{mvrm}.
}

\usage{continue(object, sweeps, discard = FALSE,...)}

\arguments{
\item{object}{An object of class "mvrm", usually a result of a call to \code{mvrm}.}
\item{sweeps}{The number of additional sweeps, maintaining the same thinning interval as specified in the original call to \code{mvrm}.}
\item{discard}{If set to true, the previous samples are discarded.}
\item{...}{other arguments.}
}

\details{The function allows the sampler to continue from the state it last stopped.}

\value{The function returns an object of class \code{mvrm}.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}
