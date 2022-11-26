\name{predict.mvrm}
\alias{predict.mvrm}

\title{Model predictions}

\description{
Provides predictions and posterior credible/prediction intervals for given feature vectors.
}

\usage{
\method{predict}{mvrm}(object, newdata, interval = c("none", "credible", "prediction"), 
                       level = 0.95, nSamples = 100, ind.preds=FALSE, ...)
}

\arguments{
\item{object}{an object of class "mvrm", usually a result of a call to \code{mvrm}.}
\item{newdata}{data frame of feature vectors to obtain predictions for. If newdata is missing, the function will use
the feature vectors in the data frame used to fit the mvrm object.}
\item{interval}{type of interval calculation.}
\item{level}{tolerance level.}
\item{nSamples}{number of samples to obtain from the posterior predictive distribution (for each sweep of the MCMC). 
Only relevant for "prediction intervals".}
\item{ind.preds}{Binary indicator. If set to TRUE the function returnsadditionally the predictions per individual MCMC sample.}
\item{...}{other arguments.}
}

\details{The function returns predictions of new responses or the means of the responses for given feature vectors.
Predictions for new responses or the means of new responses are the same. However, the two differ in the associated 
level of uncertainty: response predictions are associated with wider (prediction) intervals than mean response
predictions. To obtain prediction intervals (for new responses) the function samples from the normal distributions
with means and variances as sampled during the MCMC run.}

\value{Predictions for given covariate/feature vectors.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}
