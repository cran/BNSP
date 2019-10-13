\name{mvrm}
\alias{mvrm}

\title{Bayesian semiparametric analysis of multivariate continuous responses, with variable selection}

\description{Implements an MCMC algorithm for posterior sampling based on a semiparametric model for
continuous multivariate responses and additive models for the mean and variance  functions. The model
utilizes spike-slab priors for variable selection and regularization. See `Details' section for a full description
of the model.}

\usage{
mvrm(formula,data=list(),sweeps,burn=0,thin=1,seed,StorageDir,
c.betaPrior="IG(0.5,0.5*n*p)", pi.muPrior="Beta(1,1)", c.alphaPrior="IG(1.1,1.1)",
sigmaPrior="HN(2)", pi.sigmaPrior="Beta(1,1)", mu.RPrior="N(0,1)",
sigma.RPrior="HN(1)",corr.Model=c("common",nClust=1),DP.concPrior="Gamma(5,2)",
tuneAlpha,tuneSigma2,tuneCb,tuneCa,tuneR,tuneSigma2R,tau,FT=1,...)
}

\arguments{
\item{formula}{a formula defining the responses and the covariates in the mean and variance models e.g. \code{y1 | y2 ~ x | z} or
for smooth effects \code{y1 | y2 ~ sm(x) | sm(z)}. The package uses the extended formula notation, where the responses are
defined on the left of \code{~} and the mean and variance models on the right.}

\item{data}{a data frame.}

\item{sweeps}{total number of posterior samples, including those discarded in burn-in period
(see argument \code{burn}) and those discarded by the thinning process
(see argument \code{thin}).}

\item{burn}{length of burn-in period.}

\item{thin}{thinning parameter.}

\item{seed}{optional seed for the random generator.}

\item{StorageDir}{a required directory to store files with the posterior samples of models parameters.}

\item{c.betaPrior}{The inverse Gamma prior of \eqn{c_{\beta}}.
The default is "IG(0.5,0.5*n*p)", that is, an inverse Gamma with parameters \eqn{1/2} and
\eqn{np/2}, where \eqn{n} is the number of sampling units and \eqn{p} is the length of the response vector.}

\item{pi.muPrior}{The Beta prior of \eqn{\pi_{\mu}}.
The default is "Beta(1,1)". It can be of dimension \eqn{1}, of dimension \eqn{K} (the number of effects that enter the mean model),
or of dimension \eqn{pK}}

\item{c.alphaPrior}{The inverse Gamma prior of \eqn{c_{\alpha}}.
The default is "IG(1.1,1.1)". Half-normal priors for \eqn{\sqrt{c_{\alpha}}} are also available, declared using "HN(a)", where "a"
is a positive number. It can be of dimension \eqn{1} or \eqn{p} (the length of the multivariate response).}

\item{sigmaPrior}{The prior of \eqn{\sigma}.
The default is "HN(2)", a half-normal prior for \eqn{\sigma} with variance equal to two,
\eqn{\sigma \sim N(0,2) I[\sigma>0]}. Inverse Gamma priors for \eqn{\sigma^2} are also available,
declared using "IG(a,b)". It can be of dimension \eqn{1} or \eqn{p} (the length of the multivariate response).}

\item{pi.sigmaPrior}{The Beta prior of \eqn{\pi_{\sigma}}.
The default is "Beta(1,1)". It can be of dimension \eqn{1}, of dimension \eqn{Q} (the number of effects that enter the variance model), or of dimension \eqn{pQ}}

\item{mu.RPrior}{The normal prior for \eqn{\mu_{R}}.
The default is the standard normal distribution.}

\item{sigma.RPrior}{The half normal prior for \eqn{\sigma_{R}}.
The default is the half normal distribution with variance one.}

\item{corr.Model}{Specifies the model for the correlation matrix \eqn{R}.
The three choices supported are "common", that specifies a common correlations model,
"groupC", that specifies a grouped correlations model, and "groupV", that specifies
a grouped variables model. When the model chosen is either "groupC" or "groupV", the upper
limit on the number of clusters can also be specified, using corr.Model = c("groupC", nClust = d)
or corr.Model = c("groupV", nClust = p). If the number of clusters is left unspecified,
for the "groupV" model, it is taken to be \eqn{p}, the number of responses.
For the "groupC" model, it is taken to be \eqn{d = p(p-1)/2}, the number of free elements in the correlation matrix.}

\item{DP.concPrior}{The Gamma prior for the Dirichlet process concentration parameter.}

\item{tuneAlpha}{Starting value of the tuning parameter for sampling regression coefficients of the variance model
\eqn{\alpha}. Defaults at 5.}

\item{tuneSigma2}{Starting value of the tuning parameter for sampling variances \eqn{\sigma^2_j}.
Defaults at 1.}

\item{tuneCb}{Starting value of the tuning parameter for sampling \eqn{c_{\beta}}.
Defaults at 10.}

\item{tuneCa}{Starting value of the tuning parameter for sampling \eqn{c_{\alpha}}.
Defaults at 1.}

\item{tuneR}{Starting value of the tuning parameter for sampling correlation matrices.
Defaults at \eqn{100(p+2)}.}

\item{tuneSigma2R}{Starting value of the tuning parameter for sampling \eqn{\sigma_{R}^2}.
Defaults at 1.}

\item{tau}{The \eqn{tau} of the shadow prior. Defaults at 0.01.}

\item{FT}{Binary indicator. If set equal to 1, the Fisher's z transform of the correlations is modelled, otherwise if set equal to 0, the untransformed correlations are modelled.}

\item{...}{Other options that will be ignored.}
}

\details{Function \code{mmvrm} returns samples from the posterior distributions of the parameters of
a regression model with normally distributed multivariate responses and mean and variance functions modeled
in terms of covariates. For instance, in the presence of two responses (\eqn{y_1, y_2}) and two covariates
in the mean model (\eqn{u_1, u_2}) and two in the variance model (\eqn{w_1, w_2}), we may choose to fit
\deqn{
\mu_u = \beta_0 + \beta_1 u_1 + f_{\mu}(u_2),}
\deqn{
\log(\sigma^2_W) = \alpha_0 + \alpha_1 w_1 + f_{\sigma}(w_2),}
parametrically modelling the effects of \eqn{u_1} and \eqn{w_1} and non-parametrically modelling
the effects of \eqn{u_2} and \eqn{w_2}. Smooth functions, such as \eqn{f_{\mu}}
and \eqn{f_{\sigma}}, are represented by basis function expansion,
\deqn{
f_{\mu}(u_2) = \sum_{j} \beta_{j} \phi_{j}(u_2),
}
\deqn{
f_{\sigma}(w_2) = \sum_{j} \alpha_{j} \phi_{j}(w_2),
}
where \eqn{\phi} are the basis functions and \eqn{\beta} and \eqn{\alpha} are regression coefficients.

The variance model can equivalently be expressed as
\deqn{
\sigma^2_W = \exp(\alpha_0) \exp(\alpha_1 w_1 + f_{\sigma}(w_2)) = \sigma^2 \exp(\alpha_1 w_1 + f_{\sigma}(w_2)),}
where \eqn{\sigma^2 = \exp(\alpha_0)}. This is the parameterization that we adopt in this implementation.

Positive prior probability that the regression coefficients in the mean model are exactly
zero is achieved by defining binary variables \eqn{\gamma} that take value \eqn{\gamma=1}
if the associated coefficient \eqn{\beta \neq 0} and \eqn{\gamma = 0} if \eqn{\beta = 0}.
Indicators \eqn{\delta} that take value \eqn{\delta=1} if the associated coefficient
\eqn{\alpha \neq 0} and \eqn{\delta = 0} if \eqn{\alpha = 0} for the variance function
are defined analogously. We note that all coefficients in the mean and variance functions are
subject to selection except the intercepts, \eqn{\beta_0} and \eqn{\alpha_0}.

\emph{Prior specification:}

For the vector of non-zero regression coefficients \eqn{\beta_{\gamma}} we specify a g-prior
\deqn{
\beta_{\gamma} | c_{\beta}, \sigma^2, \gamma, \alpha, \delta \sim N(0,c_{\beta} \sigma^2
(\tilde{X}_{\gamma}^{\top} \tilde{X}_{\gamma} )^{-1}).
}
where \eqn{\tilde{X}} is a scaled version of design matrix \eqn{X} of the mean model.

For the vector of non-zero regression coefficients \eqn{\alpha_{\delta}} we specify a normal prior
\deqn{
\alpha_{\delta} | c_{\alpha}, \delta \sim N(0,c_{\alpha} I).
}

Independent priors are specified for the indicators variables \eqn{\gamma} and \eqn{\delta} as
\eqn{P(\gamma = 1 | \pi_{\mu}) = \pi_{\mu}} and
\eqn{P(\delta = 1 | \pi_{\sigma}) = \pi_{\sigma}}.
Further, Beta priors are specified for \eqn{\pi_{\mu}} and \eqn{\pi_{\sigma}}
\deqn{
\pi_{\mu} \sim Beta(c_{\mu},d_{\mu}), \pi_{\sigma} \sim Beta(c_{\sigma},d_{\sigma}).
}
We note that blocks of regression coefficients associated with distinct covariate effects
have their own probability of selection (\eqn{\pi_{\mu}} or \eqn{\pi_{\sigma}})
and this probability has its own prior distribution.

Further, we specify inverse Gamma priors for \eqn{c_{\beta}} and \eqn{c_{\alpha}}
\deqn{
c_{\beta} \sim IG(a_{\beta},b_{\beta}),
c_{\alpha} \sim IG(a_{\alpha},b_{\alpha})}

For \eqn{\sigma^2} we consider inverse Gamma and half-normal priors
\deqn{
\sigma^2 \sim IG(a_{\sigma},b_{\sigma}), |\sigma| \sim N(0,\phi^2_{\sigma}).
}

Lastly, for the elements of the correlation matrix, we specify normal distributions with mean \eqn{\mu_R} and variance \eqn{\sigma^2_R},
with the priors on these two parameters being normal and half-normal, respectively. This is the common correlations model. Further, 
the grouped correlations model can be specified. It considers a mixture of normal distributions for the means \eqn{\mu_R}. The grouped correlations model can also be specified. It clusters the variables instead of the correlations.  
}

\value{Function \code{mmvrm} returns the following:
\item{call}{the matched call.}
\item{formula}{model formula.}
\item{seed}{the seed that was used (in case replication of the results is needed).}
\item{data}{the dataset}
\item{X}{the mean model design matrix.}
\item{Z}{the variance model design matrix.}
\item{LG}{the length of the vector of indicators \eqn{\gamma}.}
\item{LD}{the length of the vector of indicators \eqn{\delta}.}
\item{mcpar}{the MCMC parameters: length of burn in period, total number of samples, thinning period.}
\item{nSamples}{total number of posterior samples}
\item{DIR}{the storage directory}
Further, function \code{mvrm} creates files where the posterior samples are written.
These files are (with all file names preceded by `BNSP.'):
\item{alpha.txt}{contains samples from the posterior of vector \eqn{\alpha}.
Rows represent posterior samples and columns represent the
regression coefficient, and they are in the same order as the columns of design matrix Z.}
\item{beta.txt}{contains samples from the posterior of vector \eqn{\beta}.
Rows represent posterior samples and columns represent the
regression coefficients, and they are in the same order as the columns of
design matrix X.}
\item{gamma.txt}{contains samples from the posterior of the vector of the indicators
\eqn{\gamma}. Rows represent posterior samples and columns represent the
indicator variables, and they are in the same order as the columns of design matrix X.}
\item{delta.txt}{contains samples from the posterior of the vector of the indicators
\eqn{\delta}. Rows represent posterior samples and columns represent the
indicator variables, and they are in the same order as the columns of design matrix Z.}
\item{sigma2.txt}{contains samples from the posterior of the error variance \eqn{\sigma^2} of each response.}
\item{cbeta.txt}{contains samples from the posterior of \eqn{c_{\beta}}.}
\item{calpha.txt}{contains samples from the posterior of \eqn{c_{\alpha}}.}
\item{R.txt}{contains samples from the posterior of the correlation matrix \eqn{R}.}
\item{theta.txt}{contains samples from the posterior of \eqn{\theta} of the shadow prior (probably not needed).}
\item{muR.txt}{contains samples from the posterior of \eqn{\mu_R}.}
\item{sigma2R.txt}{contains samples from the posterior of \eqn{\sigma^2_{R}}.}
\item{deviance.txt}{contains the deviance, minus twice the log likelihood evaluated at the sampled values of the parameters.}
In addition to the above, for models that cluster the correlations we have
\item{compAlloc.txt}{The cluster at which the correlations were allocated, \eqn{\lambda_{kl}}. These are integers from zero to the specified number of clusters minus one.}
\item{nmembers.txt}{The numbers of correlations assigned to each cluster.}
\item{DPconc.txt}{Contains samples from the posterior of the Dirichlet process concentration parameter.}
In addition to the above, for models that cluster the variables we have
\item{compAllocV.txt}{The cluster at which the variables were allocated, \eqn{\lambda_{k}}. These are integers from zero to the specified number of clusters minus one.}
\item{nmembersV.txt}{The numbers of variables assigned to each cluster.}
}

\references{
Papageorgiou, G. and Marshall, B.C. (2019). Bayesian semiparametric analysis of multivariate continuous responses, with variable selection. arXiv.

Papageorgiou, G. (2018). BNSP: an R Package for fitting Bayesian semiparametric regression models and variable selection. The R Journal, 10(2):526-548.

Chan, D., Kohn, R., Nott, D., & Kirby, C. (2006). Locally adaptive semiparametric estimation of the mean and variance functions in regression models. Journal of Computational and Graphical Statistics, 15(4), 915-936.
}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
# Fit a mean/variance regression model on the cps71 dataset from package np
require(np)
require(ggplot2)
data(cps71)
model <- logwage ~ sm(age,k=30,bs="rd") | sm(age,k=30,bs="rd")
DIR<-getwd()
\dontrun{m1 <- mvrm(formula=model,data=cps71,sweeps=10000,burn=5000,thin=2, seed=1,StorageDir=DIR)
#Print information and summarize the model
print(m1)
summary(m1)
#Summarize and plot one parameter of interest
alpha<-mvrm2mcmc(m1,"alpha")
summary(alpha)
plot(alpha)
#Obtain a plot of a term in the mean model
wagePlotOptions<-list(geom_point(data=cps71,aes(x=age,y=logwage)))
plot(x=m1,model="mean",term="sm(age)",plotOptions=wagePlotOptions)
plot(m1)
#Obtain predictions for new values of the predictor "age"
predict(m1,data.frame(age=c(21,65)),interval="credible")

# Fit a bivariate mean/variance model on the marks dataset from package ggm
# two responses: marks mechanics and vectors, and one covariate: marks on algebra 
model2 <- mechanics | vectors ~ sm(algebra,k=5) | sm(algebra,k=3)
m2 <- mvrm(formula=model2, data=marks, sweeps = 100000, burn = 50000, 
                       thin = 2, seed = 1, StorageDir = DIR)
plot(m2)
}
}

\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
