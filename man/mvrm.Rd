\name{mvrm}
\alias{mvrm}

\title{Bayesian models for normally distributed responses and semiparametric 
       mean and variance regression models}

\description{MCMC for normally distributed responses with additive model for the 
mean and variance functions achieved via spike-slab prior for variable selection. 
See `Details' section for a full description of the model.}

\usage{
mvrm(formula,data,sweeps,burn=0,thin=1,seed,StorageDir,
c.betaPrior="IG(0.5,0.5*n)",c.alphaPrior="IG(1.1,1.1)",
pi.muPrior="Beta(1,1)",pi.sigmaPrior="Beta(1,1)",sigmaPrior="HN(2)",...)
}

\arguments{
\item{formula}{a formula defining the response and the covariates in the mean and variance models e.g. \code{y ~ x | z} or 
for smooth effects \code{y ~ sm(x) | sm(z)}. The package uses the extended formula notation, where on the right of \code{~} 
we define two models: on the right of \code{|} is the mean model and on the left is the variance model.}

\item{data}{a required data frame.}

\item{sweeps}{total number of posterior samples, including those discarded in burn-in period 
(see argument \code{burn}) and those discarded by the thinning process 
(see argument \code{thin}).}

\item{burn}{length of burn-in period.}

\item{thin}{thinning parameter.}

\item{seed}{optional seed for the random generator.}

\item{StorageDir}{a required directory to store files with the posterior samples of models parameters.}

\item{c.betaPrior}{The parameters of the inverse Gamma prior of \eqn{c_{\beta}}. 
The default is "IG(0.5,0.5*n)", that is, an inverse Gamma with parameters \eqn{1/2} and 
\eqn{n/2}, where \eqn{n} is the sample size.}

\item{c.alphaPrior}{The parameters of the inverse Gamma prior of \eqn{c_{\alpha}}. 
The default is "IG(1.1,1.1)".}

\item{pi.muPrior}{The parameters of the Beta prior of \eqn{\pi_{\mu}}. 
The default is "Beta(1,1)".}

\item{pi.sigmaPrior}{The parameters of the Beta prior of \eqn{\pi_{\sigma}}. 
The default is "Beta(1,1)".}

\item{sigmaPrior}{The prior of \eqn{\sigma^2}. 
The default is "HN(2)", a half-normal with variance equal to two, that is \eqn{|\sigma| \sim N(0,2)}.
Inverse Gamma prior is also available.}

\item{...}{Other options that will be ignored.}
}

\details{Function \code{mvrm} returns samples from the posterior distributions of the parameters of 
a regression model with normally distributed responses and mean and variance functions modeled 
in terms of covariates. For instance, in the presence of two covariates in the mean model
(\eqn{u_1, u_2}) and two in the variance model (\eqn{w_1, w_2}), we may choose to fit
\deqn{
\mu_u = \beta_0 + \beta_1 u_1 + f_{\mu}(u_2),}
\deqn{
\log(\sigma^2_W) = \alpha_0 + \alpha_1 w_1 + f_{\sigma}(w_2),}
parametrically modelling the effects of \eqn{u_1} and \eqn{w_1} and non-parametrically
the effects of \eqn{u_2} and \eqn{w_2}. Smooth functions, such as \eqn{f_{\mu}}
and \eqn{f_{\sigma}}, are represented by basis function expansion. For instance
\deqn{
f_{\mu}(u_2) = \sum_{j} \beta_{j} \phi_{j}(u_2),
}
where \eqn{\phi} are the basis functions and \eqn{\beta} are the associated coefficients. 

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

Lastly, for \eqn{\sigma^2} we consider inverse Gamma and half-normal priors
\deqn{
\sigma^2 \sim IG(a_{\sigma},b_{\sigma}), |\sigma| \sim N(0,\phi^2_{\sigma}).
}
}

\value{Function \code{mvrm} returns the following:
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
\item{sigma2.txt}{contains samples from the posterior of the error variance \eqn{\sigma^2}.}
\item{cbeta.txt}{contains samples from the posterior of \eqn{c_{\beta}}.}
\item{calpha.txt}{contains samples from the posterior of \eqn{c_{\alpha}}.}
}

\references{
Chan, D., Kohn, R., Nott, D., & Kirby, C. (2006). Locally adaptive semiparametric estimation of the
mean and variance functions in regression models. Journal of Computational and Graphical Statistics,
15(4), 915-936.

Papageorgiou, G. (2018). BNSP: an R Package for fitting Bayesian regression models with semiparametric mean and variance functions. arXiv:1804.10939 [stat.OT]
}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
# Fit a mean/variance regression model on the cps71 dataset from package np
require(np)
require(ggplot2)
data(cps71)
model <- logwage ~ sm(age,nknots=30,bs="rd") | sm(age,nknots=30,bs="rd")
DIR<-getwd()
m1<-mvrm(formula=model,data=cps71,sweeps=10,burn=5,thin=1,seed=1,StorageDir=DIR)
\donttest{m2 <- mvrm(formula=model,data=cps71,sweeps=10000,burn=5000,thin=2, seed=1,StorageDir=DIR)  
#Print information and summarize the model
print(m2)
summary(m2)
#Summarize and plot one parameter of interest
alpha<-mvrm2mcmc(m2,"alpha")
summary(alpha)
plot(alpha)
#Obtain a plot of a term in the mean model
wagePlotOptions<-list(geom_point(data=cps71,aes(x=age,y=logwage)))
plot(x=m2,model="mean",term="sm(age)",plotOptions=wagePlotOptions)
#Obtain predictions for new values of the predictor "age"
predict(m2,data.frame(age=c(21,65)),interval="credible")
}
}

\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
