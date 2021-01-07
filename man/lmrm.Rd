\name{lmrm}
\alias{lmrm}

\title{Bayesian semiparametric modelling of covariance matrices for multivariate longitudinal data}

\description{Implements an MCMC algorithm for posterior sampling based on a semiparametric model for
continuous longitudinal multivariate responses. The overall model consists of 5 regression submodels and it
utilizes spike-slab priors for variable selection and function regularization. See `Details' section for a full description
of the model.}

\usage{
lmrm(formula, data = list(), id, time,
sweeps, burn = 0, thin = 1, seed, StorageDir,
c.betaPrior = "IG(0.5,0.5*n*p)", pi.muPrior = "Beta(1,1)", 
c.alphaPrior = "IG(1.1,1.1)", pi.phiPrior = "Beta(1,1)", c.psiPrior = "HN(2)",
sigmaPrior = "HN(2)", pi.sigmaPrior = "Beta(1,1)", 
corr.Model = c("common", nClust = 1), DP.concPrior = "Gamma(5,2)",
c.etaPrior = "IG(0.5,0.5*samp)", pi.nuPrior = "Beta(1,1)", 
pi.fiPrior = "Beta(1,1)", c.omegaPrior = "IG(1.1,1.1)", sigmaCorPrior = "HN(2)",
tuneCa, tuneSigma2, tuneCb, tuneAlpha, tuneSigma2R, tuneR, tuneCpsi, 
tuneCbCor, tuneOmega, tuneComega, tau, FT = 1,...)
}

\arguments{
\item{formula}{a formula defining the responses and the covariates in the 5 regression models e.g. \code{y1 | y2 ~ x | w | z | t | t} or
for smooth effects \code{y1 | y2 ~ sm(x) | sm(w) | sm(z) | sm(t) | sm(t)}. The package uses the extended formula notation, where the responses are
defined on the left of \code{~} and the models on the right.}

\item{data}{a data frame.}

\item{id}{identifiers of the individuals or other sampling units that are observed over time.}

\item{time}{a vector input that specifies the time of observation}

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
or of dimension \eqn{p K}}

\item{c.alphaPrior}{The inverse Gamma prior of \eqn{c_{\alpha}^2}.
The default is "IG(1.1,1.1)". Half-normal priors for \eqn{c_{\alpha}} are also available, declared using "HN(a)", where "a"
is a positive number. It can be of dimension \eqn{1} or \eqn{p} (the length of the multivariate response).}

\item{pi.phiPrior}{The Beta prior of \eqn{\pi_{\phi}}.
The default is "Beta(1,1)". It can be of dimension \eqn{1}, of dimension \eqn{B} (the number of effects that enter the dependence model), or of dimension \eqn{p^2 B}}

\item{c.psiPrior}{The prior of \eqn{c_{\psi}^2}.
The default is "HN(2)", a half-normal prior for \eqn{c_{\psi}} with variance equal to two,
\eqn{c_{\psi} \sim N(0,2) I[c_{\psi} > 0]}. Inverse Gamma priors for \eqn{c_{\psi}^2} are also available,
declared using "IG(a,b)". It can be of dimension \eqn{1} or \eqn{p^2} (the number of dependence models).}

\item{sigmaPrior}{The prior of \eqn{\sigma_k^2, k=1,\dots,p}.
The default is "HN(2)", a half-normal prior for \eqn{\sigma_k} with variance equal to two,
\eqn{\sigma_k \sim N(0,2) I[\sigma>0]}. Inverse Gamma priors for \eqn{\sigma_k^2} are also available,
declared using "IG(a,b)". It can be of dimension \eqn{1} or \eqn{p} (the length of the multivariate response).}

\item{pi.sigmaPrior}{The Beta prior of \eqn{\pi_{\sigma}}.
The default is "Beta(1,1)". It can be of dimension \eqn{1}, of dimension \eqn{L} (the number of effects that enter the variance model), or of dimension \eqn{pL}}

\item{corr.Model}{Specifies the model for the correlation matrices \eqn{R_t}.
The three choices supported are "common", that specifies a common correlations model,
"groupC", that specifies a grouped correlations model, and "groupV", that specifies
a grouped variables model. When the model chosen is either "groupC" or "groupV", the upper
limit on the number of clusters can also be specified, using corr.Model = c("groupC", nClust = d)
or corr.Model = c("groupV", nClust = p). If the number of clusters is left unspecified,
for the "groupV" model, it is taken to be \eqn{p}, the number of responses.
For the "groupC" model, it is taken to be \eqn{d = p(p-1)/2}, the number of free elements in the correlation matrices.}

\item{DP.concPrior}{The Gamma prior for the Dirichlet process concentration parameter.}

\item{c.etaPrior}{The inverse Gamma prior of \eqn{c_{\eta}}.
The default is "IG(0.5,0.5*samp)", that is, an inverse Gamma with parameters \eqn{1/2} and
\eqn{samp/2}, where \eqn{samp} is the number of correlations observed over time, that is $samp=M*d$ 
where $M$ is the number of unique observation time points and $d$ is the number of non-redundant elements 
of $R$.}

\item{pi.nuPrior}{The Beta prior of \eqn{\pi_{\nu}}. The default is "Beta(1,1)". It can be of dimension \eqn{1}.}

\item{pi.fiPrior}{The Beta prior of \eqn{\pi_{\varphi}}. The default is "Beta(1,1)". It can be of dimension \eqn{1}.}

\item{c.omegaPrior}{The prior of \eqn{c_{\omega}^2}.
The default is "HN(2)", a half-normal prior for \eqn{c_{\omega}} with variance equal to two,
\eqn{c_{\omega} \sim N(0,2) I[c_{\omega} > 0]}. Inverse Gamma priors for \eqn{c_{\omega}^2} are also available,
declared using "IG(a,b)". It can be of dimension \eqn{1}.}

\item{sigmaCorPrior}{The prior of \eqn{\sigma^2}.
The default is "HN(2)", a half-normal prior for \eqn{\sigma^2} with variance equal to two,
\eqn{\sigma \sim N(0,2) I[\sigma > 0]}. Inverse Gamma priors for \eqn{\sigma^2} are also available,
declared using "IG(a,b)". It can be of dimension \eqn{1}.}

\item{tuneCa}{Starting value of the tuning parameter for sampling \eqn{c_{\alpha k}, k=1,\dots,p}.
Defaults at a vector of $p$ ones. It could be of dimension \eqn{p}.}

\item{tuneSigma2}{Starting value of the tuning parameter for sampling \eqn{\sigma^2_k, k=1,\dots,p}.
Defaults at a vector of $p$ ones. It could be of dimension \eqn{p}.}

\item{tuneCb}{Starting value of the tuning parameter for sampling \eqn{c_{\beta}}.
Defaults at 100.}

\item{tuneAlpha}{Starting value of the tuning parameter for sampling regression coefficients of the variance models
\eqn{\alpha_k, k=1,\dots,p}. Defaults at a vector of 5s. It could be of dimension \eqn{L p}}

\item{tuneSigma2R}{Starting value of the tuning parameter for sampling \eqn{\sigma^2}.
Defaults at 1.}

\item{tuneR}{Starting value of the tuning parameter for sampling correlation matrices.
Defaults at \eqn{40*(p+2)^3}. Can be of dimension \eqn{1} or \eqn{M} is the number of unique observation time points.}

\item{tuneCpsi}{Starting value of the tuning parameter for sampling variances \eqn{c_{\psi}^2}.
Defaults at 5. Can be of dimension \eqn{1} or \eqn{p^2}}.

\item{tuneCbCor}{Starting value of the tuning parameter for sampling \eqn{c_{\eta}^2}.
Defaults at 10.}

\item{tuneOmega}{Starting value of the tuning parameter for sampling regression coefficients of the variance models
\eqn{\omega}. Defaults at 5.}

\item{tuneComega}{Starting value of the tuning parameter for sampling \eqn{c_{\omega}}.
Defaults at 1.}

\item{tau}{The \eqn{tau} of the shadow prior. Defaults at 0.01.}

\item{FT}{Binary indicator. If set equal to 1, the Fisher's z transform of the correlations is modelled, otherwise if set equal to 0, the untransformed correlations are modelled.}

\item{...}{Other options that will be ignored.}
}

\details{
	Function \code{lmrm} returns samples from the posterior distributions of the parameters of
    a regression model with normally distributed multivariate longitudinal responses. To describe the model, 
    let \eqn{Y_{ij} = (Y_{ij1},\dots,Y_{ijp})^{\top}} denote the vector of \eqn{p} responses
    observed on individual \eqn{i, i=1,\dots,n,} at time point \eqn{t_{ij}, j=1,\dots,n_i}. 
    The observational time points \eqn{t_{ij}} are allowed to be unequally spaced. 
    Further, let \eqn{u_{ij}} denote the covariate vector that is observed along with \eqn{Y_{ij}} and that may include time, 
    other time-dependent covariates and time-independent ones. 
    In addition, let \eqn{Y_{i}=(Y_{i1}^{\top},\dots,Y_{in_i}^{\top})^{\top}} denote the \eqn{i}th response vector. 
    With \eqn{\mu_i=E(Y_{i})} and \eqn{\Sigma_i = } cov\eqn{(Y_i)}, the model assumes multivariate normality,
    \eqn{Y_{i} \sim N(\mu_i, \Sigma_i), i=1,2,\dots,n}.
    The means \eqn{\mu_i} and covariance matrices \eqn{\Sigma_i} are modelled semiparametrically in terms of covariates.    
    For the means one can specify semiparametric models, 
    \deqn{\mu_{ijk} = \beta_{k0} + \sum_{l=1}^{K_1} u_{ijl} \beta_{kl} + \sum_{l=K_1+1}^{K} f_{\mu,k,l}(u_{ijl}).}
    This is the first of the 5 regression submodels. 
    
    To model the covariance matrix, first consider the modified Cholesky decomposition, 
    \eqn{L_i \Sigma_i L_i^{\top} = D_i,} where \eqn{L_i} is a unit block lower triangular matrix 
    and \eqn{D_i} is a block diagonal matrix, 
    \deqn{
    	\begin{array}{cc}
    		L_i 
    		= \left[ \begin{array}{cccc}
    			I & 0 & \dots & 0 \\ 
    			-\Phi_{i21} & I & \dots & 0 \\ 
    			\vdots & \vdots & \ddots & \vdots  \\
    			-\Phi_{in_i1} & -\Phi_{in_i1} & \dots & I \\ 
    		\end{array}
    		\right],
    		&
    		D_i 
    		= \left[ \begin{array}{cccc}
    			D_{1} & 0 & \dots & 0 \\ 
    			0 & D_{2} & 0 & 0 \\ 
    			\vdots & \vdots & \ddots & \vdots  \\
    			0 & 0 & 0 & D_{n_i} \\ 
    		\end{array}
    		\right],
    	\end{array}
    }
    
    For modelling \eqn{D_{ij}, i=1,\dots,n, j=1,\dots,n_i} in terms of covariates,   
    first we separate the variances and the correlations \eqn{D_{ij} = S_{ij}^{1/2} R_{ij} S_{ij}^{1/2}}.
    It is easy to model matrix \eqn{S_{ij}} in terms of covariates as the only requirement on its diagonal 
    elements is that they are nonnegative, 
    \deqn{\log \sigma^2_{ijk} = \alpha_{k0} + \sum_{l=1}^{L_1} w_{ijl} \alpha_{kl} + \sum_{l=L_1+1}^{L} f_{\sigma,k,l}(w_{ijl})}
    This is the second of the 5 regression submodels. 
    
    For \eqn{\phi_{ijklm}}, the \eqn{(l,m)} element of \eqn{\Phi_{ijk}, l,m=1,\dots,p}, one can specify semiparametric models
    \deqn{
     \phi_{ijklm} = \psi_{lm0} + \sum_{b=1}^{B_1} v_{ijkb} \psi_{lmb} + \sum_{b=B_1+1}^{B} f_{\phi,l,m,b}(v_{ijkb}) 
    }
    This is the third of the 5 regression submodels. 
    
    The elements of the correlations matrices \eqn{R_{ij}} are modelled in terms of covariate time only, hence they are denoted by \eqn{R_t}. 
    Subject to the positive definiteness constraint, the elements of \eqn{R_t} are modelled using a normal distribution
    with location and scale parameters, \eqn{\mu_{ct}} and \eqn{\sigma^2_{ct}}, modelled as     
    \deqn{
     \mu_{ct} = \eta_0 + f_{\mu}(t), 
    }
    \deqn{ 
      \log \sigma^2_{ct} = \omega_0 + f_{\sigma}(t), 
    }
    and these are the last 2 of the 5 submodels. 
    
}

\value{Function \code{lmrm} returns samples from the posteriors of the model parameters.}

\references{Papageorgiou, G. (2020). Bayesian semiparametric modelling of covariance matrices for multivariate longitudinal data. arXiv:2012.09833.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
	# Fit a joint mean-covariance model on the simulated dataset simD2
	require(ggplot2)
	data(simD2)
	model <- Y1 | Y2 ~ time | sm(time) | sm(lag) | sm(time) | 1 
	# the above defines the responses and the regression models on the left and 
	# right of "~", respectively 
	# the first model, for the mean, is a linear function of time, this is sufficient as 
	# the 2 responses have constant mean.
	# the second model, for the variances, is a smooth function of time
	# the third model, for the dependence structure, is a smooth function of lag, 
	# that lmrm figures out and it does not need to be computed by the user
	# the fourth model, for location of the correlations, is a smooth function of time
	# the fifth model, for scale of the correlations, is just an intercept model
	\dontrun{
	m1 <- lmrm(formula = model, corr.Model = c("common", nClust = 1), data = simD2,
		       id = id, time = time, sweeps = 2500, burn = 500, thin = 2, 
		       StorageDir = getwd(), seed = 1)
	plot(m1)
	}
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
