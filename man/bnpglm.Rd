\name{bnpglm}
\alias{bnpglm}

\title{Bayesian nonparametric generalized linear models}

\description{Fits Dirichlet Process mixtures of joint response-covariate models, where the covariates are
continuous while the discrete response is represented utilizing continuous latent variables. See `Details'
section for a full model description.}

\usage{
bnpglm(formula,family=poisson,data,offset,sampler="slice",WorkingDir,
       ncomp,sweeps,burn,seed,V,Vdf,Mu.nu,Sigma.nu,Mu.mu,Sigma.mu,
       Alpha.gamma,Beta.gamma,Alpha.alpha,Beta.alpha,Turnc.alpha,
       Xpred,offsetPred,...)
}

\arguments{
\item{formula}{a formula defining the response and the covariates e.g. \code{y ~ x}.}

\item{family}{a description of the distribution of the response variable. \code{family=poisson} and
              \code{family=binomial} are the current options.}

\item{data}{an optional data frame, list or environment (or object coercible by `as.data.frame' to a data frame)
            containing the variables in the model.  If not found in `data', the variables are taken from
            `environment(formula)'.}

\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor
              during fitting.  This should be `NULL' or a numeric vector of length equal to the sample size.
              One `offset' term can be included in the formula, and if more are
              required, their sum should be used.}

\item{sampler}{the MCMC algorithm to be utilized. The two options are \code{sampler="slice"} which
               implements the slice sampler (Walker, 2007; Papaspiliopoulos, 2008) and
               \code{sampler="truncated"} which proceeds
               by truncating the countable mixture at \code{ncomp} components (see argument \code{ncomp}).}

\item{WorkingDir}{a directory to store files with the posterior samples of all models parameters. If a
                  directory is not provided, files are not written.}

\item{ncomp}{number of mixture components. Defines where the countable mixture of densities [in (1) below] is truncated.
             Even if  \code{sampler="slice"} is chosen, \code{ncomp} needs to be specified as it is used in the initialization process.}

\item{sweeps}{integer number of posterior samples required.}

\item{burn}{integer number of samples in the burn-in period. In the files described below, \code{sweeps-burn}
            samples will be written.}

\item{seed}{optional seed for the random generator.}

\item{V}{scale matrix \eqn{V} of the prior Wishart distribution assigned to precision matrix \eqn{T_h}.
         See `Details' section.}

\item{Vdf}{degrees of freedom Vdf of the prior Wishart distribution assigned to precision matrix \eqn{T_h}.
           See `Details' section.}

\item{Mu.nu}{prior mean \eqn{\mu_{\nu}} of the covariance vector \eqn{\nu_h}. See `Details' section.}

\item{Sigma.nu}{prior covariance matrix \eqn{\Sigma_{\nu}} of \eqn{\nu_h}. See `Details' section.}

\item{Mu.mu}{prior mean \eqn{\mu_{\mu}} of the mean vector \eqn{\mu_h}. See `Details' section.}

\item{Sigma.mu}{prior covariance matrix \eqn{\Sigma_{\mu}} of \eqn{\mu_h}. See `Details' section.}

\item{Alpha.gamma}{shape parameter \eqn{\alpha_{\gamma}} of the Gamma prior assigned to the Poisson rate \eqn{\gamma_h}.
                   See `Details' section.}

\item{Beta.gamma}{rate parameter \eqn{\beta_{\gamma}} of the Gamma prior assigned to the Poisson rate \eqn{\gamma_h}.
                  See `Details' section.}

\item{Alpha.alpha}{shape parameter \eqn{\alpha_{\alpha}} of the Gamma prior assigned to the concentration parameter \eqn{\alpha}.
                   See `Details' section.}

\item{Beta.alpha}{rate parameter \eqn{\beta_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{Turnc.alpha}{truncation point \eqn{c_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{Xpred}{A design matrix the rows of which include the covariates \eqn{x} for which the conditional distribution
             of \eqn{Y|x,D} (where \eqn{D} denotes the data) is calculated. These are treated as `new' covariates i.e.
             they do not contribute to the likelihood.}

\item{offsetPred}{The offset term associated with the new covariates \code{Xpred}. \code{offsetPred} is a vector of
                  length equal to the rows of \code{Xpred}. If \code{family=poisson}, its entries are the associated
                  Poisson offsets. If \code{family=binomial}, its entries are the Binomial number of trials.}

\item{...}{Other options that will be ignored.}
}

\details{Function \code{bnpglm} returns samples from the posterior distributions of the parameters of the model:
\deqn{
f(y_i,x_i) = \sum_{h=0}^{\infty} \pi_h f(y_i,x_i|\theta_h), \hspace{80pt} (1) }
where \eqn{y_i} is a univariate response with distribution that belongs in the exponential family,
\eqn{x_i} is a \eqn{p}-dimensional vector of continuous covariates, and \eqn{\pi_h, h \geq 1,} are obtained according to Sethuraman's (1994) stick-breaking construction: \eqn{\pi_1 = v_1}, and for \eqn{l \geq 2, \pi_l = v_l \prod_{j=1}^{l-1} (1-v_j)}, where \eqn{v_k} are iid samples \eqn{v_k \sim}Beta\eqn{(1,\alpha), k \geq 1.}

The discrete responses \eqn{y_i} are represented as discretized versions of continuous latent variables \eqn{y_i^*}. Observed discrete and continuous latent variables are connected by:
\deqn{
y_{i} = q \iff c_{h,q-1} < y^*_{i} < c_{h,q},}
where the cut-points are obtained as: \eqn{c_{h,-1} = -\infty},
while for \eqn{q \geq 0}, \eqn{c_{h,q} = c_{q}(\gamma_{h},H_i) = \Phi^{-1}\{F(q;\gamma_{h},H_i)\}.} Here \eqn{\Phi(.)}
is the cumulative distribution function (cdf) of a standard normal variable,
and \eqn{F(.;\gamma,H)} denotes an appropriate cdf. For instance, for modeling count data,
\eqn{F(.;\gamma,H)} denotes the cdf of a Poisson\eqn{(H\gamma)} variable, where \eqn{H} denotes the offset term,
while for modeling Binomial data, \eqn{F(.;\gamma,H)} denotes the cdf of a
Binomial\eqn{(H,\gamma)} variable, where \eqn{H} denotes the number of trials. Further, latent variables are assumed to independently follow a \eqn{y_i^* \sim N(0,1)} distribution, where the mean and variance are restricted to be zero and one as they are non-identifiable by the data.

Joint vectors \eqn{(y_i^{*},x_{i})} are modeled utilizing Gaussian distributions. Then, with \eqn{\theta_h} denoting model parameters associated with the \eqn{h}th component, the joint density \eqn{f(y_{i},x_{i}|\theta_h)} takes the form
\deqn{
f(y_{i},x_{i}|\theta_h) = \int_{c_{i,y_i-1}}^{c_{i,y_i}} N_{p+1}(y_{i}^{*},x_{i}|\mu_{h},C_h) dy_{i}^{*},}
where \eqn{\mu_h} and \eqn{C_h} denote the mean vector and covariance matrix, respectively.

The joint distribution of the latent variable \eqn{y_i^{*}} and the covariates \eqn{x_{i}} is
\deqn{
(y_{i}^{*},x_{i}^T)^T|\theta_h \sim N_{p+1}\left(
\begin{array}{ll}
\left(
\begin{array}{l}
0 \\
\mu_h \\
\end{array}
\right),
 &
C_h=\left[
\begin{array}{ll}
1 &  \nu_h^T \\
\nu_h & \Sigma_h \\
\end{array}
\right]
\end{array}\right),
}
where \eqn{\nu_h} denotes the vector of covariances cov\eqn{(y_{i}^{*},x_{i}|\theta_h)}.
Sampling from the posterior of constrained covariance matrix \eqn{C_h}
is done using methods similar to those of McCulloch et al. (2000).
Specifically, the conditional \eqn{x_{i}|y_{i}^{*} \sim
N_{p}(\mu_h+y_{i}^{*}\nu_h, B_h = \Sigma_h - \nu_h \nu_h^T)} simplifies matters as there
are no constraints on matrix \eqn{B_h} (other than positive definiteness).
Given priors for \eqn{B_h} and \eqn{\nu_h}, it is easy to sample from their posteriors, and thus obtain samples from the posterior of \eqn{\Sigma_h=B_h+\nu_h \nu_h^T}.

\emph{Specification of the prior distributions:}
\enumerate{
\item Define \eqn{T_h=B_h^{-1} = (\Sigma_{h} - \nu_h \nu_h^T)^{-1}, h \geq 1}.
We specify that a priori \eqn{T_h \sim} Wishart\eqn{_{p}(V,}Vdf\eqn{)}, where \eqn{V} is a \eqn{p \times p} scale matrix and Vdf is a scalar degrees of freedom parameter. Default values
are: \eqn{V = I_{p}/p} and Vdf\eqn{=p}, however, these can be changed using arguments \code{V} and
Vdf.

\item The assumed prior is \eqn{\nu_h \sim N_p(\mu_{\nu},\Sigma_{\nu}), h \geq 1}, with default
vaules \eqn{\mu_{\nu}=0} and \eqn{\Sigma_{\nu} = I_{p}}. Arguments \code{Mu.nu} and \code{Sigma.nu} allow the user to change the default values.

\item A priori \eqn{\mu_{h} \sim N_p(\mu_{\mu},\Sigma_{\mu}), h \geq 1}.
Here the default values are \eqn{\mu_{\mu} = \bar{x}} where \eqn{\bar{x}} denotes the sample
mean of the covariates, and \eqn{\Sigma_{\mu} = D} where \eqn{D} denotes a diagonal matrix with diagonal elements equal to the square of the observed range of the covariates. Arguments \code{Mu.mu} and \code{Sigma.mu} allow the user to change the default values.

\item For count data, with \code{family=poisson}, a priori we take
\eqn{\gamma_{h} \sim} Gamma\eqn{(\alpha_{\gamma},\beta_{\gamma}), h \geq 1}.
The default values are \eqn{\alpha_{\gamma}=1.0,\beta_{\gamma}=0.1}, that define a Gamma distribution
with mean \eqn{\alpha_{\gamma}/\beta_{\gamma}=10} and variance \eqn{\alpha_{\gamma}/\beta_{\gamma}^2=100.}
For binomial data, with \code{family=binomial}, a priori we take
\eqn{\gamma_{h} \sim} Beta\eqn{(\alpha_{\gamma},\beta_{\gamma}), h \geq 1}.
The default values are \eqn{\alpha_{\gamma}=1.0,\beta_{\gamma}=1.0}, that define a uniform distribution.
Users can alter the default using using arguments \code{Alpha.gamma} and \code{Beta.gamma}.

\item The concentration parameter \eqn{\alpha} is assigned a Gamma\eqn{(\alpha_{\alpha},\beta_{\alpha})}
prior over the range \eqn{(c_{\alpha},\infty)}, that is,
\eqn{f(\alpha) \propto \alpha^{\alpha_{\alpha}-1} \exp\{-\alpha \beta_{\alpha}\} I[\alpha > c_{\alpha}]},
where \eqn{I[.]} is the indicator function. The default values are \eqn{\alpha_{\alpha}=2.0, \beta_{\alpha}=4.0},
and \eqn{c_{\alpha}=0.25}. Users can alter the default using using arguments \code{Alpha.alpha}, \code{Beta.alpha} and
\code{Turnc.alpha}.
}
}

\value{Function \code{bnpglm} returns the following:
\item{call}{the matched call.}
\item{seed}{the seed that was used (in case replication of the results is needed).}
\item{meanReg}{if \code{Xpred} is specified, the function returns the conditional expectation of the response given each new covariate \eqn{x}.}
\item{medianReg}{if \code{Xpred} is specified, the function returns the conditional median of the response given each new covariate \eqn{x}.}
\item{q1Reg}{if \code{Xpred} is specified, the function returns the conditional first quantile of the response given each new covariate \eqn{x}.}
\item{q3Reg}{if \code{Xpred} is specified, the function returns the conditional third quantile of the response given each new covariate \eqn{x}.}
\item{modeReg}{if \code{Xpred} is specified, the function returns the conditional mode of the response given each new covariate \eqn{x}.}
Further, function \code{bnpglm} creates files where the posterior samples are written. These files are (with all file names
preceded by `BNSP.'):
\item{Th.txt}{this file contains samples from the posteriors of the \eqn{p \times p} precision matrices \eqn{T_h, h=1,2,\dots,ncomp}. The file is arranged in \code{(sweeps-burn)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{T_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of
\eqn{T_h^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{T_h^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{T_h^{(sw)}}.}
\item{Sigmah.txt}{this file contains samples from the posteriors of the \eqn{p \times p} covariance  matrices \eqn{\Sigma_h, h=1,2,\dots,ncomp}. The file is arranged in \code{(sweeps-burn)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\Sigma_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of
\eqn{\Sigma_h^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{\Sigma_h^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{\Sigma_h^{(sw)}}.}
\item{SigmahI.txt}{this file contains samples from the posteriors of the \eqn{p \times p} precision   matrices \eqn{\Sigma_h^{-1}, h=1,2,\dots,ncomp}. The file is arranged in \code{(sweeps-burn)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{(\Sigma_h^{-1})^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{(\Sigma_h^{-1})^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{(\Sigma_h^{-1})^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{(\Sigma_h^{-1})^{(sw)}}.}
\item{nuh.txt}{this file contains samples from the posteriors of the \eqn{p}-dimensional covariance vectors \eqn{\nu_h, h=1,2,\dots,ncomp}. The file is arranged in \code{(sweeps-burn)*ncomp} lines and \eqn{p} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\nu_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{\nu_h^{(sw)}} are written in the columns of the file.}
\item{muh.txt}{this file contains samples from the posteriors of the \eqn{p}-dimensional mean vectors  \eqn{\mu_h, h=1,2,\dots,ncomp}. The file is arranged in \code{(sweeps-burn)*ncomp} lines and \eqn{p} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\mu_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{\mu_h^{(sw)}} are written in the columns of the file.}
\item{gammah.txt}{this file contains samples from the posteriors of the mean parameters  \eqn{\gamma_h}, \eqn{h=1,2,\dots,ncomp}. The file is arranged in
\code{(sweeps-burn)*ncomp} lines and one column. Sweeps write in the file \code{ncomp} lines representing samples \eqn{\gamma_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. }
\item{alpha.txt}{this file contains samples from the posterior of the concentration parameters  \eqn{\alpha}. The file is arranged in \code{(sweeps-burn)} lines and one column, each line including one posterior sample.}
\item{compAlloc.txt}{this file contains the allocations or configurations obtained at each iteration of the sampler. It consists of \code{(sweeps-burn)} lines, that represent the posterior samples, and \eqn{n} columns, that represent the sampling units. Entries in this file range from 0 to \eqn{ncomp-1}.  }
\item{nmembers.txt}{this file contains \code{(sweeps-burn)} lines and \code{ncomp} columns, where the lines represent posterior samples while the columns represent the components or clusters. The entries represent the number of sampling units allocated to the components.  }
\item{Updated.txt}{this file contains \code{(sweeps-burn)} lines with the number of components updated at each iteration of the sampler.}
\item{PD.txt}{this file contains samples from the posterior conditional distribution \eqn{Y|x,D} described in \code{Xpred}.
The file has \code{(sweeps-burn)*npred} lines, where \code{npred} is the number of rows in \code{Xpred}. That is,
at each iteration of the sampler, one line for each `new' covariate vector \eqn{x} is written. The columns of the file represent
the possible values of \eqn{Y}, staring from zero and continuing to a max number.}}



\references{
McCulloch, R. E., Polson, N. G., & Rossi, P. E. (2000). A Bayesian analysis of the multinomial probit
model with fully identified parameters. Journal of Econometrics, 99(1), 173-193.

Papageorgiou, G., Richardson, S. and Best, N. (2014). Bayesian nonparametric models for spatially indexed data of mixed type.

Papaspiliopoulos, O. (2008). A note on posterior sampling from Dirichlet mixture models. Technical report,
University of Warwick.

Sethuraman, J. (1994). A constructive definition of Dirichlet priors. Statistica Sinica, 4, 639-650.

Walker, S. G. (2007). Sampling the Dirichlet mixture model with slices. Communications in Statistics
Simulation and Computation, 36(1), 45-54.
}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
# Bayesian nonparametric GLM with Binomial response Y and one predictor X
data(simD)
pred<-seq(with(simD,min(X))+0.1,with(simD,max(X))-0.1,length.out=30)
npred<-length(pred)
# fit1 and fit2 define the same model but with different numbers of
# components and posterior samples
fit1 <- bnpglm(cbind(Y,(E-Y))~X, family=binomial, data=simD, ncomp=30, sweeps=150,
               burn=100, Xpred=pred, offsetPred=rep(30,npred))
\donttest{fit2 <- bnpglm(cbind(Y,(E-Y))~X, family=binomial, data=simD, ncomp=50, sweeps=5000,
               burn=1000, Xpred=pred, offsetPred=rep(30,npred))
plot(X,Y/E)
lines(pred,fit2$medianReg,col=3,lwd=2)}
}

\keyword{models}
\keyword{regression}
