\name{bnpglm}
\alias{bnpglm}

\title{Bayesian nonparametric generalized linear models}

\description{Fits Dirichlet process mixtures of joint response-covariate models, where the covariates are
continuous while the discrete responses are represented utilizing continuous latent variables. See `Details'
section for a full model description.}

\usage{
bnpglm(formula,family,data,offset,sampler="slice",StorageDir,
       ncomp,sweeps,burn,thin=1,seed,prec,V,Vdf,Mu.nu,Sigma.nu,
       Mu.mu,Sigma.mu,Alpha.xi,Beta.xi,Alpha.alpha,Beta.alpha,Turnc.alpha,
       Xpred,offsetPred,...)
}

\arguments{
\item{formula}{a formula defining the response and the covariates e.g. \code{y ~ x}.}

\item{family}{a description of the distribution of the response variable. Currently five options are
              supported: "poisson", "binomial", "negative binomial", "beta binomial" and "generalized poisson".}

\item{data}{an optional data frame, list or environment (or object coercible by `as.data.frame' to a data frame)
            containing the variables in the model.  If not found in `data', the variables are taken from
            `environment(formula)'.}

\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor
              during fitting.  This should be `NULL' or a numeric vector of length equal to the sample size.
              One `offset' term can be included in the formula, and if more are
              required, their sum should be used.}

\item{sampler}{the MCMC algorithm to be utilized. The two options are \code{sampler="slice"} which
               implements a slice sampler (Walker, 2007; Papaspiliopoulos, 2008) and
               \code{sampler="truncated"} which proceeds
               by truncating the countable mixture at \code{ncomp} components (see argument \code{ncomp}).}

\item{StorageDir}{a directory to store files with the posterior samples of all models parameters. If a
                  directory is not provided, files are created in the current directory and removed when the sampler completes.}

\item{ncomp}{number of mixture components. Defines where the countable mixture of densities [in (1) below] is truncated.
             Even if  \code{sampler="slice"} is chosen, \code{ncomp} needs to be specified as it is used in the initialization process.}

\item{sweeps}{total number of posterior samples, including those discarded in burn-in period (see argument \code{burn}) 
              and those discarded by the thinning process (see argument \code{thin}).}

\item{burn}{length of burn-in period.}

\item{thin}{thinning parameter.}
            
\item{seed}{optional seed for the random generator.}

\item{prec}{precision parameter. Updating of the parameters of the response distribution requires a 
            Metropolis - Hastings step, with proposal distributions centered at current values and 
            with precision equal to this argument. It can be of length one ("poisson" and "binomial"
            families) or of length two ("negative binomial", "beta binomial" and "generalized-poisson").}

\item{V}{optional scale matrix \eqn{V} of the prior Wishart distribution assigned to precision matrix \eqn{T_h}.
         See `Details' section.}

\item{Vdf}{optional degrees of freedom Vdf of the prior Wishart distribution assigned to precision matrix \eqn{T_h}.
           See `Details' section.}

\item{Mu.nu}{optional prior mean \eqn{\mu_{\nu}} of the covariance vector \eqn{\nu_h}. See `Details' section.}

\item{Sigma.nu}{optional prior covariance matrix \eqn{\Sigma_{\nu}} of \eqn{\nu_h}. See `Details' section.}

\item{Mu.mu}{optional prior mean \eqn{\mu_{\mu}} of the mean vector \eqn{\mu_h}. See `Details' section.}

\item{Sigma.mu}{optional prior covariance matrix \eqn{\Sigma_{\mu}} of \eqn{\mu_h}. See `Details' section.}

\item{Alpha.xi}{an optional parameter that depends on the specified family.
                \enumerate{
                \item If \code{family="poisson"}, this argument is parameter \eqn{\alpha_{\xi}} of the prior of the Poisson rate:
                \eqn{\xi \sim} Gamma(\eqn{\alpha_{\xi},\beta_{\xi}}).
                
                \item If \code{family="binomial"}, this argument is parameter \eqn{\alpha_{\xi}} of the prior of the Binomial probability:
                \eqn{\xi \sim} Beta(\eqn{\alpha_{\xi},\beta_{\xi}}).
                
                \item If \code{family="netative binomial"},
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Negative Binomial pmf.
                
                \item If \code{family="beta binomial"}, 
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Beta Binomial pmf.
                
                \item If \code{family="generalized-poisson"}, this argument is a two-dimensional vector 
                that includes parameters \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Normal(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Generalized Poisson pmf.
                }
                See `Details' section.}

\item{Beta.xi}{an optional parameter that depends on the specified family.
               \enumerate{
               \item If \code{family="poisson"}, this argument is parameter \eqn{\beta_{\xi}} of the prior of the Poisson rate:
               \eqn{\xi \sim} Gamma(\eqn{\alpha_{\xi},\beta_{\xi}}).
               
               \item If \code{family="binomial"}, this argument is parameter \eqn{\beta_{\xi}} of the prior of the Binimial probability:
               \eqn{\xi \sim} Beta(\eqn{\alpha_{\xi},\beta_{\xi}}).
               
               \item If \code{family="netative binomial"},
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Negative Binomial pmf.
                
                \item If \code{family="beta binomial"}, 
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Beta Binomial pmf.
                
                \item If \code{family="generalized-poisson"}, this argument is a two-dimensional vector 
                that includes parameters \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Normal(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Generalized Poisson pmf.
               }
               See `Details' section.}

\item{Alpha.alpha}{optional shape parameter \eqn{\alpha_{\alpha}} of the Gamma prior assigned to the concentration parameter \eqn{\alpha}.
                   See `Details' section.}

\item{Beta.alpha}{optional rate parameter \eqn{\beta_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{Turnc.alpha}{optional truncation point \eqn{c_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{Xpred}{an optional design matrix the rows of which include the covariates \eqn{x} for which the conditional distribution
             of \eqn{Y|x,D} (where \eqn{D} denotes the data) is calculated. These are treated as `new' covariates i.e.
             they do not contribute to the likelihood. The matrix shouldn't include a column of 1's.}

\item{offsetPred}{the offset term associated with the new covariates \code{Xpred}. \code{offsetPred} is a vector of
                  length equal to the rows of \code{Xpred}. If \code{family} is one of \code{poisson} or 
                  \code{negative binomial} or \code{generalized poisson}, its entries are the associated
                  Poisson offsets. If \code{family} is one of \code{binomial} or \code{beta binomial}, 
                  its entries are the Binomial number of trials. If \code{offsetPred} is missing, it is taken to be the
                  mean of \code{offset}, rounded to the nearest integer.}

\item{...}{Other options that will be ignored.}
}

\details{Function \code{bnpglm} returns samples from the posterior distributions of the parameters of the model:
\deqn{
f(y_i,x_i) = \sum_{h=0}^{\infty} \pi_h f(y_i,x_i|\theta_h), \hspace{80pt} (1) }
where \eqn{y_i} is a univariate discrete response,
\eqn{x_i} is a \eqn{p}-dimensional vector of continuous covariates, and \eqn{\pi_h, h \geq 1,} are obtained according to
Sethuraman's (1994) stick-breaking construction:
\eqn{\pi_1 = v_1}, and for \eqn{l \geq 2, \pi_l = v_l \prod_{j=1}^{l-1} (1-v_j)}, where \eqn{v_k} are iid samples
\eqn{v_k \sim}Beta \eqn{(1,\alpha), k \geq 1.}

The discrete responses \eqn{y_i} are represented as discretized versions of continuous latent variables \eqn{y_i^*}.
Observed discrete and continuous latent variables are connected by:
\deqn{
y_{i} = q \iff c_{i,q-1} < y^*_{i} < c_{i,q},}
where the cut-points are obtained as: \eqn{c_{i,-1} = -\infty},
while for \eqn{q \geq 0}, \eqn{c_{i,q} = c_{q}(\lambda_{i}) = \Phi^{-1}\{F(q;\lambda_i)\}.} Here \eqn{\Phi(.)}
is the cumulative distribution function (cdf) of a standard normal variable
and \eqn{F()} denotes an appropriate cdf. Further, latent variables are assumed to
independently follow a \eqn{N(0,1)} distribution, where the mean and variance are restricted to be zero and one as
they are non-identifiable by the data. Choices for \eqn{F()} are described next.

For counts, currently three options are supported. First, \eqn{F(.;\lambda_i)} can be specified as the 
cdf of a Poisson\eqn{(H_i \xi_h)} variable. Here \eqn{\lambda_i=(\xi_h,H_i)^T, \xi_h} denotes the Poisson rate 
associated with cluster \eqn{h}, and \eqn{H_i} the offset term associated with sampling unit \eqn{i}.
Second, \eqn{F(.;\lambda_i)} can be specified as the negative binomial cdf, where \eqn{\lambda=(\xi_{1h},\xi_{2h},H_i)^T}.
This option allows for overdispersion within each cluster. With this added flexibility the model usually requires 
fewer effective clusters. Third, \eqn{F(.;\lambda_i)} can be specified as the Generalized Poisson cdf, where, again, 
\eqn{\lambda=(\xi_{1h},\xi_{2h},H_i)^T}. This option allows for both over- and under-dispersion within each
cluster. See Consul and Famoye (1992) for more details. 

For Binomial data, currently two options are supported. First, \eqn{F(.;\lambda_i)} may be taken to be the cdf of a
Binomial\eqn{(H_i,\xi_h)} variable, with \eqn{\xi_h} denoting the success probability of cluster \eqn{h} and \eqn{H_i} 
the number of trials associated with sampling unit \eqn{i}.
Second, \eqn{F(.;\lambda_i)} may be specified to be the beta-binomial cdf, where \eqn{\lambda=(\xi_{1h},\xi_{2h},H_i)^T}.

Joint vectors \eqn{(y_i^{*},x_{i})} are modeled utilizing Gaussian distributions. Then, with \eqn{\theta_h} denoting model parameters
associated with the \eqn{h}th cluster, the joint density \eqn{f(y_{i},x_{i}|\theta_h)} takes the form
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

\item The assumed prior for \eqn{\nu_h} is \eqn{N_p(\mu_{\nu},\Sigma_{\nu}), h \geq 1}, with default
vaules \eqn{\mu_{\nu}=0} and \eqn{\Sigma_{\nu} = I_{p}}. Arguments \code{Mu.nu} and \code{Sigma.nu} allow the user to change the default values.

\item A priori \eqn{\mu_{h} \sim N_p(\mu_{\mu},\Sigma_{\mu}), h \geq 1}.
Here the default values are \eqn{\mu_{\mu} = \bar{x}} where \eqn{\bar{x}} denotes the sample
mean of the covariates, and \eqn{\Sigma_{\mu} = D} where \eqn{D} denotes a diagonal matrix with diagonal elements equal to the square of the observed range of the covariates. Arguments \code{Mu.mu} and \code{Sigma.mu} allow the user to change the default values.

\item For count data, with \code{family="poisson"}, a priori we take
\eqn{\xi_{h} \sim} Gamma\eqn{(\alpha_{\xi},\beta_{\xi}), h \geq 1}.
The default values are \eqn{\alpha_{\xi}=1.0,\beta_{\xi}=0.1}, that define a Gamma distribution
with mean \eqn{\alpha_{\xi}/\beta_{\xi}=10} and variance \eqn{\alpha_{\xi}/\beta_{\xi}^2=100.}
Defaults can be altered using arguments \code{Alpha.xi} and \code{Beta.xi}.
For binomial data, with \code{family="binomial"}, a priori we take
\eqn{\xi_{h} \sim} Beta\eqn{(\alpha_{\xi},\beta_{\xi})}, \eqn{h \geq 1}.
The default values are \eqn{\alpha_{\xi}=1.0,\beta_{\xi}=1.0}, that define a uniform distribution.
Defaults can be altered using arguments \code{Alpha.xi} and \code{Beta.xi}.
For both count data with \code{family="negative binomial"} and binomial data with \code{family="beta binomial"}, 
a priori we take \eqn{\xi_{jh} \sim} Gamma\eqn{(\alpha_{j\xi},\beta_{j\xi})}, \eqn{j=1,2, h \geq 1}.
The default values are \eqn{\alpha_{j\xi}=1.0,\beta_{j\xi}=0.1}.
Default values for \eqn{\{\alpha_{j\xi}: j=1,2\}} can be altered using argument \code{Alpha.xi},
and default values for \eqn{\{\beta_{j\xi}: j=1,2\}} can be altered using argument \code{Beta.xi}.
For count data with \code{family="generalized poisson"}, 
a priori we take \eqn{\xi_{1h} \sim} Gamma\eqn{(\alpha_{1\xi},\beta_{1\xi})}, 
and \eqn{\xi_{2h} \sim} Normal\eqn{(\alpha_{2\xi},\beta_{2\xi})}.
The default values are \eqn{\alpha_{j\xi}=1.0,\beta_{j\xi}=0.1}.
Default values for \eqn{\{\alpha_{j\xi}: j=1,2\}} can be altered using argument \code{Alpha.xi},
and default values for \eqn{\{\beta_{j\xi}: j=1,2\}} can be altered using argument \code{Beta.xi}.


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
\item{Th.txt}{this file contains samples from the posteriors of the \eqn{p \times p} precision matrices \eqn{T_h, h=1,2,\dots,ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{T_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of
\eqn{T_h^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{T_h^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{T_h^{(sw)}}.}
\item{Sigmah.txt}{this file contains samples from the posteriors of the \eqn{p \times p} covariance  matrices \eqn{\Sigma_h, h=1,2,\dots,ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\Sigma_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of
\eqn{\Sigma_h^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{\Sigma_h^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{\Sigma_h^{(sw)}}.}
\item{SigmahI.txt}{this file contains samples from the posteriors of the \eqn{p \times p} precision   matrices \eqn{\Sigma_h^{-1}, h=1,2,\dots,ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p^2} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{(\Sigma_h^{-1})^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{(\Sigma_h^{-1})^{(sw)}} are written in the columns of the file: the entries in the first \eqn{p} columns of the file are those in the first column (or row) of \eqn{(\Sigma_h^{-1})^{(sw)}}, while the entries in the last \eqn{p} columns of the file are those in the last  column (or row) of \eqn{(\Sigma_h^{-1})^{(sw)}}.}
\item{nuh.txt}{this file contains samples from the posteriors of the \eqn{p}-dimensional covariance vectors \eqn{\nu_h, h=1,2,\dots,ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\nu_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{\nu_h^{(sw)}} are written in the columns of the file.}
\item{muh.txt}{this file contains samples from the posteriors of the \eqn{p}-dimensional mean vectors  \eqn{\mu_h, h=1,2,\dots,ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p} columns. In more detail, each sweep creates \code{ncomp} lines representing samples \eqn{\mu_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. The elements of \eqn{\mu_h^{(sw)}} are written in the columns of the file.}
\item{xih.txt}{this file contains samples from the posteriors of parameters  \eqn{\xi_h}, \eqn{h=1,2,\dots,ncomp}.
The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and one or two columns, depending on the number of parameters in the selected \eqn{F(.;\lambda)}. Sweeps write in the file \code{ncomp} lines representing samples \eqn{\xi_h^{(sw)}, h=1,\dots,ncomp}, where superscript \eqn{sw} represents a particular sweep. }
\item{alpha.txt}{this file contains samples from the posterior of the concentration parameters  \eqn{\alpha}. The file is arranged in \code{(sweeps-burn)/thin} lines and one column, each line including one posterior sample.}
\item{compAlloc.txt}{this file contains the allocations or configurations obtained at each iteration of the sampler. It consists of \code{(sweeps-burn)/thin} lines, that represent the posterior samples, and \eqn{n} columns, that represent the sampling units. Entries in this file range from 0 to \eqn{ncomp-1}.  }
\item{nmembers.txt}{this file contains \code{(sweeps-burn)/thin} lines and \code{ncomp} columns, where the lines represent posterior samples while the columns represent the components or clusters. The entries represent the number of sampling units allocated to the components.  }
\item{Updated.txt}{this file contains \code{(sweeps-burn)/thin} lines with the number of components updated at each iteration of the sampler.}
%\item{PD.txt}{this file contains samples from the posterior conditional distribution \eqn{Y|x,D} described in \code{Xpred}.
%The file has \code{(sweeps-burn)*npred} lines, where \code{npred} is the number of rows in \code{Xpred}. That is,
%at each iteration of the sampler, one line for each `new' covariate vector \eqn{x} is written. The columns of the file represent
%the possible values of \eqn{Y}, staring from zero and continuing to a max number.}
}



\references{
Consul, P. C. & Famoye, G. C. (1992). Generalized Poisson regression model. Communications in 
Statistics - Theory and Methods, 1992, 89-109.

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
# components and posterior samples. They both use a slice sampler
# and parameter prec=200 achieves optimal acceptance rate, about 22%.
fit1 <- bnpglm(cbind(Y,(E-Y))~X, family="binomial", data=simD, ncomp=30, sweeps=150,
               burn=100, sampler="slice", prec1=c(200), Xpred=pred, offsetPred=rep(30,npred))
\donttest{fit2 <- bnpglm(cbind(Y,(E-Y))~X, family="binomial", data=simD, ncomp=50, sweeps=5000,
               burn=1000, sampler="slice", prec=c(200), Xpred=pred, offsetPred=rep(30,npred))
plot(with(simD,X),with(simD,Y)/with(simD,E))
lines(pred,fit2$medianReg,col=3,lwd=2)}
}

\keyword{nonparametric}
\keyword{cluster}




