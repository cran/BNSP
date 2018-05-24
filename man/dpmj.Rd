\name{dpmj}
\alias{dpmj}

\title{Dirichlet process mixtures of joint models}

\description{Fits Dirichlet process mixtures of joint response-covariate models, where the covariates are
of mixed type while the discrete responses are represented utilizing continuous latent variables. See `Details'
section for a full model description and Papageorgiou (2018) for all technical details.}

\usage{
dpmj(formula, Fcdf, data, offset, sampler = "truncated", Xpred, offsetPred,
     StorageDir, ncomp, sweeps, burn, thin = 1, seed, H, Hdf, d, D,
     Alpha.xi, Beta.xi, Alpha.alpha, Beta.alpha, Trunc.alpha, ...)
}

\arguments{
\item{formula}{a formula defining the response and the covariates e.g. \code{y ~ x}.}

\item{Fcdf}{a description of the kernel of the response variable. Currently five options are
              supported: 1. "poisson", 2. "negative binomial", 3. "generalized poisson", 4. "binomial" and 
              5. "beta binomial". The first three kernels are used for count data analysis, where the third
              kernel allows for both over- and under-dispersion relative to the Poisson distribution.   
              The last two kernels are used for binomial data analysis. See `Details' section for some of the 
              kernel details.}

\item{data}{an optional data frame, list or environment (or object coercible by `as.data.frame' to a data frame)
            containing the variables in the model.  If not found in `data', the variables are taken from
            `environment(formula)'.}

\item{offset}{this can be used to specify an a priori known component to be included in the model. 
              This should be `NULL' or a numeric vector of length equal to the sample size.
              One `offset' term can be included in the formula, and if more are
              required, their sum should be used.}

\item{sampler}{the MCMC algorithm to be utilized. The two options are \code{sampler = "slice"} which
               implements a slice sampler (Walker, 2007; Papaspiliopoulos, 2008) and
               \code{sampler = "truncated"} which proceeds
               by truncating the countable mixture at \code{ncomp} components (see argument \code{ncomp}).}

\item{Xpred}{an optional design matrix the rows of which include the values of the covariates \eqn{x} for which the conditional distribution
             of \eqn{Y|x,D} (where \eqn{D} denotes the data) is calculated. These are treated as `new' covariates i.e.
             they do not contribute to the likelihood. The matrix shouldn't include a column of 1's. NA's can be 
             included to obtain averaged effects.}

\item{offsetPred}{the offset term associated with the new covariates \code{Xpred}. It is of dimension one i.e.
                  the same offset term is used for all rows of \code{Xpred}. If \code{Fcdf} is one of "poisson" or 
                  "negative binomial" or "generalized poisson", then \code{offsetPred} is the Poisson offset term. 
                  If \code{Fcdf} is one of "binomial" or "beta binomial", then \code{offsetPred} is the number of Binomial trials. 
                  If \code{offsetPred} is missing, it is taken to be the mean of \code{offset}, rounded to the nearest integer.}
               
\item{StorageDir}{a directory to store files with the posterior samples of models parameters and other quantities 
                  of interest. If a directory is not provided, files are created in the current directory and removed 
                  when the sampler completes.}

\item{ncomp}{number of mixture components. It defines where the countable mixture of densities [in (1) below] is truncated.
             Even if  \code{sampler="slice"} is chosen, \code{ncomp} needs to be specified as it is used in the initialization process.}

\item{sweeps}{total number of posterior samples, including those discarded in burn-in period (see argument \code{burn}) 
              and those discarded by the thinning process (see argument \code{thin}).}

\item{burn}{length of burn-in period.}

\item{thin}{thinning parameter.}
            
\item{seed}{optional seed for the random generator.}

\item{H}{optional scale matrix of the Wishart-like prior assigned to the restricted covariance matrices \eqn{\Sigma_h^*}.
         See `Details' section.}

\item{Hdf}{optional degrees of freedom of the prior Wishart-like prior assigned to the restricted covariance matrices \eqn{\Sigma_h^*}.
           See `Details' section.}

\item{d}{optional prior mean of the mean vector \eqn{\mu_h}. See `Details' section.}

\item{D}{optional prior covariance matrix of the mean vector \eqn{\mu_h}. See `Details' section.}

\item{Alpha.xi}{an optional parameter that depends on the specified \code{Fcdf} argument.
                \enumerate{
                \item If \code{Fcdf = "poisson"}, this argument is parameter \eqn{\alpha_{\xi}} of the prior of the Poisson rate:
                \eqn{\xi \sim} Gamma(\eqn{\alpha_{\xi},\beta_{\xi}}).
                
                \item If \code{Fcdf = "negative binomial"},
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Negative Binomial pmf.
                
                \item If \code{Fcdf = "generalized poisson"}, this argument is a two-dimensional vector 
                that includes parameters \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } N(\eqn{\alpha_{2\xi},\beta_{2\xi})I[\xi_2 \in R_{\xi_2}]}, where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Generalized Poisson pmf. Parameter \eqn{\xi_2} is restricted in the range
                \eqn{R_{\xi_2} = (0.05,\infty)} as it is a dispersion parameter. 
                
                \item If \code{Fcdf = "binomial"}, this argument is parameter \eqn{\alpha_{\xi}} of the prior of the Binomial probability:
                \eqn{\xi \sim} Beta(\eqn{\alpha_{\xi},\beta_{\xi}}).
                
                \item If \code{Fcdf = "beta binomial"}, 
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\alpha_{1\xi}} and \eqn{\alpha_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Beta Binomial pmf.
                }
                See `Details' section.}

\item{Beta.xi}{an optional parameter that depends on the specified family.
                \enumerate{
                \item If \code{Fcdf = "poisson"}, this argument is parameter \eqn{\beta_{\xi}} of the prior of the Poisson rate:
                \eqn{\xi \sim} Gamma(\eqn{\alpha_{\xi},\beta_{\xi}}).
               
                \item If \code{Fcdf = "negative binomial"},
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Negative Binomial pmf.
                
                \item If \code{Fcdf = "generalized poisson"}, this argument is a two-dimensional vector 
                that includes parameters \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Normal(\eqn{\alpha_{2\xi},\beta_{2\xi})I[\xi_2 \in R_{\xi_2}]}, where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Generalized Poisson pmf. 
                Parameter \eqn{\xi_2} is restricted in the range
                \eqn{R_{\xi_2} = (0.05,\infty)} as it is a dispersion parameter.
                Note that \eqn{\beta_{2\xi}} is a standard deviation.                          
                
                \item If \code{Fcdf = "binomial"}, this argument is parameter \eqn{\beta_{\xi}} of the prior of the Binomial probability:
                \eqn{\xi \sim} Beta(\eqn{\alpha_{\xi},\beta_{\xi}}).
               
                \item If \code{Fcdf = "beta binomial"}, 
                this argument is a two-dimensional vector that includes parameters 
                \eqn{\beta_{1\xi}} and \eqn{\beta_{2\xi}} of the priors:
                \eqn{\xi_1 \sim } Gamma(\eqn{\alpha_{1\xi},\beta_{1\xi}}) and
                \eqn{\xi_2 \sim } Gamma(\eqn{\alpha_{2\xi},\beta_{2\xi}}), where \eqn{\xi_1} and \eqn{\xi_2}
                are the two parameters of the Beta Binomial pmf.
               }
               See `Details' section.}

\item{Alpha.alpha}{optional shape parameter \eqn{\alpha_{\alpha}} of the Gamma prior assigned to the concentration parameter \eqn{\alpha}.
                   See `Details' section.}

\item{Beta.alpha}{optional rate parameter \eqn{\beta_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{Trunc.alpha}{optional truncation point \eqn{c_{\alpha}} of the Gamma prior assigned to concentration parameter \eqn{\alpha}.
                  See `Details' section.}

\item{...}{Other options that will be ignored.}
}

\details{
Function \code{dpmj} returns samples from the posterior distributions of the parameters of the model:
\deqn{
f(y_i,x_i) = \sum_{h=1}^{\infty} \pi_h f(y_i,x_i|\theta_h), \hspace{200pt} (1) }
where \eqn{y_i} is a univariate discrete response,
\eqn{x_i} is a \eqn{p}-dimensional vector of mixed type covariates, and \eqn{\pi_h, h \geq 1,} are obtained according to
Sethuraman's (1994) stick-breaking construction:
\eqn{\pi_1 = v_1}, and for \eqn{l \geq 2, \pi_l = v_l \prod_{j=1}^{l-1} (1-v_j)}, where \eqn{v_k} are iid samples
\eqn{v_k \sim}Beta \eqn{(1,\alpha), k \geq 1.}

Let \eqn{Z} denote a discrete variable (response or covariate). It is represented as discretized version of a continuous 
latent variable \eqn{Z^*}.
Observed discrete \eqn{Z} and continuous latent variable \eqn{Z^*} are connected by:
\deqn{
z = q \iff c_{q-1} < z^* < c_{q}, q=0,1,2,\dots,}
where the cut-points are obtained as: \eqn{c_{-1} = -\infty},
while for \eqn{q \geq 0}, \eqn{c_{q} = c_{q}(\lambda) = \Phi^{-1}\{F(q;\lambda)\}.} Here \eqn{\Phi(.)}
is the cumulative distribution function (cdf) of a standard normal variable
and \eqn{F()} denotes an appropriate cdf. Further, latent variables are assumed to
independently follow a \eqn{N(0,1)} distribution, where the mean and variance are restricted to be zero and one as
they are non-identifiable by the data. Choices for \eqn{F()} are described next.

For counts, three options are supported. First, \eqn{F(.;\lambda_i)} can be specified as the 
cdf of a Poisson\eqn{(H_i \xi_h)} variable. Here \eqn{\lambda_i=(\xi_h,H_i)^T, \xi_h} denotes the Poisson rate 
associated with cluster \eqn{h}, and \eqn{H_i} the offset term associated with sampling unit \eqn{i}.
Second, \eqn{F(.;\lambda_i)} can be specified as the negative binomial cdf, where \eqn{\lambda_i=
(\xi_{1h},\xi_{2h},H_i)^T}. This option allows for overdispersion within each cluster relative to the 
Poisson distribution. Third, \eqn{F(.;\lambda_i)} can be specified as the Generalized Poisson cdf, where, again, 
\eqn{\lambda_i=(\xi_{1h},\xi_{2h},H_i)^T}. This option allows for both over- and under-dispersion within each
cluster. 

For Binomial data, two options are supported. First, \eqn{F(.;\lambda_i)} may be taken to be the cdf of a
Binomial\eqn{(H_i,\xi_h)} variable, where \eqn{\xi_h} denotes the success probability of cluster \eqn{h} and \eqn{H_i} 
the number of trials associated with sampling unit \eqn{i}.
Second, \eqn{F(.;\lambda_i)} may be specified to be the beta-binomial cdf, where \eqn{\lambda=(\xi_{1h},\xi_{2h},H_i)^T}.

The special case of Binomial data is treated as 
\deqn{
Z = 0 \iff z^* < 0, z^* \sim N(\mu_z^{*},1).} 

Details on all kernels are provided in the two tables below. The first table provides the probability mass functions
and the mean in the presence of an offset term (which may be taken to be one). The column `Sample' indicates
for which parameters the routine provides posterior samples. The second table provides information on the assumed priors
along with the default values of the parameters of the prior distributions and it also indicates the 
function arguments that allow the user to alter these.  

\tabular{llcll}{
Kernel  \tab PMF \tab Offset \tab Mean \tab Sample \cr
Poisson \tab \eqn{\exp(-H\xi) (H\xi)^y /y!} \tab \eqn{H}      \tab \eqn{H \xi}  \tab \eqn{\xi} \cr 
Negative Binomial \tab \eqn{\frac{\Gamma(y+\xi_1)}{\Gamma(\xi_1)\Gamma(y+1)}(\frac{\xi_2}{H+\xi_2})^{\xi_1}(\frac{H}{H+\xi_2})^{y}} 
\tab \eqn{H} \tab \eqn{H \xi_1/\xi_2} \tab \eqn{\xi_1, \xi_2}\cr 
Generalized Poisson \tab \eqn{\xi_1 \{\xi_1+(\xi_2-1)y\}^{y-1} \xi_2^{-y} \times} \tab \eqn{H} \tab \eqn{H\xi_1} \tab  \eqn{\xi_1,\xi_2}\cr
\tab \eqn{ ~~ \exp\{-[\xi_1+(\xi_2-1)y]/\xi_2\}/y!} \tab \tab \tab \cr   
%Hyper Poisson \tab \eqn{\frac{1}{_1F_1(1,\xi_2,\xi_3)} \frac{\xi_3^y}{(\xi_2)_y}} \tab \eqn{H} \tab
%\eqn{H \xi_1 = \xi_3 - }  \tab \eqn{\xi_1,\xi_2} \cr                                                                 
%\tab \tab \tab \eqn{ ~~ (\xi_2-1) \frac{_1F_1(1,\xi_2,\xi_3)-1}{_1F_1(1,\xi_2,\xi_3)}}  \tab \cr                                                                 
%CTP \tab \eqn{f_0 \frac{(\xi_3+\xi_4 i)_y (\xi_3-\xi_4 i)_y}{(\xi_2)_y y!}} \tab \eqn{H} \tab 
%\eqn{H \xi_1 = \frac{\xi_3^2+\xi_4^2}{\xi_2-2\xi_3-1}} \tab \eqn{\xi_1, \xi_2, \xi_3}\cr
%COM-Poisson \tab \eqn{\frac{\xi_3^y}{Z(\xi_2,\xi_3)(y!)^{\xi_2}}} \tab 
%\eqn{H} \tab \eqn{H \xi_1 = \xi_3 \frac{\partial \log(Z)}{\partial \xi_3}} \tab \eqn{\xi_1,\xi_2}\cr                                                                 
Binomial \tab \eqn{{N \choose y} \xi^y (1-\xi)^{N-y}} \tab \eqn{N} \tab \eqn{N \xi} \tab \eqn{\xi}\cr 
Beta Binomial \tab \eqn{{N \choose y} \frac{{Beta}{(y+\xi_1,N-y+\xi_2)}}{{Beta}{(\xi_1,\xi_2)}}} 
\tab \eqn{N} \tab \eqn{N \xi_1/(\xi_1+\xi_2)} \tab \eqn{\xi_1,\xi_2}\cr
}

\tabular{lll}{
Kernel  \tab Priors \tab Default Values \cr
Poisson \tab \eqn{\xi \sim} Gamma\eqn{(\alpha_{\xi},\beta_{\xi})} \tab Alpha.xi = 1.0, Beta.xi = 0.1 \cr
Negative Binomial \tab \eqn{\xi_i \sim} Gamma\eqn{(\alpha_{\xi_i},\beta_{\xi_i}), i=1,2} \tab Alpha.xi = c(1.0,1.0), Beta.xi = c(0.1,0.1) \cr
Generalized Poisson \tab \eqn{\xi_1 \sim} Gamma\eqn{(\alpha_{\xi_1},\beta_{\xi_1})} \tab \cr
                    \tab \eqn{\xi_2 \sim} N\eqn{(\alpha_{\xi_2},\beta_{\xi_2})I[\xi_2 > 0.05]} \tab Alpha.xi = c(1.0,1.0), Beta.xi = c(0.1,1.0) \cr
                    \tab where \eqn{\beta_{\xi_2}} denotes st.dev. \tab \cr
%Hyper Poisson \tab \eqn{\xi_i \sim} Gamma\eqn{(\alpha_{\xi_i},\beta_{\xi_i}), i=1,2} \tab Alpha.xi = c(1.0,0.5), Beta.xi = c(0.1,0.5) \cr
%CTP \tab \eqn{\xi_i \sim} Gamma\eqn{(\alpha_{\xi_i},\beta_{\xi_i}), i=1,2} \tab \cr
%    \tab \eqn{\xi_3 \sim TN(\alpha_{\xi_3},\beta_{\xi_3})} (\eqn{\beta_{\xi_3} \equiv} st.dev.) \tab Alpha.xi = c(1.0,1.0,0.0) \cr
%    \tab TN: truncated normal \tab Beta.xi = c(0.1,0.1,100.0) \cr  
%COM-Poisson \tab \eqn{\xi_i \sim} Gamma\eqn{(\alpha_{\xi_i},\beta_{\xi_i}), i=1,2} \tab Alpha.xi = c(1.0,0.5), Beta.xi = c(0.1,0.5) \cr
Binomial \tab \eqn{\xi \sim} Beta\eqn{(\alpha_{\xi},\beta_{\xi})} \tab Alpha.xi = 1.0, Beta.xi = 1.0 \cr
Beta Binomial \tab \eqn{\xi_i \sim} Gamma\eqn{(\alpha_{\xi_i},\beta_{\xi_i}), i=1,2} \tab Alpha.xi = c(1.0,1.0), Beta.xi = c(0.1,0.1) \cr
}

%\tabular{ll}{
%Kernel  \tab Notes\cr
%Generalized Poisson \tab \eqn{\xi_1 > 0} is the mean and \eqn{\xi_2 > 1/2} is a dispersion parameter. When \eqn{\xi_2 = 1}, \cr
%                    \tab  the pmf reduces to the Poisson. Parameter values \eqn{\xi_2 > 1} suggest over- \cr 
%                    \tab dispersion and parameter values \eqn{1/2 < \xi_2  < 1} suggest under-dispersion  \cr
%                    \tab relative to the Poisson. \cr
%Hyper Poisson \tab \eqn{\xi_1 > 0} is the mean and \eqn{\xi_2 > 0} is a dispersion parameter. When \eqn{\xi_2 = 1}, \cr
%              \tab  the pmf reduces to the Poisson. When \eqn{\xi_2 > 1} the pmf is over-dispersed \cr 
%              \tab and when \eqn{\xi_2  < 1} the pmf is under-dispersed relative to the Poisson. \cr
%COM-Poisson \tab The mean is \eqn{\xi_1 (> 0)} and the variance approximately \eqn{\xi_1/\xi_2}, so\cr
%            \tab similar comments as for the hyper Poisson hold.\cr
%CTPD \tab Things are a bit more complex here. See Rodriguez-Avi et al. (2004) for the details. \cr
%}

Let \eqn{z_i = (y_i,x_{i}^T)^T} denote the joint vector of observed continuous and discrete variables and \eqn{z_i^*}
the corresponding vector of continuous observed and latent variables. With \eqn{\theta_h} denoting model parameters
associated with the \eqn{h}th cluster, the joint density \eqn{f(z_{i}|\theta_h)} takes the form
\deqn{
f(z_i|\theta_h) = \int_{R(y)} \int_{R(x_{d})} N_{q}(z_i^*;\mu^*_h,\Sigma^*_h) dx_{d}^{*} dy^{*},}
where \deqn{
\begin{array}{ll}
\mu^*_h = \left(
\begin{array}{l}
0 \\
\mu_h \\
\end{array}
\right),
 &
\Sigma^*_h=\left[
\begin{array}{ll}
C_h &  \nu_h^T \\
\nu_h & \Sigma_h \\
\end{array}
\right]
\end{array},
}
where \eqn{C_h} is the covariance matrix of the latent continuous variables and it has  
diagonal elements equal to one i.e. it is a correlation matrix.

In addition to the priors defined in the table above, we specify the following: 
\enumerate{
\item 
The restricted covariance matrix \eqn{\Sigma^*_h} is assigned a prior distribution that is based on the Wishart
distribution with degrees of freedom set by default to dimension of matrix plus two and diagonal scale matrix,
with the sub-matrix that corresponds to discrete variables taken to be the identity matrix and with sub-matrix 
that corresponds to continuous variables having entries equal to 1/8 of the square of 
the observed data range. Default values can be changed using arguments \code{H} and \code{Hdf}.

\item
The prior on \eqn{\mu_h}, the non-zero part of \eqn{\mu_h^*}, is taken to be multivariate normal \eqn{\mu_h \sim N(d,D)}. 
The mean \eqn{d} is taken to be equal to the center of the dataset. The covariance matrix \eqn{D} is taken to be diagonal. 
Its elements that correspond to continuous variables are set equal to 1/8 of the square of the observed data range while the 
elements that correspond to binary variables are set equal to 5.  
Arguments \code{Mu.mu} and \code{Sigma.mu} allow the user to change the default values.

\item 
The concentration parameter \eqn{\alpha} is assigned a Gamma\eqn{(\alpha_{\alpha},\beta_{\alpha})}
prior over the range \eqn{(c_{\alpha},\infty)}, that is,
\eqn{f(\alpha) \propto \alpha^{\alpha_{\alpha}-1} \exp\{-\alpha \beta_{\alpha}\} I[\alpha > c_{\alpha}]},
where \eqn{I[.]} is the indicator function. The default values are \eqn{\alpha_{\alpha}=2.0, \beta_{\alpha}=5.0},
and \eqn{c_{\alpha}=0.25}. Users can alter the default using using arguments \code{Alpha.alpha}, \code{Beta.alpha} and
\code{Turnc.alpha}.
}
}

\value{
Function \code{dpmj} returns the following:
\item{call}{the matched call.}

\item{seed}{the seed that was used (in case replication of the results is needed).}

\item{meanReg}{if \code{Xpred} is specified, the function returns the posterior mean of the conditional expectation of the response \eqn{y} 
given each new covariate \eqn{x}.}

\item{medianReg}{if \code{Xpred} is specified, the function returns the posterior mean of the conditional 50\% quantile of the response \eqn{y} 
given each new covariate \eqn{x}.}

\item{q1Reg}{if \code{Xpred} is specified, the function returns the posterior mean of the conditional 25\% quantile of the response \eqn{y} 
given each new covariate \eqn{x}.}

\item{q3Reg}{if \code{Xpred} is specified, the function returns the posterior mean of the conditional 75\% quantile of the response \eqn{y} 
given each new covariate \eqn{x}.}

\item{modeReg}{if \code{Xpred} is specified, the function returns the posterior mean of the conditional mode of the response \eqn{y} 
given each new covariate \eqn{x}.}

\item{denReg}{if \code{Xpred} is specified, the function returns the posterior mean conditional density of the response \eqn{y} 
given each new covariate \eqn{x}. Results are presented in a matrix the rows of which correspond to the different \eqn{x}s.}

\item{denVar}{if \code{Xpred} is specified, the function returns the posterior variance of the conditional density of the response \eqn{y} 
given each new covariate \eqn{x}. Results are presented in a matrix the rows of which correspond to the different \eqn{x}s.}

Further, function \code{dpmj} creates files where the posterior samples are written. These files are (with all file names
preceded by `BNSP.'):

\item{alpha.txt}{this file contains samples from the posterior of the concentration parameters  \eqn{\alpha}. 
The file is arranged in \code{(sweeps-burn)/thin} lines and one column, each line including one posterior sample.}

\item{compAlloc.txt}{this file contains the allocations to clusters obtained during posterior sampling. 
It consists of \code{(sweeps-burn)/thin} lines, that represent the posterior samples, and \eqn{n} columns, that represent the sampling units. 
Clusters are represented by integers ranging from 0 to \code{ncomp}-1.}

\item{MeanReg.txt}{this file contains the conditional means of the response \eqn{y} given covariates \eqn{x} obtained during posterior sampling. 
The rows represent the \code{(sweeps-burn)/thin} posterior samples. The columns represent the various covariate values \eqn{x} for which the 
means are obtained.}

\item{MedianReg.txt}{this file contains the 50\% conditional quantile of the response \eqn{y} given covariates \eqn{x} obtained 
during posterior sampling. The rows represent the \code{(sweeps-burn)/thin} posterior samples. The columns represent the various covariate 
values \eqn{x} for which the medians are obtained.}

\item{muh.txt}{this file contains samples from the posteriors of the \eqn{p}-dimensional mean vectors  \eqn{\mu_h, h=1,2,\dots},\code{ncomp}. 
The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{p} columns. In more detail, sweeps create \code{ncomp} 
lines representing samples \eqn{\mu_h^{(sw)}, h=1,\dots,}\code{ncomp}, where superscript \eqn{sw} represents a particular sweep. 
The elements of \eqn{\mu_h^{(sw)}} are written in the columns of the file.}

\item{nmembers.txt}{this file contains \code{(sweeps-burn)/thin} lines and \code{ncomp} columns, where the lines represent posterior 
samples while the columns represent the components or clusters. The entries represent the number of sampling units allocated to each 
component.}

\item{Q05Reg.txt}{this file contains the 5\% conditional quantile of the response \eqn{y} given covariates \eqn{x} obtained during 
posterior sampling. The rows represent the \code{(sweeps-burn)/thin} posterior samples. The columns represent the various covariate values 
\eqn{x} for which the quantiles are obtained.}

\item{Q10Reg.txt}{as above, for the 10\% conditional quantile.}

\item{Q15Reg.txt}{as above, for the 15\% conditional quantile.}

\item{Q20Reg.txt}{as above, for the 20\% conditional quantile.}

\item{Q25Reg.txt}{as above, for the 25\% conditional quantile.}

\item{Q75Reg.txt}{as above, for the 75\% conditional quantile.}

\item{Q80Reg.txt}{as above, for the 80\% conditional quantile.}

\item{Q85Reg.txt}{as above, for the 85\% conditional quantile.}

\item{Q90Reg.txt}{as above, for the 90\% conditional quantile.}

\item{Q95Reg.txt}{as above, for the 95\% conditional quantile.}

\item{Sigmah.txt}{this file contains samples from the posteriors of the \eqn{q \times q} restricted covariance  
matrices \eqn{\Sigma_h^*, h=1,2,\dots,}\code{ncomp}. The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and \eqn{q^2} 
columns. In more detail, sweeps create \code{ncomp} lines representing samples \eqn{\Sigma_h^{(sw)}, h=1,\dots,}\code{ncomp}, where 
superscript \eqn{sw} represents a particular sweep. The elements of \eqn{\Sigma_h^{(sw)}} are written in the columns of the file.}

\item{xih.txt}{this file contains samples from the posteriors of parameters  \eqn{\xi_h}, \eqn{h=1,2,\dots,}\code{ncomp}. 
The file is arranged in \code{((sweeps-burn)/thin)*ncomp} lines and one or two columns, depending on the number of parameters in the 
selected \code{Fcdf}. Sweeps write in the file \code{ncomp} lines representing samples \eqn{\xi_h^{(sw)}, h=1,\dots,}\code{ncomp}, 
where superscript \eqn{sw} represents a particular sweep.}

\item{Updated.txt}{this file contains \code{(sweeps-burn)/thin} lines with the number of components updated at each iteration of the sampler
(relevant for slice sampling).}
}

\references{
Consul, P. C. & Famoye, G. C. (1992). Generalized Poisson regression model. Communications in 
Statistics - Theory and Methods, 1992, 89-109.

Papageorgiou, G. (2018). Bayesian density regression for discrete outcomes. arXiv:1603.09706v3 [stat.ME]. 

Papaspiliopoulos, O. (2008). A note on posterior sampling from Dirichlet mixture models. Technical report,
University of Warwick.

Sethuraman, J. (1994). A constructive definition of Dirichlet priors. Statistica Sinica, 4, 639-650.

Walker, S. G. (2007). Sampling the Dirichlet mixture model with slices. Communications in Statistics
Simulation and Computation, 36(1), 45-54.
}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
#Bayesian nonparametric joint model with binomial response Y and one predictor X
data(simD)
pred<-seq(with(simD,min(X))+0.1,with(simD,max(X))-0.1,length.out=30)
npred<-length(pred)
# fit1 and fit2 define the same model but with different numbers of
# components and posterior samples
fit1 <- dpmj(cbind(Y,(E-Y))~X, Fcdf="binomial", data=simD, ncomp=30, sweeps=120,
             burn=100, sampler="truncated", Xpred=pred, offsetPred=30)
\donttest{fit2 <- dpmj(cbind(Y,(E-Y))~X, Fcdf="binomial", data=simD, ncomp=50, sweeps=5000,
               burn=1000, sampler="truncated", Xpred=pred, offsetPred=30)
plot(with(simD,X),with(simD,Y)/with(simD,E))
lines(pred,fit2$medianReg/30,col=3,lwd=2)}
}

\keyword{nonparametric}
\keyword{cluster}
