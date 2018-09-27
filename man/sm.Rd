\name{sm}
\alias{sm}

\title{Smooth terms in mvrm formulae}

\description{Function used to define smooth effects in the mean and variance formulae of function \code{mvrm}. 
The function is used internally to construct the design matrices.}

\usage{sm(...,k=10,knots=NULL,bs="rd")}

\arguments{
\item{...}{one or two covariates that the smooth term is a function of. If two covariates are used, 
they may be both continuous or one continuous and one discrete. Discrete variables should be defined as \code{factor}
in the \code{data} argument of the calling \code{mvrm} function.}
\item{k}{the number of knots to be utilized in the basis function expansion.}
\item{knots}{the knots to be utilized in the basis function expansion.}
\item{bs}{a two letter character indicating the basis functions to be used. Currently, the options are
\code{"rd"} that specifies radial basis functions and is available for univariate and bivariate smooths,
and \code{"pl"} that specifies thin plate splines that are available for univariate smooths.}
}

\details{Use this function within calls to function \code{mvrm} to specify smooth terms in the mean and/or variance function
of the regression model.

Univariate radial basis functions with \eqn{q} basis functions or \eqn{q-1} knots are defined by
\deqn{
\mathcal{B}_1 = \left\{\phi_{1}(u)=u , \phi_{2}(u)=||u-\xi_{1}||^2 \log\left(||u-\xi_{1}||^2\right), \dots, 
\phi_{q}(u)=||u-\xi_{q-1}||^2 \log\left(||u-\xi_{q-1}||^2\right)\right\},
}
where \eqn{||u||} denotes the Euclidean norm of \eqn{u} and \eqn{\xi_1,\dots,\xi_{q-1}} are the knots that  
are chosen as the quantiles of the observed values of explanatory variable \eqn{u}, 
with \eqn{\xi_1=\min(u_i), \xi_{q-1}=\max(u_i)} and the remaining knots chosen as equally spaced quantiles between
\eqn{\xi_1} and \eqn{\xi_{q-1}}. 


Thin plate splines are defined by 
\deqn{
\mathcal{B}_2 = \left\{\phi_{1}(u)=u , \phi_{2}(u)=(u-\xi_{1})_{+}, \dots, \phi_{q}(u)=(u-\xi_{q})_{+}\right\},
}
where \eqn{(a)_+ = \max(a,0)}. 
 
Radial basis functions for bivariate smooths are defined by
\deqn{
\mathcal{B}_3 = \left\{u_1,u_2,\phi_{3}(u)=||u-\xi_{1}||^2 \log\left(||u-\xi_{1}||^2\right), \dots, 
\phi_{q}(u)=||u-\xi_{q-1}||^2 \log\left(||u-\xi_{q-1}||^2\right)\right\}.
} 
}

\value{Specifies the design matrices of an \code{mvrm} call}


\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{#see \code{mvrm} example}

\seealso{\code{\link{mvrm}}}

\keyword{models} 
\keyword{smooth} 
\keyword{regression}
