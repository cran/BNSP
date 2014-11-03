\name{dnorm.pois}
\alias{dnorm.pois}

\title{Bivariate Normal-Poisson distribution}

\description{Normal-Poisson probability density function.}

\usage{
dnorm.pois(x,y,mu,Sigma,rate,E)
}

\arguments{
\item{x}{Real valued sample from a univariate Normal distribution.}

\item{y}{Integer valued sample from a univariate Poisson distribution.}

\item{mu}{Two-dimentional vector with means, denoted as \eqn{\mu} below.}

\item{Sigma}{Two-by-two covariance matrix, denoted as \eqn{\Sigma} below.}

\item{rate}{Mean of the Poisson variable, denoted as \eqn{\lambda} below. }

\item{E}{Offset term, denoted as \eqn{E} below.}

}

\value{Function \code{dnorm.pois} returns the joint probability density function of correlated Normal and Poisson 
random variables
\deqn{
f(x,y|\theta) = \int_{c_{y-1}}^{c_{y}} N(x,y^{*}|\mu,\Sigma) dy^{*},}
where \eqn{y^*} denotes a continuous random variable that determines the count according to the rule
\deqn{Y = y \iff c_{y-1} < y^* < c_{y}.}
Cut-points \eqn{c_y} are defined by
\deqn{c_{y} = \Phi^{-1}(F(y;E\lambda)),}
where \eqn{\Phi()} is the Normal cdf, \eqn{F()} the Poisson cdf, \eqn{E} denotes an offset term 
and \eqn{\lambda} is the mean of the Poisson.
Further, \eqn{\mu} and \eqn{\Sigma} denote the mean vector and covariance matrix of \eqn{(x,y^{*})}.

The integral is evaluated using the univariate Normal cdf
\deqn{
f(y,x|\theta) = 
N(x|\mu_{1},\sqrt\Sigma_{11})\{\Phi(\frac{c_y-E(y^*|x)}{sd(y^*|x)})-\Phi(\frac{c_{y-1}-E(y^*|x)}{sd(y^*|x)})\},}
where \eqn{E(y^*|x) = \mu_{2}+\Sigma_{12}(x-\mu_1)/\Sigma_{11}} and
\eqn{sd(y^*|x) = \sqrt{\Sigma_{22}-\Sigma_{12}^2 / \Sigma_{11}}}.
}

\note{The mean \eqn{\mu_2} and variannce \eqn{\Sigma_{22}} of \eqn{y^*} are usually set to zero and one.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
#When the covariance matrix is diagonal dnorm.pois is equal to the product of dnorm and dpois
mu<-c(0,0)
cov.mat<-matrix(c(1,0.0,0.0,1),ncol=2,nrow=2)
dnorm.pois(0,5,mu=mu,Sigma=cov.mat,rate=3,E=2)
dnorm(0,0,1)*dpois(5,6)
#Otherwise not equal
mu<-c(0,0)
cov.mat<-matrix(c(1,-0.8,-0.8,1),ncol=2,nrow=2)
dnorm.pois(0,5,mu=mu,Sigma=cov.mat,rate=3,E=2)}

\keyword{distribution}
