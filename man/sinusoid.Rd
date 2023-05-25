\name{sinusoid}
\alias{sinusoid}

\title{Sinusoid terms in mvrm formulae}

\description{Function used to define sinusoidal curves in the mean formula of function \code{mvrm}. 
The function is used internally to construct design matrices.}

\usage{sinusoid(..., harmonics = 1, amplitude = 1, period = 0, periodRange = NULL, 
breaks = NULL, knots = NULL)}

\arguments{
\item{...}{a single covariate that the sinusoid term is a function of.} 
\item{harmonics}{an integer value that denotes the number of sins and cosines to be utilized in the representation of a sinusoidal curve.}
\item{amplitude}{a positive integer. If set equal to one, it denotes a fixed amplitude. Otherwise, if set to an integer that is greater than one, it denotes the number of knots to be utilized in the representation of the time-varying amplitude.}
\item{period}{the period of the sinusoidal wave. Values less than or equal to zero signify that the period is unknown. 
Positive values signify that the period is known and fixed.}
\item{periodRange}{a vector of length two with the range of possible period values. It is required when the period is unknown.}
\item{breaks}{the growth break points.}
\item{knots}{the knots to be utilized in the representation of the time-varying amplitude. Relevant only when \code{amplitude} is greater than 1.}
}

\details{Use this function within calls to function \code{mvrm} to specify sinusoidal waves in the mean function
of a regression model.

Consider the sinusoidal curve 
\deqn{
y_t = \beta_0 + A(t) \sin(2\pi t/\omega+\varphi) + \epsilon_t,}
where \eqn{y_t} is the response at time \eqn{t}, \eqn{\beta_0} is an intercept term, \eqn{A(t)} is a time-varying amplitude, \eqn{\varphi \in [0,2\pi]} is the phase shift parameter, \eqn{\omega} is the period taken to be known, and
\eqn{\epsilon_t} is the error term.

The period \eqn{\omega} is defined by the argument \code{period}. 

The time-varying amplitude is represented using \eqn{A(t) = \sum_{j=1}^{K} \beta_{Aj} \phi_{Aj}(t)}, 
where \eqn{K}, the number of knots, is defined by argument \code{amplitude}. If \code{amplitude = 1}, then
the amplitude is taken to be fixed: \eqn{A(t)=A}.  

Further, \eqn{\sin(2\pi t/\omega+\varphi)} is represented utilizing
\eqn{\sin(2\pi t/\omega+\varphi) = \sum_{k=1}^{L} a_k \sin(2k\pi t/\omega) + b_k \cos(2k\pi t/\omega)},
where \eqn{L}, the number of harmonics, is defined by argument \code{harmonics}. 
}

\value{Specifies the design matrices of an \code{mvrm} call}


\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
# Simulate and fit a sinusoidal curve
# First load releveant packages
require(BNSP)
require(ggplot2)
require(gridExtra)
require(Formula)
# Simulate the data 
mu <- function(u) {cos(0.5 * u) * sin(2 * pi * u + 1)}
set.seed(1)
n <- 100
u <- sort(runif(n, min = 0, max = 2*pi))
y <- rnorm(n, mu(u), 0.1)
data <- data.frame(y, u)
# Define the model and call function \code{mvrm} that perfomes posterior sampling for the given 
# dataset and defined model 
model <- y ~ sinusoid(u, harmonics = 2, amplitude = 20, period = 1)
\dontrun{
m1 <- mvrm(formula = model, data = data, sweeps = 10000, burn = 5000, thin = 2, seed = 1, 
           StorageDir = getwd())
# Plot
x1 <- seq(min(u), max(u), length.out = 100)
plotOptionsM <- list(geom_line(aes_string(x = x1, y = mu(x1)), col = 2, alpha = 0.5, lty = 2),
                     geom_point(data = data, aes(x = u, y = y)))
plot(x = m1, term = 1, plotOptions = plotOptionsM, intercept = TRUE, 
     quantiles = c(0.005, 0.995), grid = 100, combine = 1)
}
}

\seealso{\code{\link{mvrm}}}

\keyword{models} 
\keyword{smooth} 
\keyword{regression}
