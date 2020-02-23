\name{chol}
\alias{chol}

\title{The Cholesky and modified Cholesky decompositions}

\description{
Computes the Cholesky factorization and modified Cholesky factorizations of a real symmetric positive-definite square matrix.
}

\usage{chol(x, mod = TRUE, p = 1, ...)}

\arguments{
\item{x}{A symmetric, positive-definite matrix.}
\item{mod}{Defaults to TRUE. With this choice, the function returns the modified Cholesky decomposition. When mod = FALSE, the function returns the usual Cholesky
           decomposition.}
\item{p}{Relevant only when \code{mod = TRUE}. It determines the size of the blocks of the block diagonal matrix.}
\item{...}{other arguments.}
}

\details{The function computes the modified Cholesky decomposition of a real symmetric positive-definite square matrix \eqn{\Sigma}. This is given by
         \deqn{L \Sigma L^{\top} = D,}
         where \eqn{L} is a lower tringular matrix with ones on its main diagonal and D is a block diagonal matrix with block size determined by argument \code{p}.}

\value{The function returns matrices \eqn{L} and \eqn{D}.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}

\examples{
Sigma <- matrix(c(1.21,0.18,0.13,0.41,0.06,0.23,
                  0.18,0.64,0.10,-0.16,0.23,0.07,
                  0.13,0.10,0.36,-0.10,0.03,0.18,
                  0.41,-0.16,-0.10,1.05,-0.29,-0.08,
                  0.06,0.23,0.03,-0.29,1.71,-0.10,
                  0.23,0.07,0.18,-0.08,-0.10,0.36),6,6)
LD <- chol(Sigma)
L <- LD$L
D <- LD$D
round(L,5)
round(D,5)
solve(L) \%*\% D \%*\% solve(t(L))
LD <- chol(Sigma, p = 2)
L <- LD$L
D <- LD$D
round(L, 5)
round(D, 5)
solve(L) \%*\% D \%*\% solve(t(L))
}

\seealso{The default function from \pkg{base}, \code{\link[base]{chol}}}
