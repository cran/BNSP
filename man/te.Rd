\name{te}
\alias{te}

\title{mgcv constructor \code{te}}

\description{Provides interface between mgcv::te and BNSP. \code{te(...)} calls 
\code{mgcv::smoothCon(mgcv::te(...),...}}

\usage{te(..., data, knots = NULL, absorb.cons = FALSE, scale.penalty = TRUE, 
n = nrow(data), dataX = NULL, null.space.penalty = FALSE, sparse.cons = 0, 
diagonal.penalty = FALSE, apply.by = TRUE, modCon = 0, k = NA, bs = "cr", 
m = NA, d = NA, by = NA, fx = FALSE, np = TRUE, xt = NULL, id = NULL, 
sp = NULL, pc = NULL)}

\arguments{
\item{...}{a list of variables. See mgcv::te}
\item{data}{see mgcv::smoothCon}
\item{knots}{see mgcv::knots}
\item{absorb.cons}{see mgcv::smoothCon}
\item{scale.penalty}{see mgcv::smoothCon}
\item{n}{see mgcv::smoothCon}
\item{dataX}{see mgcv::smoothCon}
\item{null.space.penalty}{see mgcv::smoothCon}
\item{sparse.cons}{see mgcv::smoothCon}
\item{diagonal.penalty}{see mgcv::smoothCon}
\item{apply.by}{see mgcv::smoothCon}
\item{modCon}{see mgcv::smoothCon}
\item{k}{see mgcv::te}
\item{bs}{see mgcv::te}
\item{m}{see mgcv::te}
\item{d}{see mgcv::te}
\item{by}{see mgcv::te}
\item{fx}{see mgcv::te}
\item{np}{see mgcv::te}
\item{xt}{see mgcv::te}
\item{id}{see mgcv::te}
\item{sp}{see mgcv::te}
\item{pc}{see mgcv::te}
}

\details{The most relevant arguments for BNSP users are the list of variables \code{...}, \code{knots}, \code{absorb.cons}, \code{bs}, and \code{by}.}

\value{A design matrix that specifies a smooth term in a model.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}
