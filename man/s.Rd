\name{s}
\alias{s}

\title{mgcv constructor \code{s}}

\description{Provides interface between mgcv::s and BNSP. \code{s(...)} calls 
\code{mgcv::smoothCon(mgcv::s(...),...}}

\usage{s(..., data, knots = NULL, absorb.cons = FALSE, scale.penalty = TRUE, 
n = nrow(data), dataX = NULL, null.space.penalty = FALSE, sparse.cons = 0, 
diagonal.penalty = FALSE, apply.by = TRUE, modCon = 0, k = -1, fx = FALSE, 
bs = "tp", m = NA, by = NA, xt = NULL, id = NULL, sp = NULL, pc = NULL)}

\arguments{
\item{...}{a list of variables. See mgcv::s}
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
\item{k}{see mgcv::s}
\item{fx}{see mgcv::s}
\item{bs}{see mgcv::s}
\item{m}{see mgcv::s}
\item{by}{see mgcv::s}
\item{xt}{see mgcv::s}
\item{id}{see mgcv::s}
\item{sp}{see mgcv::s}
\item{pc}{see mgcv::s}
}

\details{The most relevant arguments for BNSP users are the list of variables \code{...}, \code{knots}, \code{absorb.cons}, \code{bs}, and \code{by}.}

\value{A design matrix that specifies a smooth term in a model.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}
