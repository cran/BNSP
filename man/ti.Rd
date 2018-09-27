\name{ti}
\alias{ti}

\title{mgcv constructor \code{ti}}

\description{Provides interface between mgcv::ti and BNSP. \code{ti(...)} calls 
\code{mgcv::smoothCon(mgcv::ti(...),...}}

\usage{ti(..., data, knots = NULL, absorb.cons = FALSE, scale.penalty = TRUE, 
n = nrow(data), dataX = NULL, null.space.penalty = FALSE, sparse.cons = 0, 
diagonal.penalty = FALSE, apply.by = TRUE, modCon = 0, k = NA, bs = "cr", 
m = NA, d = NA, by = NA, fx = FALSE, np = TRUE, xt = NULL, id = NULL, 
sp = NULL, mc = NULL, pc = NULL)}

\arguments{
\item{...}{a list of variables. See mgcv::ti}
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
\item{k}{see mgcv::ti}
\item{bs}{see mgcv::ti}
\item{m}{see mgcv::ti}
\item{d}{see mgcv::ti}
\item{by}{see mgcv::ti}
\item{fx}{see mgcv::ti}
\item{np}{see mgcv::ti}
\item{xt}{see mgcv::ti}
\item{id}{see mgcv::ti}
\item{sp}{see mgcv::ti}
\item{mc}{see mgcv::ti}
\item{pc}{see mgcv::ti}
}

\details{The most relevant arguments for BNSP users are the list of variables \code{...}, \code{knots}, \code{absorb.cons}, \code{bs}, and \code{by}.}

\value{A design matrix that specifies a smooth term in a model.}

\author{Georgios Papageorgiou \email{gpapageo@gmail.com}}
