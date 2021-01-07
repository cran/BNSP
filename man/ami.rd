\name{ami}
\docType{data}
\alias{ami}

\title{Amitriptyline dataset from Johnson and Wichern}

\description{Amitriptyline is a prescription antidepressant. The dataset consists of measurements on 17 patients who had over-dosed on amitriptyline.}

\usage{data(ami)}

\format{
  A data frame containing 17 rows and 7 columns. The columns represent
    \describe{
    \item{\code{tot}}{total blood plasma level.}
    \item{\code{ami}}{amount of amitriptyline found in the plasma.}
    \item{\code{gen}}{gender (1 for female).}
    \item{\code{amt}}{amount of the drug taken.}
    \item{\code{pr}}{PR wave measurement.}
    \item{\code{bp}}{diastolic blood pressure.}
    \item{\code{qrs}}{QRS wave measurement.}
  }
}

\source{Johnson, R. A., and Wichern, D. W. (2007), Applied Multivariate Statistical Analysis, Essex: Pearson, page 426.}

\references{
Johnson, R. A., and Wichern, D. W. (2007). \emph{Applied Multivariate Statistical Analysis}, Essex: Pearson.}

\keyword{datasets}
