\name{als}
\alias{als}
\docType{data}
\title{
Example data set with a longitudinal variable
}
\description{
A simulated data set with a longitudinal variable and a survival variable
}
\usage{data("als")}
\format{
  A data frame with 312 observations on the following 6 variables.
  \describe{
    \item{\code{subjects}}{A subject indicator}
    \item{\code{trt}}{A treatment indicator}
    \item{\code{surv}}{The patient's survival time or last follow up time}
    \item{\code{csurv}}{An indicator that the patient died}
    \item{\code{times}}{The time that the longitudinal variable is observed}
    \item{\code{ALSFRS}}{The longitudinal variable}
  }
}
\keyword{datasets}
