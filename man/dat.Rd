\name{dat}
\alias{dat}
\docType{data}
\title{
Example dataset with progression and survival
}
\description{
This dataset provides a simulated example for the pssm and NP method of calculating the winRatio and it's components over time
}
\usage{data("dat")}
\format{
  A data frame with 200 observations on the following 7 variables.
  \describe{
    \item{\code{tprog0}}{Last time the patient was observed to be free of progression}
    \item{\code{tprog1}}{The first time the patient was observed to have progressed}
    \item{\code{cdeath}}{An indicator that the patient has died 1=death}
    \item{\code{tdeath}}{The time of death or last follow up time}
    \item{\code{rx}}{Treatment indicator(0 or 1)}
    \item{\code{cprog}}{An indicator of progression}
    \item{\code{tprogression}}{Exact time of progression or last follow up time}
  }
}
\keyword{datasets}
