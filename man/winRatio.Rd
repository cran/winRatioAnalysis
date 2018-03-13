\name{winRatio}
\alias{winRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyze the winRatio as a function of time}
\description{Calulates the win ratio of the data and additional fits a model to data separately for each treatment group and then calculates the win-Ratio as a function of follow-up time. The data must be one of three types, progression and survival data where progression is interval censored and progression always proceeds survival, dual time to event data where both events have a common censoring time, and longitudinal and survival data where the longitudinal variable follows a random slopes model and the survival data has a frailty that is equal to the value of $u+b t$ where $u$ is the random intercept and $b$ is the random slope. }
\usage{
winRatio(survivalObject,treatmentVariable,treatmentCodes,data,
                method=c('pssm','NP'),
                plotTimeUnit=NULL,
                secondSurvivalObject=NULL,
                longitudinalVariable=NULL,
                   timeVar=NULL,
                   subjectId=NULL,
                   plotPoints=NULL,
                   nominalTimes=NULL,
                   pssmIntervals=3,
                   integrationIntervals=1,
                   type='l',
                   mult=100,
                   findValue=function(times,values,c)
                     stats::approx(times,values,xout=c,method="constant",
                            rule=2,f=0)$y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survivalObject}{An survival object for the primary survival variable of the form \code{Surv(time,indicator)}}

  \item{treatmentVariable}{A character variable containing the name of the treatment variable}
  \item{treatmentCodes}{ A vector containing the treatment codes, control followed by active
}
 \item{data}{A data frame that the variable names refer to}
 \item{method}{Method of model fitting one of 'pssm','NP', if the data is longitudinal this is ignored and a simulation is used}
\item{plotTimeUnit}{If plots are desired enter the plot time unit eg. \code{'days','months','years'},
if not present no plots will be produced}
\item{secondSurvivalObject}{A survival object for the secondary time to event variable of the form \code{Surv(time,indicator)},
  or \code{Surv(time0,time1,type='interval')}}

  \item{longitudinalVariable}{A character variable containing the name of the longitudinal variable if present}
  \item{timeVar}{A character variable containing the name of the time variable in the longitudinal model}
  \item{subjectId}{A character variable containing the name of the subject identifier}

  \item{plotPoints}{A numeric vector of follow up times to calculate the Win Ratio estimator and related statistics, if NULL the follow up in the data is used }

  \item{nominalTimes}{The times that the longitudinal variable is measured. Note the data may not contain these exact times but the simulation used to calculate the win ratio over time will generate data for these times}
  \item{pssmIntervals}{How many intervals to used if pssm is used to calculate the joint model}
  \item{integrationIntervals}{Interval to use for trapozoidal rule}
  \item{type}{If NP is used whether to use the upper or lower probability}
  \item{mult}{The sample size multiplier for the simulation when their is longitudinal data. \code{mult=100} would produce a sample size that was 100 times that of the orginal sample.}
\item{findValue}{A function of a numeric vector of times, a numeric vector of corresponding values of the longitudinal variable at those times and a scalar c, considered a time. The function interpolates the values to give a value at time c.  The default is last value carried forward.}

}
\value{
A winratio object which is implemented as a data frame with columms:

\item{time }{The follow-up times, The first row it is the median follow up in the data,
subsequent rows are the times in \code{plotPoints}}
\item{firstP10}{The probability that the treatment group survival is better than the control group survival}
\item{firstP10}{The probability that the control group survival is better than the treatment group survival}
\item{secondP10}{The probability that the treatment group secondary endpoint is better than the control group secondary endpoint}
\item{secondP01}{The probability that the control group secondary endpoint is better than the treatment group secondary endpoint}
\item{winRatio}{The win ratio or probability that a treatment group does better than then a control group patient divided by the converse}
\item{winRatioSE}{Standard Error of the Win Ratio}
}
\references{Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a composite outcome based on prioritized components. Biostatistics. 2015 Sep 8;17(1):178-87.
}
\author{David A. Schoenfeld
}
\note{
Code for calculating the standard error of the win ratio was contributed by
Ionut Bebu, Ph.D. and Lohn Lachin Ph.D.
}


\examples{

winRatio(Surv(tdeath,cdeath),
treatmentVariable='rx', treatmentCodes = c(1, 0),data=dat,
secondSurvivalObject=Surv(tprog0,tprog1,type='interval'),pssmIntervals=1,
method = "pssm", plotPoints =3, integrationIntervals=1)

winRatio(Surv(tdeath,cdeath),
treatmentVariable='rx', treatmentCodes = c(1, 0),data=dat,
secondSurvivalObject=Surv(tprogression,cprog),
method = "NP", plotPoints =3)


\dontrun{
winRatio(Surv(surv,csurv),
longitudinalVariable = "ALSFRS",
timeVar="times",subject="subjects",
treatmentVariable='trt',
treatmentCodes = c(0, 1),data=als,method = "NP",
nominalTimes=c(0,6,12),
plotPoints =6,mult=1)
#          time firstP10 firstP01 secondP10 secondP01 winRatio WinRatioSE
#Value 11.93747  0.24500 0.313125  0.199375  0.242500 1.250352  0.3007245
#plots  6.00100  0.49375 0.468750  0.160625  0.200625 1.022923         NA
}
}


\keyword{htest}
\keyword{survival}
\keyword{nonparametric}
