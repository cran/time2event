\name{time2data}
\alias{time2data}
\title{
Classification-based quantitative SILAC analysis
}
\description{
finding differentially expressed proteins using classification methods
}
\usage{
time2data(tvar,tcov,data,na.time=c("remove","censor"),verbose=FALSE,weights=NULL)
}
\arguments{
  \item{tvar}{
a vector of time-to-event response variables' names as a pair of 'time to event' and 'status'.
}
  \item{tcov}{
a vector of time-to-event covariates' names as a pair of 'time to event' and 'status'.
}
  \item{data}{
a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument. The variables in \code{tvar} and \code{tcov} must be part of data.
}
  \item{na.time}{
a missing-data filter function for time-to-event covariates. The option 'remove' will remove all the data with 'NA', while the option 'censor' will treat the missing data as censored and then replace with the logest time. Default is 'remove'.
}
  \item{verbose}{
logical. Should R report extra information on progress? Default is 'FALSE'.
}
  \item{weights}{
vector of case weights. For a thorough discussion of these see the book by Therneau and Grambsch. Default is 'NULL'.
}
}
\details{
It will restructure the data according to time-to-event covariates. After converting, the time-to-event response variable will be formulated as interval censored data and its structure will be 'start', 'end', and 'status'.
}
\value{
\item{data}{restructured data by time-to-event data as covariates. The time-to-event response variable will be structured as interval censored data as the form of 'start','end', and 'status'.}
\item{wt}{restructed weights}
}
\references{
S. Kim (2016). time2event: an R package for the analysis of event time data with time-to-event data as covariates. Wayne State University/Karmanos Cancer Institute. Manuscript.
}
\author{
Seongho Kim
}
\examples{
  data(pegvhd)

  # convert to data with time-to-event data as covariates
  # os with cgvhd
  tos1data = time2data(c("os.t","os.s"),c("gvhd.t","gvhd.s","pe.t","pe.s"),pegvhd)$data

  data(bmtelder)

  # convert to data with time-to-event data as covariates
  # os with cgvhd
  tos2data = time2data(c("os.t","os.s"),c("cgvhd.t","cgvhd.s"),bmtelder)$data

  # nrm with cgvhd
  tnrm2data = time2data(c("nrm.t","nrm.s"),c("cgvhd.t","cgvhd.s"),bmtelder)$data
}
