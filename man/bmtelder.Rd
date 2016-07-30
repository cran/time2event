\name{bmtelder}
\docType{data}
\title{Bone marror transplant for old patients}
\alias{bmtelder}
\usage{bmtelder}
\description{
This hypothetical data set gives the performance of allogenic bone marror transplantation for old patients. This data set is provided to illustrate the analysis with time-to-event data as covariates, in particular, competing risk regression.
}
\format{This data frame consists of 50 individuals and the following columns:
\itemize{
  \item{\code{os.t}}{Time from the date of transplantation to the date of death. Days}
  \item{\code{os.s}}{Status. 1 = death and 0 = censored}
  \item{\code{nrm.t}}{Time from the date of transplantation to the date of death without relapse, with relapse as competing risk. Days}
  \item{\code{nrm.s}}{Status. 2 = relapse, 1 = death, and 0 = censored}
  \item{\code{cgvhd.t}}{Time from the date of transplantation to the date of chronic graft-verse-host disease (GVHD). Days}
  \item{\code{cgvhd.s}}{Status. 1 = chronic GVHD and 0 = censored}
  \item{\code{cond}}{Conditioning regimen. Reduced intensity and full intensity}
  \item{\code{donor}}{donor to host matching. Allo Related (allogeneic related) and Allo Unrelated (allogeneic unrelated)}
}
}
\keyword{datasets}
