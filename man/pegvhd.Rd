\name{pegvhd}
\docType{data}
\title{Plueral effusion and graft-verse-host disease after transplantation}
\alias{pegvhd}
\usage{pegvhd}
\description{
  This hypothetical data set is the study for the plueral effusion and graft-verse-host disease (GVHD) after allogenic transplantation. This data set is provided to illustrate the analysis with time-to-event data as covariates.
}
\format{This data frame consists of 50 individuals (25 with plueral effusion and 25 without plueral effusion) and the following columns:
\itemize{
  \item{\code{os.t}}{Time from the date of transplantation to the date of death. Days}
  \item{\code{os.s}}{Status. 1 = death and 0 = censored}
  \item{\code{gvhd.t}}{Time from the date of transplantation to the date of GVHD (either acute GVHD or chronic GVHD, whichever is first). Days}
  \item{\code{gvhd.s}}{Status. 1 = GVHD and 0 = censored}
  \item{\code{pe.t}}{Time from the date of transplantation to the date of plueral effusion. Days}
  \item{\code{pe.s}}{Status. 1 = plueral effusion and 0 = censored}
  \item{\code{age}}{Age}
  \item{\code{sex}}{Gender}
}
}
\references{
Modi D, Jang H, Kim S, Deol A, Ayash L, Bhutani D, Lum LG, Ratanatharathorn V, Manasa R, Mellert K, Uberti JP. (2016) Incidence, etiology, and outcome of pleural effusions in allogeneic hematopoietic stem cell transplantation. \emph{American Journal of Hematology}, in press.
}
\keyword{datasets}
