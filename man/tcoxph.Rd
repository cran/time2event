\name{tcoxph}
\alias{tcoxph}
\title{
Fit Proportional Hazards Regression Model with time-to-event data as covariates.
}
\description{
Fits a Cox proportional hazards regression model with time-to-event data as covariates.
}
\usage{
tcoxph(formula, na.time=c("remove","censor"),verbose=FALSE,
                     data, weights, subset, na.action, init, control,
                     ties = c("efron", "breslow", "exact"), singular.ok = TRUE,
                     robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties,
                     ...)
}
\arguments{
  \item{formula}{
a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the Surv function. In case that time-to-event data are covariates, use the wrapper \code{time()} to indicate the time-to-event data as covariates.
}
  \item{na.time}{
a missing-data filter function for time-to-event covariates. The option 'remove' will remove all the data with 'NA', while the option 'censor' will treat the missing data as censored and then replace with the logest time. Default is 'remove'.
}
  \item{verbose}{
logical. Should R report extra information on progress? Default is 'FALSE'.
}
  \item{data}{
a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument.
}
  \item{weights}{
vector of case weights. For a thorough discussion of these see the book by Therneau and Grambsch.
}
  \item{subset}{
expression indicating which subset of the rows of data should be used in the fit. All observations are included by default.
}
  \item{na.action}{
a missing-data filter function. This is applied to the model.frame after any subset argument has been used. Default is \code{options()}\$na.action.
}
  \item{init}{
vector of initial values of the iteration. Default initial value is zero for all variables.
}
  \item{control}{
Object of class coxph.control specifying iteration limit and other control options. Default is \code{coxph.control(...)}.
}
  \item{ties}{
a character string specifying the method for tie handling. If there are no tied death times all the methods are equivalent. Nearly all Cox regression programs use the Breslow method by default, but not this one. The Efron approximation is used as the default here, it is more accurate when dealing with tied death times, and is as efficient computationally. The ``exact partial likelihood'' is equivalent to a conditional logistic model, and is appropriate when the times are a small set of discrete values. See further below.
}
  \item{singular.ok}{
logical value indicating how to handle collinearity in the model matrix. If TRUE, the program will automatically skip over columns of the X matrix that are linear combinations of earlier columns. In this case the coefficients for such columns will be NA, and the variance matrix will contain zeros. For ancillary calculations, such as the linear predictor, the missing coefficients are treated as zeros.
}
  \item{robust}{
this argument has been deprecated, use a cluster term in the model instead. (The two options accomplish the same goal -- creation of a robust variance -- but the second is more flexible).
}
  \item{model}{
logical value: if TRUE, the model frame is returned in component model.
}
  \item{x}{
logical value: if TRUE, the x matrix is returned in component x.
}
  \item{y}{
logical value: if TRUE, the response vector is returned in component y.
}
  \item{tt}{
optional list of time-transform functions.
}
  \item{method}{
alternate name for the ties argument.
}
  \item{...}{
Other arguments will be passed to \code{coxph.control}
}
}
\details{
The funciton \code{tcoxph} is an extention of the function \code{\link{coxph}} for time-to-event covariates. If the model has no time-to-event covariates, \code{tcoxph} will print the warning sign 'No time-varying covariate!!!' and then do exactly the same procedure as \code{coxph} does. If the model has time-to-event covariates, the time-to-event covaraites should be wrapped with \code{time()} by placing the right-hand side of a ~ operator. In particular, the wrapper \code{time(a1,b1,a2,b2,a3,b3,...)} will be used with time-to-event covariates, where \eqn{ai} and \eqn{bi}, \eqn{i=1,2,...} are time-to-event and status, respectively. See \code{\link{coxph}} for other details.
}
\value{
the same object as that of class \code{\link{coxph}} representing the fit. See \code{\link{coxph.object}} for details.
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

  # no time-varying analysis with 'coxph' and 'comp.risk'
  os1r = coxph(Surv(os.t,os.s)~gvhd.s+pe.s+age+sex,data=pegvhd)

  # time-varying analysis with 'coxph' and 'comp.risk'
  nt.os1r = coxph(Surv(start,end,os.s)~gvhd.s+pe.s+age+sex,data=tos1data)

  # time-varying analysis with 'tcoxph' and 'tcomp.risk'
  t.os1r = tcoxph(Surv(os.t,os.s)~time(gvhd.t,gvhd.s)+time(pe.t,pe.s)+age+sex
  			,data=pegvhd)

  data(bmtelder)

  # convert to data with time-to-event data as covariates
  # os with cgvhd
  tos2data = time2data(c("os.t","os.s"),c("cgvhd.t","cgvhd.s"),bmtelder)$data

  # no time-varying analysis with 'coxph'
  os2r = coxph(Surv(os.t,os.s)~cgvhd.s+cond+donor,data=bmtelder)

  # time-varying analysis with 'coxph'
  nt.os2r = coxph(Surv(start,end,os.s)~cgvhd.s+cond+donor,data=tos2data)

  # time-varying analysis with 'tcoxph'
  t.os2r = tcoxph(Surv(os.t,os.s)~time(cgvhd.t,cgvhd.s)+cond+donor,data=bmtelder)

}
