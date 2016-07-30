\name{tcomp.risk}
\alias{tcomp.risk}
\title{
Competings Risks Regression with time-to-event data as covariates.
}
\description{
Fits a semiparametric model for the cause-specific quantitie with time-to-event data as covariates.
}
\usage{
tcomp.risk(formula, na.time=c("remove","censor"), verbose=FALSE,
      data = sys.parent(), cause, times = NULL,
    	Nit = 50, clusters = NULL, est = NULL, fix.gamma = 0, gamma = 0,
    	n.sim = 0, weighted = 0, model = "fg", detail = 0, interval = 0.01,
    	resample.iid = 1, cens.model = "KM", cens.formula = NULL,
    	time.pow = NULL, time.pow.test = NULL, silent = 1, conv = 1e-06,
    	weights = NULL, max.clust = 1000, n.times = 50, first.time.p = 0.05,
    	estimator = 1, trunc.p = NULL, cens.weights = NULL, admin.cens = NULL,
    	conservative = 1, monotone = 0, step = NULL)
}
\arguments{
  \item{formula}{
a formula object, with the response on the left of a '~' operator, and the terms on the right. The response must be a survival object as returned by the `Event' function. The status indicator is not important here. Time-invariant regressors are specified by the wrapper \code{const()}, and cluster variables (for computing robust variances) by the wrapper \code{cluster()}. In case that time-to-event data are covariates, use the wrapper \code{time()} to indicate the time-to-event data as covariates.
}
  \item{na.time}{
a missing-data filter function for time-to-event covariates. The option 'remove' will remove all the data with 'NA', while the option 'censor' will treat the missing data as censored and then replace with the logest time. Default is 'remove'.
}
  \item{verbose}{
logical. Should R report extra information on progress? Default is 'FALSE'.
}
  \item{data}{
a data.frame with the variables.
}
  \item{cause}{
For competing risk models specificies which cause we consider.
}
  \item{times}{
specifies the times at which the estimator is     considered. Defaults to all the times where an event of interest occurs, with the first 10 percent or max 20 jump points removed for numerical stability in simulations.
}
  \item{Nit}{
number of iterations for Newton-Raphson algorithm.
}
  \item{clusters}{
specifies cluster structure, for backwards compability.
}
  \item{est}{
possible starting value for nonparametric component of model.
}
  \item{fix.gamma}{
to keep gamma fixed, possibly at 0.
}
  \item{gamma}{
starting value for constant effects.
}
  \item{n.sim}{
number of simulations in resampling.
}
  \item{weighted}{
Not implemented. To compute a variance weighted version of the test-processes used for testing time-varying effects.
}
  \item{model}{
"additive", "prop"ortional, "rcif", or "logistic".
}
  \item{detail}{
if 0 no details are printed during iterations, if 1 details are given.
}
  \item{interval}{
specifies that we only consider timepoints where the Kaplan-Meier of the censoring distribution is larger than this value.
}
  \item{resample.iid}{
to return the iid decomposition, that can be used to construct confidence bands for predictions
}
  \item{cens.model}{
specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"
}
  \item{cens.formula}{
specifies the regression terms used for the regression model for chosen regression model. When cens.model is specified, the default is to use the same design as specified for the competing risks model.
}
  \item{time.pow}{
specifies that the power at which the time-arguments is transformed, for each of the arguments of the \code{const()} terms, default is 1 for the additive model and 0 for the proportional model.
}
  \item{time.pow.test}{
specifies that the power the time-arguments is transformed for each of the arguments of the \code{non-const()} terms. This is relevant for testing if a coefficient function is consistent with the specified form \eqn{A_l(t)=beta_l t^time.pow.test(l)}. Default is 1 for the additive model and 0 for the proportional model.
}
  \item{silent}{
if 0 information on convergence problems due to non-invertible derviates of scores are printed.
}
  \item{conv}{
gives convergence criterie in terms of sum of absolute change of parameters of model
}
  \item{weights}{
weights for estimating equations.
}
  \item{max.clust}{
sets the total number of i.i.d. terms in i.i.d. decompostition. This can limit the amount of memory used by coarsening the clusters. When NULL then all clusters are used. Default is 1000 to save memory and time.
}
  \item{first.time.p}{
first point for estimation is pth percentile of cause jump times.
}
  \item{n.times}{
only uses 50 points for estimation, if NULL then uses all points, subject to p.start condition.
}
  \item{estimator}{
default estimator is 1.
}
  \item{trunc.p}{
truncation weight for delayed entry, P(T > entry.time | Z_i), typically Cox model.
}
  \item{cens.weights}{
censoring weights can be given here rather than calculated using the KM, cox or aalen models.
}
  \item{admin.cens}{
censoring times for the administrative censoring
}
  \item{conservative}{
set to 0 to compute correct variances based on censoring weights, default is conservative estimates that are much quicker.
}
  \item{monotone}{
monotone=0, uses estimating equations
montone 1 uses
}
  \item{step}{
step size for Fisher-Scoring algorithm.
}
}
\details{
The funciton \code{tcomp.risk} is an extention of the function \code{comp.risk} for time-to-event covariates. If the model has no time-to-event covariates, \code{tcomp.risk} will print the warning sign 'No time-varying covariate!!!' and then do exactly the same procedure as \code{comp.risk} does. If the model has time-to-event covariates, the time-to-event covaraites should be wrapped with \code{time()} by placing the right-hand side of a ~ operator. In particular, the wrapper \code{time(a1,b1,a2,b2,a3,b3,...)} will be used with time-to-event covariates, where \eqn{ai} and \eqn{bi}, \eqn{i=1,2,...} are time-to-event and status, respectively. See \code{comp.risk} for other details.
}
\value{
returns the same object as that of \code{comp.risk()}. See \code{comp.risk()} for details
}
\references{
S. Kim (2016). time2event: an R package for the analysis of event time data with time-to-event data as covariates. Wayne State University/Karmanos Cancer Institute. Manuscript.
}
\author{
Seongho Kim
}
\examples{
  data(bmtelder)

  # convert to data with time-to-event data as covariates
  # nrm with cgvhd
  tnrm2data = time2data(c("nrm.t","nrm.s"),c("cgvhd.t","cgvhd.s"),bmtelder)$data

  # no time-varying analysis with 'comp.risk'
  set.seed(3927)
  cr2r = comp.risk(Event(nrm.t,nrm.s)~cgvhd.s+cond+donor,data=bmtelder,
  			cause=1,resample.iid=1,n.sim=1000,model="additive")
  cr2r.pred = predict(cr2r,X=1)
  plot(cr2r.pred)

  # time-varying analysis with 'comp.risk'
  set.seed(3927)
  nt.cr2r = comp.risk(Event(start,end,nrm.s)~cgvhd.s+cond+donor,data=tnrm2data,
  			cause=1,resample.iid=1,n.sim=1000,model="additive")
  nt.cr2r.pred = predict(nt.cr2r,X=1)
  plot(nt.cr2r.pred)

  # time-varying analysis with 'tcomp.risk'
  set.seed(3927)
  t.cr2r = tcomp.risk(Event(nrm.t,nrm.s)~time(cgvhd.t,cgvhd.s)+cond+donor,data=bmtelder,
  			cause=1,resample.iid=1,n.sim=1000,model="additive")
  t.cr2r.pred = predict(t.cr2r,X=1)
  plot(t.cr2r.pred)
}

