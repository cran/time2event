# time-varying covariates testing
time2data <- function(tvar,tcov,data,na.time=c("remove","censor"),verbose=FALSE,weights=NULL){

  cov.class <- function(var,iclass){
    cvar = as.character(var)
    if(iclass=="numeric"){
      return(as.numeric(cvar))
    }else if(iclass=="factor"){
      return(factor(cvar))
    }else if(iclass=="integer"){
      return(as.integer(cvar))
    }else if(iclass=="character"){
      return(as.character(cvar))
    }else{
      return(var)
    }
  }

  na.time = match.arg(na.time)

  vnames = colnames(data)
  tvar.pos = match(tvar,vnames)
  tcov.pos = match(tcov,vnames)
  cnames = vnames[-c(tvar.pos,tcov.pos)]

  rlen = dim(data)[1]

  data.cov = data[,-c(tvar.pos,tcov.pos)]
  clen  = length(cnames)
  if(clen==1){
    tclass = class(data.cov)
    data.cov = data.frame(data.cov)
    colnames(data.cov) = cnames
    data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
  }

  tte = as.numeric(as.character(data[,tvar.pos[1]]))
  tts = as.numeric(as.character(data[,tvar.pos[2]]))

  tce = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(1,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)
  tcs = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(2,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)

  if(na.time=="remove"){
    na.t = which(is.na(tts))
    na.t = unique(c(na.t,which(is.na(tte))))
    na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tce),"na.action"))))
    na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tcs),"na.action"))))
    if(length(na.t)>0){
      tte = tte[-na.t]
      tts = tts[-na.t]
      tce = tce[-na.t,]
      tcs = tcs[-na.t,]
      data.cov = data.cov[-na.t,]
      if(clen==1){
        tclass = class(data.cov)
        data.cov = data.frame(data.cov)
        colnames(data.cov) = cnames
        data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
      }

      rlen = dim(data.cov)[1]

      tce = matrix(tce,rlen,length(tcov.pos)/2)
      tcs = matrix(tcs,rlen,length(tcov.pos)/2)
    }
  }else{ # censor, take the longest time
    max.t = max(tte,na.rm=T)
    na.t = which(is.na(tts))
    tts[na.t] = 0
    na.t = which(is.na(tte))
    tte[na.t] = max.t
    for(i in 1:dim(tce)[2]){
      max.t = max(tce[,i],na.rm=T)
      na.t = which(is.na(tcs[,i]))
      tcs[na.t,i] = 0
      na.t = which(is.na(tce[,i]))
      tce[na.t,i] = max.t
    }
  }

  ntd = c()
  ncd = c()
  for(i in 1:rlen){
    tmp.cs = c(tcs[i,])
    tmp2.cevent = which(tmp.cs!=0)
    tmp2.ct = c(tce[i,tmp2.cevent])
    tmp2.cct = which(tmp2.ct<tte[i])
    if(length(tmp2.cct)==0){
      tmp.tte = tte[i]
      if(tmp.tte<=0)
        tmp.tte = 0
      ntd = rbind(ntd,c(0,tmp.tte,tts[i],rep(0,length(tmp.cs))))
      ncd = c(ncd,i)
    }else{
      tmp3.cevent = tmp2.cevent[tmp2.cct]
      tmp3.ct = tmp2.ct[tmp2.cct]
      tmp3.cct = which(tmp3.ct>0)
      if(length(tmp3.cct)==0){
        tmpcs = rep(0,length(tmp.cs))
        tmpcs[tmp3.cevent] = 1
        ntd = rbind(ntd,c(0,tte[i],tts[i],tmpcs))
        ncd = c(ncd,i)
      }else{
        tmp3.cct0 = which(tmp3.ct<=0)
        tmpcs = rep(0,length(tmp.cs))
        if(length(tmp3.cct0)>0){
          tmp4.revent = tmp3.cevent[tmp3.cct0]
          tmpcs[tmp4.revent] = 1
        }
        tmp4.cevent = tmp3.cevent[tmp3.cct]
        tmp4.ct = tmp3.ct[tmp3.cct]
        tmp.stime = sort(unique(c(0,tmp4.ct,tte[i])))
        for(j in 1:(length(tmp.stime)-1)){
          if(j<(length(tmp.stime)-1)){
            tmpcs2 = tmpcs
            tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
            tmpcs2[tmp4.cevent[tmp.pos.event]] = 1
            ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],0,tmpcs2))
            ncd = c(ncd,i)
          }else{
            tmpcs2 = tmpcs
            tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
            tmpcs2[tmp4.cevent[tmp.pos.event]] = 1
            ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],tts[i],tmpcs2))
            ncd = c(ncd,i)
          }
        }
      }
    }
  }

  ntd = data.frame(ntd)
  colnames(ntd) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)])
  ndata.cov = data.cov[ncd,]

  dtimecov = data.frame(ntd,ndata.cov)
  colnames(dtimecov) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)],cnames)
  if(verbose){
    print(data)
    print(utils::str(data))
    print(utils::str(dtimecov))
  }

  if(is.null(weights)){
    wt = NULL
  }else{
    wt = weights[ncd]
  }

  list(data=dtimecov,wt=wt)
}

# time-varying cox analysis
tcoxph <- function (formula, na.time=c("remove","censor"),verbose=FALSE,
                     data, weights, subset, na.action, init, control,
                     ties = c("efron", "breslow", "exact"), singular.ok = TRUE,
                     robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties,
                     ...)
{
  time2data <- function(tvar,tcov,data,na.time=c("remove","censor"),verbose=FALSE,weights=NULL){

    cov.class <- function(var,iclass){
      cvar = as.character(var)
      if(iclass=="numeric"){
        return(as.numeric(cvar))
      }else if(iclass=="factor"){
        return(factor(cvar))
      }else if(iclass=="integer"){
        return(as.integer(cvar))
      }else if(iclass=="character"){
        return(as.character(cvar))
      }else{
        return(var)
      }
    }

    na.time = match.arg(na.time)

    vnames = colnames(data)
    tvar.pos = match(tvar,vnames)
    tcov.pos = match(tcov,vnames)
    cnames = vnames[-c(tvar.pos,tcov.pos)]

    rlen = dim(data)[1]

    data.cov = data[,-c(tvar.pos,tcov.pos)]
    clen  = length(cnames)
    if(clen==1){
      tclass = class(data.cov)
      data.cov = data.frame(data.cov)
      colnames(data.cov) = cnames
      data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
    }

    tte = as.numeric(as.character(data[,tvar.pos[1]]))
    tts = as.numeric(as.character(data[,tvar.pos[2]]))

    tce = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(1,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)
    tcs = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(2,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)

    if(na.time=="remove"){
      na.t = which(is.na(tts))
      na.t = unique(c(na.t,which(is.na(tte))))
      na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tce),"na.action"))))
      na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tcs),"na.action"))))
      if(length(na.t)>0){
        tte = tte[-na.t]
        tts = tts[-na.t]
        tce = tce[-na.t,]
        tcs = tcs[-na.t,]
        data.cov = data.cov[-na.t,]
        if(clen==1){
          tclass = class(data.cov)
          data.cov = data.frame(data.cov)
          colnames(data.cov) = cnames
          data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
        }

        rlen = dim(data.cov)[1]

        tce = matrix(tce,rlen,length(tcov.pos)/2)
        tcs = matrix(tcs,rlen,length(tcov.pos)/2)
      }
    }else{ # censor, take the longest time
      max.t = max(tte,na.rm=T)
      na.t = which(is.na(tts))
      tts[na.t] = 0
      na.t = which(is.na(tte))
      tte[na.t] = max.t
      for(i in 1:dim(tce)[2]){
        max.t = max(tce[,i],na.rm=T)
        na.t = which(is.na(tcs[,i]))
        tcs[na.t,i] = 0
        na.t = which(is.na(tce[,i]))
        tce[na.t,i] = max.t
      }
    }

    ntd = c()
    ncd = c()
    for(i in 1:rlen){
      tmp.cs = c(tcs[i,])
      tmp2.cevent = which(tmp.cs!=0)
      tmp2.ct = c(tce[i,tmp2.cevent])
      tmp2.cct = which(tmp2.ct<tte[i])
      if(length(tmp2.cct)==0){
        tmp.tte = tte[i]
        if(tmp.tte<=0)
          tmp.tte = 0
        ntd = rbind(ntd,c(0,tmp.tte,tts[i],rep(0,length(tmp.cs))))
        ncd = c(ncd,i)
      }else{
        tmp3.cevent = tmp2.cevent[tmp2.cct]
        tmp3.ct = tmp2.ct[tmp2.cct]
        tmp3.cct = which(tmp3.ct>0)
        if(length(tmp3.cct)==0){
          tmpcs = rep(0,length(tmp.cs))
          tmpcs[tmp3.cevent] = 1
          ntd = rbind(ntd,c(0,tte[i],tts[i],tmpcs))
          ncd = c(ncd,i)
        }else{
          tmp3.cct0 = which(tmp3.ct<=0)
          tmpcs = rep(0,length(tmp.cs))
          if(length(tmp3.cct0)>0){
            tmp4.revent = tmp3.cevent[tmp3.cct0]
            tmpcs[tmp4.revent] = 1
          }
          tmp4.cevent = tmp3.cevent[tmp3.cct]
          tmp4.ct = tmp3.ct[tmp3.cct]
          tmp.stime = sort(unique(c(0,tmp4.ct,tte[i])))
          for(j in 1:(length(tmp.stime)-1)){
            if(j<(length(tmp.stime)-1)){
              tmpcs2 = tmpcs
              tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
              tmpcs2[tmp4.cevent[tmp.pos.event]] = 1
              ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],0,tmpcs2))
              ncd = c(ncd,i)
            }else{
              tmpcs2 = tmpcs
              tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
              tmpcs2[tmp4.cevent[tmp.pos.event]] = 1
              ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],tts[i],tmpcs2))
              ncd = c(ncd,i)
            }
          }
        }
      }
    }

    ntd = data.frame(ntd)
    colnames(ntd) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)])
    ndata.cov = data.cov[ncd,]

    dtimecov = data.frame(ntd,ndata.cov)
    colnames(dtimecov) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)],cnames)
    if(verbose){
      print(data)
      print(utils::str(data))
      print(utils::str(dtimecov))
    }

    if(is.null(weights)){
      wt = NULL
    }else{
      wt = weights[ncd]
    }

    list(data=dtimecov,wt=wt)
  }
    if (missing(control))
    control <- coxph.control(...)
  if(missing(na.action)){
    na.action = options()$na.action
  }
  if(missing(init)){
    init = NULL
  }
  if(missing(subset)){
    subset = NULL
  }

  tmp.call <- stats::terms(formula)
  tmp.labels <- attr(tmp.call,"term.labels")
  tmp.cpos <- grep("time",tmp.labels)
  if(length(tmp.cpos)==0){
    warning("No time-varying covariate!!!\n")
  }else{
    tmp.mcall <- as.character(tmp.call[[2]])
    if(length(tmp.mcall)>3){
      stop("In this moment, tcoxph does not support the interval-censoring case!!!\n")
    }
    tmp.tte <- tmp.mcall[2]
    tmp.tts <- tmp.mcall[3]
    tmp.main = c(tmp.tte,tmp.tts)

    tmp.covs <- c()
    for(i in tmp.cpos){
      tmp.c1 <- stats::terms(stats::as.formula(paste(tmp.labels[i],"~1",sep="")))
      tmp.c2 <- as.character(tmp.c1[[2]])
      if(length(tmp.c2)!=3){
        stop("Time-varying covariates should be as 'time(time,status)'!!!")
      }
      tmp.covs <- c(tmp.covs,tmp.c2[2:3])
      tmp.labels[i] <- tmp.c2[3]
    }

    if(missing(weights))
      weights <- NULL

    time.cox <- time2data(tvar=tmp.main,tcov=tmp.covs,data=data,na.time=na.time,verbose=verbose,weights=weights)
    formula <- stats::as.formula(paste("Surv(start,end,",tmp.main[2],")~",paste(tmp.labels,collapse="+",sep=""),sep=""))
    weights <- time.cox$wt
    data <- time.cox$data
  }

  pos.time.data <- 1
  assign(".time.data",data, envir = as.environment(pos.time.data))

  ties <- match.arg(ties)
  Call <- match.call()
  Call[[2]] = stats::as.formula(formula)
  Call[[3]] = as.name(".time.data")
  indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                names(Call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  special <- c("strata", "cluster", "tt")
  temp$formula <- if (missing(data))
    stats::terms(formula, special)
  else stats::terms(formula, special, data = data)

  if (is.R())
    m <- eval(temp, parent.frame())
  else m <- eval(temp, sys.parent())
  if (nrow(m) == 0)
    stop("No (non-missing) observations")
  Terms <- stats::terms(m)
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L))
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx ==
                                                                  0L]), domain = NA)
  }
  if (missing(control))
    control <- coxph.control(...)
  Y <- stats::model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type != "right" && type != "counting")
    stop(paste("Cox model doesn't support \"", type, "\" survival data",
               sep = ""))
  weights <- model.weights(m)
  data.n <- nrow(Y)
  strats <- attr(Terms, "specials")$strata
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$terms) > 0)
      Terms2 <- Terms[-stemp$terms]
    else Terms2 <- Terms
    if (length(stemp$vars) == 1)
      strata.keep <- m[[stemp$vars]]
    else strata.keep <- strata(m[, stemp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
  }
  else Terms2 <- Terms
  timetrans <- attr(Terms, "specials")$tt
  if (length(timetrans)) {
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (missing(tt) || is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt))
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function)))
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(strats, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(strats[sorted]) !=
                                           0)))
      }
      if (storage.mode(Y) != "double")
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(strats, -Y[, 2], Y[, 3])
        sort.start <- order(strats, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(strats[sort.end]) !=
                                       0))
      }
      if (storage.mode(Y) != "double")
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount2, Y, as.integer(sort.start -
                                                  1L), as.integer(sort.end - 1L), as.integer(newstrat))
      tindex <- counts$index
    }
    m <- m[tindex, ]
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- "right"
    strats <- factor(rep(1:length(counts$nrisk), counts$nrisk))
    weights <- model.weights(m)
    for (i in 1:ntrans) m[[timetrans$var[i]]] <- (tt[[i]])(m[[timetrans$var[i]]],
                                                           Y[, 1], strats, weights)
  }
  offset <- model.offset(m)
  if (is.null(offset) | all(offset == 0))
    offset <- rep(0, nrow(m))
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) {
    robust <- TRUE
    tempc <- untangle.specials(Terms2, "cluster", 1:10)
    ord <- attr(Terms2, "order")[tempc$terms]
    if (any(ord > 1))
      stop("Cluster can not be used in an interaction")
    cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
    Terms2 <- Terms2[-tempc$terms]
  }
  else {
    if (!missing(robust))
      warning("The robust option is depricated")
    else robust <- FALSE
  }
  attr(Terms2, "intercept") <- 1
  X <- stats::model.matrix(Terms2, m)
  Xatt <- attributes(X)
  if (is.R()) {
    assign <- lapply(attrassign(X, Terms2)[-1], function(x) x -
                       1)
    xlevels <- .getXlevels(Terms2, m)
    contr.save <- attr(X, "contrasts")
  }
  else {
    assign <- lapply(attr(X, "assign")[-1], function(x) x -
                       1)
    xvars <- as.character(attr(Terms2, "variables"))
    xvars <- xvars[-attr(Terms2, "response")]
    if (length(xvars) > 0) {
      xlevels <- lapply(m[xvars], levels)
      xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
      if (length(xlevels) == 0)
        xlevels <- NULL
    }
    else xlevels <- NULL
    contr.save <- attr(X, "contrasts")
  }
  X <- X[, -1, drop = F]
  if (missing(init))
    init <- NULL
  pterms <- sapply(m, inherits, "coxph.penalty")
  if (any(pterms)) {
    pattr <- lapply(m[pterms], attributes)
    pname <- names(pterms)[pterms]
    ord <- attr(Terms, "order")[match(pname, attr(Terms,
                                                  "term.labels"))]
    if (any(ord > 1))
      stop("Penalty terms cannot be in an interaction")
    pcols <- assign[match(pname, names(assign))]
    fit <- coxpenal.fit(X, Y, strats, offset, init = init,
                        control, weights = weights, method = method, row.names(m),
                        pcols, pattr, assign)
  }
  else {
    if (method == "breslow" || method == "efron") {
      if (type == "right")
        fitter <- get("coxph.fit")
      else fitter <- get("agreg.fit")
    }
    else if (method == "exact") {
      if (type == "right")
        fitter <- get("coxexact.fit")
      else fitter <- get("agexact.fit")
    }
    else stop(paste("Unknown method", method))
    fit <- fitter(X, Y, strats, offset, init, control, weights = weights,
                  method = method, row.names(m))
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    if (is.R())
      class(fit) <- "coxph"
    else oldClass(fit) <- "coxph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse = " "))
      if (singular.ok)
        warning(msg)
      else stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    fit$assign <- assign
    if (is.R())
      class(fit) <- fit$method
    else oldClass(fit) <- fit$method[1]
    if (robust) {
      fit$naive.var <- fit$var
      fit$method <- method
      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
      if (length(strats))
        fit2$strata <- strats
      if (length(cluster)) {
        temp <- residuals.coxph(fit2, type = "dfbeta",
                                collapse = cluster, weighted = TRUE)
        if (is.null(init))
          fit2$linear.predictors <- 0 * fit$linear.predictors
        else fit2$linear.predictors <- c(X %*% init)
        temp0 <- residuals.coxph(fit2, type = "score",
                                 collapse = cluster, weighted = TRUE)
      }
      else {
        temp <- residuals.coxph(fit2, type = "dfbeta",
                                weighted = TRUE)
        fit2$linear.predictors <- 0 * fit$linear.predictors
        temp0 <- residuals.coxph(fit2, type = "score",
                                 weighted = TRUE)
      }
      fit$var <- t(temp) %*% temp
      u <- apply(as.matrix(temp0), 2, sum)
      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u,
                                 control$toler.chol)$test
    }
    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (is.null(init))
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
 	      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta],
                                    temp, control$toler.chol)$test
    }
    na.action <- attr(m, "na.action")
    if (length(na.action))
      fit$na.action <- na.action
    if (model) {
      if (length(timetrans)) {
        m[[".surv."]] <- Y
        m[[".strata."]] <- strats
        stop("Time transform + model frame: code incomplete")
      }
      fit$model <- m
    }
    if (x) {
      Xatt$dim <- attr(X, "dim")
      Xatt$dimnames <- attr(X, "dimnames")
      Xatt$assign <- Xatt$assign[-1]
      attributes(X) <- Xatt
      fit$x <- X
      if (length(strats)) {
        if (length(timetrans))
          fit$strata <- strats
        else fit$strata <- strata.keep
      }
    }
    if (y)
      fit$y <- Y
  }
  if (!is.null(weights) && any(weights != 1))
    fit$weights <- weights
  names(fit$means) <- names(fit$coefficients)
  fit$formula <- formula(Terms)
  if (length(xlevels) > 0)
    fit$xlevels <- xlevels
  fit$contrasts <- contr.save
  if (any(offset != 0))
    fit$offset <- offset
  fit$call <- Call
  fit$method <- method
  fit
}
environment(tcoxph) <- environment(coxph)

tcomp.risk <- function (formula, na.time=c("remove","censor"), verbose=FALSE,
	data = sys.parent(), cause, times = NULL,
    	Nit = 50, clusters = NULL, est = NULL, fix.gamma = 0, gamma = 0,
    	n.sim = 0, weighted = 0, model = "fg", detail = 0, interval = 0.01,
    	resample.iid = 1, cens.model = "KM", cens.formula = NULL,
    	time.pow = NULL, time.pow.test = NULL, silent = 1, conv = 1e-06,
    	weights = NULL, max.clust = 1000, n.times = 50, first.time.p = 0.05,
    	estimator = 1, trunc.p = NULL, cens.weights = NULL, admin.cens = NULL,
    	conservative = 1, monotone = 0, step = NULL)
{
  time2data <- function(tvar,tcov,data,na.time=c("remove","censor"),verbose=FALSE,weights=NULL){
    cov.class <- function(var,iclass){
      cvar = as.character(var)
      if(iclass=="numeric"){
        return(as.numeric(cvar))
      }else if(iclass=="factor"){
        return(factor(cvar))
      }else if(iclass=="integer"){
        return(as.integer(cvar))
      }else if(iclass=="character"){
        return(as.character(cvar))
      }else{
        return(var)
      }
    }

    na.time = match.arg(na.time)

    vnames = colnames(data)
    tvar.pos = match(tvar,vnames)
    tcov.pos = match(tcov,vnames)
    cnames = vnames[-c(tvar.pos,tcov.pos)]

    rlen = dim(data)[1]

    data.cov = data[,-c(tvar.pos,tcov.pos)]
    clen  = length(cnames)
    if(clen==1){
      tclass = class(data.cov)
      data.cov = data.frame(data.cov)
      colnames(data.cov) = cnames
      data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
    }

    tte = as.numeric(as.character(data[,tvar.pos[1]]))
    tts = as.numeric(as.character(data[,tvar.pos[2]]))

    tce = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(1,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)
    tcs = matrix(as.numeric(as.character(c(unlist(data[,tcov.pos[seq(2,length(tcov.pos),by=2)]])))),rlen,length(tcov.pos)/2)

    if(na.time=="remove"){
      na.t = which(is.na(tts))
      na.t = unique(c(na.t,which(is.na(tte))))
      na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tce),"na.action"))))
      na.t = unique(c(na.t,as.numeric(attr(stats::na.omit(tcs),"na.action"))))
      if(length(na.t)>0){
        tte = tte[-na.t]
        tts = tts[-na.t]
        tce = tce[-na.t,]
        tcs = tcs[-na.t,]
        data.cov = data.cov[-na.t,]
        if(clen==1){
          tclass = class(data.cov)
          data.cov = data.frame(data.cov)
          colnames(data.cov) = cnames
          data.cov[,cnames] = cov.class(data.cov[,cnames],tclass)
        }

        rlen = dim(data.cov)[1]

        tce = matrix(tce,rlen,length(tcov.pos)/2)
        tcs = matrix(tcs,rlen,length(tcov.pos)/2)
      }
    }else{ # censor, take the longest time
      max.t = max(tte,na.rm=T)
      na.t = which(is.na(tts))
      tts[na.t] = 0
      na.t = which(is.na(tte))
      tte[na.t] = max.t
      for(i in 1:dim(tce)[2]){
        max.t = max(tce[,i],na.rm=T)
        na.t = which(is.na(tcs[,i]))
        tcs[na.t,i] = 0
        na.t = which(is.na(tce[,i]))
        tce[na.t,i] = max.t
      }
    }

    ntd = c()
    ncd = c()
    for(i in 1:rlen){
      tmp.cs = c(tcs[i,])
      tmp2.cevent = which(tmp.cs!=0)
      tmp2.ct = c(tce[i,tmp2.cevent])
      tmp2.cct = which(tmp2.ct<tte[i])
      if(length(tmp2.cct)==0){
        tmp.tte = tte[i]
        if(tmp.tte<=0)
          tmp.tte = 0
        ntd = rbind(ntd,c(0,tmp.tte,tts[i],rep(0,length(tmp.cs))))
        ncd = c(ncd,i)
      }else{
        tmp3.cevent = tmp2.cevent[tmp2.cct]
        tmp3.ct = tmp2.ct[tmp2.cct]
        tmp3.cct = which(tmp3.ct>0)
        if(length(tmp3.cct)==0){
          tmpcs = rep(0,length(tmp.cs))
          tmpcs[tmp3.cevent] = 1
          ntd = rbind(ntd,c(0,tte[i],tts[i],tmpcs))
          ncd = c(ncd,i)
        }else{
          tmp3.cct0 = which(tmp3.ct<=0)
          tmpcs = rep(0,length(tmp.cs))
          if(length(tmp3.cct0)>0){
            tmp4.revent = tmp3.cevent[tmp3.cct0]
            tmpcs[tmp4.revent] = 1
          }
          tmp4.cevent = tmp3.cevent[tmp3.cct]
          tmp4.ct = tmp3.ct[tmp3.cct]
          tmp.stime = sort(unique(c(0,tmp4.ct,tte[i])))
          for(j in 1:(length(tmp.stime)-1)){
            if(j<(length(tmp.stime)-1)){
              tmpcs2 = tmpcs
              tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
              tmpcs2[tmp4.cevent[tmp.pos.event]] = 1 ######
              ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],0,tmpcs2))
              ncd = c(ncd,i)
            }else{
              tmpcs2 = tmpcs
              tmp.pos.event = which(tmp4.ct<tmp.stime[j+1])
              tmpcs2[tmp4.cevent[tmp.pos.event]] = 1 ######
              ntd = rbind(ntd,c(tmp.stime[j],tmp.stime[j+1],tts[i],tmpcs2))
              ncd = c(ncd,i)
            }
          }
        }
      }
    }

    ntd = data.frame(ntd)
    colnames(ntd) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)])
    ndata.cov = data.cov[ncd,]

    dtimecov = data.frame(ntd,ndata.cov)
    colnames(dtimecov) = c("start","end",tvar[2],tcov[seq(2,length(tcov),by=2)],cnames)
    if(verbose){
      print(data)
      print(utils::str(data))
      print(utils::str(dtimecov))
    }

    if(is.null(weights)){
      wt = NULL
    }else{
      wt = weights[ncd]
    }

    list(data=dtimecov,wt=wt)
  }

  if (!missing(cause)) {
        if (length(cause) != 1)
            stop("Argument cause has new meaning since \n   timereg version 1.8.4., it now specifies the cause of interest, see help(comp.risk) for details.")
    }

  	sktmp.call <- stats::terms(formula)
  	sktmp.labels <- attr(sktmp.call,"term.labels")
  	sktmp.cpos <- grep("time",sktmp.labels)
  	if(length(sktmp.cpos)==0){
    		warning("No time-varying covariate!!!\n")
  	}else{
    		sktmp.mcall <- as.character(sktmp.call[[2]])
    		if(length(sktmp.mcall)>3){
      		stop("In this moment, tcomp.risk does not support the interval-censoring case!!!\n")
    		}
    		sktmp.tte <- sktmp.mcall[2]
    		sktmp.tts <- sktmp.mcall[3]
    		sktmp.main = c(sktmp.tte,sktmp.tts)

    		sktmp.covs <- c()
    		for(i in sktmp.cpos){
      		sktmp.c1 <- stats::terms(stats::as.formula(paste(sktmp.labels[i],"~1",sep="")))
      		sktmp.c2 <- as.character(sktmp.c1[[2]])
      		if(length(sktmp.c2)!=3){
        			stop("Time-varying covariates should be as 'time(time,status)'!!!")
      		}
      		sktmp.covs <- c(sktmp.covs,sktmp.c2[2:3])
      		sktmp.labels[i] <- sktmp.c2[3]
    		}

    		time.comp.risk <- time2data(tvar=sktmp.main,tcov=sktmp.covs,data=data,na.time=na.time,verbose=verbose,weights=NULL)
    		formula <- stats::as.formula(paste("Event(start,end,",sktmp.main[2],")~",paste(sktmp.labels,collapse="+",sep=""),sep=""))
    		data <- time.comp.risk$data
  	}

  	pos.time.data <- 1
  	assign(".time.data",data, envir = as.environment(pos.time.data))

    trans <- switch(model, additive = 1, prop = 2, logistic = 3,
        rcif = 4, rcif2 = 5, fg = 6, logistic2 = 7)
    line <- 0
    cause.call <- causeS <- cause
    m <- match.call(expand.dots = FALSE)
	m[[2]] <- stats::as.formula(formula)
	m[[3]] <- as.name(".time.data")
    m$gamma <- m$times <- m$n.times <- m$cause <- m$Nit <- m$weighted <- m$n.sim <- m$model <- m$detail <- m$cens.model <- m$time.pow <- m$silent <- m$step <- m$cens.formula <- m$interval <- m$clusters <- m$resample.iid <- m$monotone <- m$time.pow.test <- m$conv <- m$weights <- m$max.clust <- m$first.time.p <- m$trunc.p <- m$cens.weights <- m$admin.cens <- m$fix.gamma <- m$est <- m$conservative <- m$estimator <- NULL
    if ((trans == 2 || trans == 3 || trans == 7) && is.null(step))
        step <- 0.5
    if (is.null(step))
        step <- 1
    special <- c("const", "cluster")

    if (missing(data)) {
        Terms <- stats::terms(formula, special)
    }
    else {
        Terms <- stats::terms(formula, special, data = data)
    }
    m$formula <- Terms
    if (substr(as.character(m$formula)[2], 1, 4) == "Hist") {
        stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as \n       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
    }
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (NROW(m) == 0)
        stop("No (non-missing) observations")
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    event.history <- stats::model.extract(m, "response")
    if (class(event.history) != "Event") {
        stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as \n       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
    }
    model.type <- "competing.risks"
    cens.code <- attr(event.history, "cens.code")
    if (ncol(event.history) == 2) {
        time2 <- eventtime <- event.history[, 1]
        status <- delta <- event.history[, 2]
        entrytime <- rep(0, length(time2))
        left <- 0
    }
    else {
        time2 <- eventtime <- event.history[, 2]
        status <- delta <- event.history[, 3]
        entrytime <- event.history[, 1]
        left <- 1
        if (max(entrytime) == 0)
            left <- 0
    }
    event <- (status == cause)
    if (sum(event) == 0)
        stop("No events of interest in data\n")
    if (n.sim == 0)
        sim <- 0
    else sim <- 1
    antsim <- n.sim
    des <- read.design(m, Terms)
    X <- des$X
    Z <- des$Z
    npar <- des$npar
    px <- des$px
    pz <- des$pz
    covnamesX <- des$covnamesX
    covnamesZ <- des$covnamesZ
    if (nrow(X) != nrow(data))
        stop("Missing values in design matrix not allowed\n")
    if (is.diag(t(X) %*% X) == TRUE)
        stratum <- 1
    else stratum <- 0
    if (is.null(clusters)) {
        clusters <- des$clusters
    }
    if (is.null(clusters)) {
        cluster.call <- clusters
        clusters <- 0:(nrow(X) - 1)
        antclust <- nrow(X)
    }
    else {
        cluster.call <- clusters
        antclust <- length(unique(clusters))
        clusters <- as.integer(factor(clusters, labels = 1:antclust)) -
            1
    }
    coarse.clust <- FALSE
    if ((!is.null(max.clust)))
        if (max.clust < antclust) {
            coarse.clust <- TRUE
            qq <- unique(stats::quantile(clusters, probs = seq(0, 1,
                by = 1/max.clust)))
            qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)
            clusters <- as.integer(qqc) - 1
            max.clusters <- length(unique(clusters))
            antclust <- max.clust
        }
    pxz <- px + pz
    if (is.null(times)) {
        timesc <- sort(unique(eventtime[event == 1]))
        if (!is.null(n.times)) {
            if (length(timesc) > n.times)
                times <- stats::quantile(timesc, prob = seq(first.time.p,
                  1, length = n.times))
            else times <- timesc
        }
        else {
            times <- timesc
            times <- times[times > stats::quantile(timesc, prob = first.time.p)]
        }
    }
    else times <- sort(times)
    n <- nrow(X)
    ntimes <- length(times)
    if (npar == TRUE) {
        Z <- matrix(0, n, 1)
        pg <- 1
        fixed <- 0
    }
    else {
        fixed <- 1
        pg <- pz
    }
    if (is.null(weights) == TRUE)
        weights <- rep(1, n)
    if (!is.null(admin.cens))
        estimator <- 3
    Gcxe <- 1
    ordertime <- order(eventtime)
    if (estimator == 1 || estimator == 2) {
        if (is.null(cens.weights)) {
            if (cens.model == "KM") {
                if (left == 1)
                  ud.cens <- survfit(Surv(entrytime, eventtime,
                    delta == cens.code) ~ +1)
                else ud.cens <- survfit(Surv(eventtime, delta ==
                  cens.code) ~ +1)
                Gfit <- cbind(ud.cens$time, ud.cens$surv)
                Gfit <- rbind(c(0, 1), Gfit)
                Gcx <- Cpred(Gfit, eventtime, strict = TRUE)[,
                  2]
                Gcxe <- Cpred(Gfit, entrytime, strict = TRUE)[,
                  2]
                Gcxe[Gcxe == 0] <- 1
                if (!is.null(trunc.p))
                  Gcx <- Gcx/Gcxe
                Gctimes <- Cpred(Gfit, times, strict = TRUE)[,
                  2]
            }
            else if (cens.model == "stratKM") {
                XZ <- stats::model.matrix(cens.formula, data = data)
                strata <- as.factor(XZ)
                Gcx <- pred.stratKM(data, time = eventtime, cause = delta,
                  strata = strata)
                if (!is.null(trunc.p))
                  Gcx <- Gcx/Gcxe
                Gctimes <- Cpred(Gfit, times)[, 2]
            }
            else if (cens.model == "cox") {
                if (!is.null(cens.formula)) {
                  XZ <- stats::model.matrix(cens.formula, data = data)
                  if (sum(XZ[, 1]) == nrow(XZ))
                    XZ <- as.matrix(XZ[, -1])
                }
                else {
                  if (npar == TRUE)
                    XZ <- X[, -1]
                  else XZ <- cbind(X, Z)[, -1]
                }
                if (left == 1)
                  ud.cens <- coxph(Surv(entrytime, eventtime,
                    delta == cens.code) ~ XZ)
                else ud.cens <- coxph(Surv(eventtime, delta ==
                  cens.code) ~ XZ)
                baseout <- basehaz(ud.cens, centered = FALSE)
                baseout <- cbind(baseout$time, baseout$hazard)
                Gcx <- Cpred(baseout, eventtime, strict = TRUE)[,
                  2]
                Gcxe <- Cpred(baseout, entrytime, strict = TRUE)[,
                  2]
                Gcxe[Gcxe == 0] <- 1
                RR <- exp(as.matrix(XZ) %*% stats::coef(ud.cens))
                Gcx <- exp(-Gcx * RR)
                Gcxe <- exp(-Gcxe * RR)
                Gfit <- rbind(c(0, 1), cbind(eventtime, Gcx))
                if (!is.null(trunc.p))
                  Gcx <- Gcx/Gcxe
                Gctimes <- Cpred(Gfit, times, strict = TRUE)[,
                  2]
            }
            else if (cens.model == "aalen") {
                if (!is.null(cens.formula)) {
                  XZ <- stats::model.matrix(cens.formula, data = data)
                }
                else {
                  if (npar == TRUE)
                    XZ <- X
                  else XZ <- cbind(X, Z)
                }
                if (left == 1)
                  ud.cens <- aalen(Surv(entrytime, eventtime,
                    delta == cens.code) ~ -1 + XZ + cluster(clusters),
                    n.sim = 0, residuals = 0, robust = 0, silent = 1)
                else ud.cens <- aalen(Surv(eventtime, delta ==
                  cens.code) ~ -1 + XZ + cluster(clusters), n.sim = 0,
                  residuals = 0, robust = 0, silent = 1)
                Gcx <- Cpred(ud.cens$cum, eventtime, strict = TRUE)[,
                  -1]
                Gcx <- exp(-apply(Gcx * XZ, 1, sum))
                Gcx[Gcx > 1] <- 1
                Gcx[Gcx < 0] <- 1
                Gcxe <- Cpred(ud.cens$cum, entrytime, strict = TRUE)[,
                  2]
                Gcxe[Gcxe == 0] <- 1
                if (!is.null(trunc.p))
                  Gcx <- Gcx/Gcxe
                Gfit <- rbind(c(0, 1), cbind(eventtime, Gcx))
                Gctimes <- Cpred(Gfit, times, strict = TRUE)[,
                  2]
            }
            else stop("Unknown censoring model")
            cens.weights <- Gcx
            if ((min(Gcx[event == 1]) < 1e-05) && (silent ==
                0)) {
                cat("Censoring dist. approx zero for some points, summary cens:\n")
                print(summary(Gcx))
            }
        }
        else {
            if (length(cens.weights) != n)
                stop("censoring weights must have length equal to nrow in data\n")
            Gcx <- cens.weights
            ord2 <- order(time2)
            Gctimes <- Cpred(cbind(time2[ord2], weights[ord2]),
                times)
        }
    }
    else {
        if (length(admin.cens) != n)
            stop("censoring weights must have length equal to nrow in data\n")
        Gcx <- admin.cens
        Gctimes <- rep(1, length(times))
    }
    if (left == 1 & is.null(trunc.p) & is.null(cens.weights)) {
        stop("For left-truncated data call prep.comp.risk\n call with weights and cens.weights\n")
        n = length(time2)
        prec.factor <- 100
        prec <- .Machine$double.eps * prec.factor
        surv.trunc <- survfit(Surv(-time2, -entrytime + prec,
            rep(1, n)) ~ 1)
        trunc.dist <- summary(surv.trunc)
        trunc.dist$time <- rev(-trunc.dist$time)
        trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
        Lfit <- Cpred(cbind(trunc.dist$time, trunc.dist$surv),
            time2)
        Lw <- Lfit[, 2]
        weights <- 1/((Lw) * Gcx)
        weights[delta == cens.code] <- 0
        Gcx <- rep(1, n)
    }
    if (is.null(trunc.p))
        trunc.p <- rep(1, n)
    if (length(trunc.p) != n)
        stop("truncation weights must have same length as data\n")
    if (resample.iid == 1) {
        biid <- double(ntimes * antclust * px)
        gamiid <- double(antclust * pg)
    }
    else {
        gamiid <- biid <- NULL
    }
    ps <- px
    betaS <- rep(0, ps)
    if (is.null(est)) {
        est <- matrix(0 + 0.1, ntimes, px + 1)
        est[, 1] <- times
    }
    else {
        est <- as.matrix(est)
    }
    if (nrow(est) != length(times))
        est <- Cpred(est, times)
    hess <- matrix(0, ps, ps)
    var <- score <- matrix(0, ntimes, ps + 1)
    if (sum(gamma) == 0)
        gamma <- rep(0, pg)
    gamma2 <- rep(0, ps)
    test <- matrix(0, antsim, 3 * ps)
    testOBS <- rep(0, 3 * ps)
    unifCI <- c()
    testval <- c()
    rani <- -round(stats::runif(1) * 10000)
    Ut <- matrix(0, ntimes, ps + 1)
    simUt <- matrix(0, ntimes, 50 * ps)
    var.gamma <- matrix(0, pg, pg)
    pred.covs.sem <- 0
    if (is.null(time.pow) == TRUE & model == "prop")
        time.pow <- rep(0, pg)
    if (is.null(time.pow) == TRUE & model == "fg")
        time.pow <- rep(0, pg)
    if (is.null(time.pow) == TRUE & model == "additive")
        time.pow <- rep(1, pg)
    if (is.null(time.pow) == TRUE & model == "rcif")
        time.pow <- rep(0, pg)
    if (is.null(time.pow) == TRUE & model == "rcif2")
        time.pow <- rep(0, pg)
    if (is.null(time.pow) == TRUE & model == "logistic")
        time.pow <- rep(0, pg)
    if (is.null(time.pow) == TRUE & model == "logistic2")
        time.pow <- rep(0, pg)
    if (length(time.pow) != pg)
        time.pow <- rep(time.pow[1], pg)
    if (is.null(time.pow.test) == TRUE & model == "prop")
        time.pow.test <- rep(0, px)
    if (is.null(time.pow.test) == TRUE & model == "fg")
        time.pow.test <- rep(0, px)
    if (is.null(time.pow.test) == TRUE & model == "additive")
        time.pow.test <- rep(1, px)
    if (is.null(time.pow.test) == TRUE & model == "rcif")
        time.pow.test <- rep(0, px)
    if (is.null(time.pow.test) == TRUE & model == "rcif2")
        time.pow.test <- rep(0, px)
    if (is.null(time.pow.test) == TRUE & model == "logistic")
        time.pow.test <- rep(0, px)
    if (is.null(time.pow.test) == TRUE & model == "logistic2")
        time.pow.test <- rep(0, px)
    if (length(time.pow.test) != px)
        time.pow.test <- rep(time.pow.test[1], px)
    if (ntimes > 1)
        silent <- c(silent, rep(0, ntimes - 1))
    ssf <- step
    out <- .C("itfit", as.double(times), as.integer(ntimes),
        as.double(eventtime), as.integer(cens.code), as.integer(status),
        as.double(Gcx), as.double(X), as.integer(n), as.integer(px),
        as.integer(Nit), as.double(betaS), as.double(score),
        as.double(hess), as.double(est), as.double(var), as.integer(sim),
        as.integer(antsim), as.integer(rani), as.double(test),
        as.double(testOBS), as.double(Ut), as.double(simUt),
        as.integer(weighted), as.double(gamma), as.double(var.gamma),
        as.integer(fixed), as.double(Z), as.integer(pg), as.integer(trans),
        as.double(gamma2), as.integer(cause), as.integer(line),
        as.integer(detail), as.double(biid), as.double(gamiid),
        as.integer(resample.iid), as.double(time.pow), as.integer(clusters),
        as.integer(antclust), as.double(time.pow.test), as.integer(silent),
        as.double(conv), as.double(weights), as.double(entrytime),
        as.double(trunc.p), as.integer(estimator), as.integer(fix.gamma),
        as.integer(stratum), as.integer(ordertime - 1), as.integer(conservative),
        as.double(ssf), as.double(Gctimes), as.double(rep(0,
            pg)), as.double(matrix(0, pg, pg)), as.integer(monotone),
        PACKAGE = "timereg")
    ssf <- out[[51]]
    gamma <- matrix(out[[24]], pg, 1)
    var.gamma <- matrix(out[[25]], pg, pg)
    Dscore.gamma <- matrix(out[[54]], pg, pg)
    gamma2 <- matrix(out[[30]], ps, 1)
    rownames(gamma2) <- covnamesX
    conv <- list(convp = out[[41]], convd = out[[42]])
    if (fixed == 0)
        gamma <- NULL
    if (resample.iid == 1) {
        biid <- matrix(out[[34]], ntimes, antclust * px)
        if (fixed == 1)
            gamiid <- matrix(out[[35]], antclust, pg)
        else gamiid <- NULL
        B.iid <- list()
        for (i in (0:(antclust - 1)) * px) {
            B.iid[[i/px + 1]] <- matrix(biid[, i + (1:px)], ncol = px)
            colnames(B.iid[[i/px + 1]]) <- covnamesX
        }
        if (fixed == 1)
            colnames(gamiid) <- covnamesZ
    }
    else B.iid <- gamiid <- NULL
    if (sim == 1) {
        simUt <- matrix(out[[22]], ntimes, 50 * ps)
        UIt <- list()
        for (i in (0:49) * ps) UIt[[i/ps + 1]] <- as.matrix(simUt[,
            i + (1:ps)])
        Ut <- matrix(out[[21]], ntimes, ps + 1)
        test <- matrix(out[[19]], antsim, 3 * ps)
        testOBS <- out[[20]]
        supUtOBS <- apply(abs(as.matrix(Ut[, -1])), 2, max)
        p <- ps
        for (i in 1:(3 * p)) testval <- c(testval, pval(test[,
            i], testOBS[i]))
        for (i in 1:p) unifCI <- as.vector(c(unifCI, percen(test[,
            i], 0.95)))
        pval.testBeq0 <- as.vector(testval[1:p])
        pval.testBeqC <- as.vector(testval[(p + 1):(2 * p)])
        pval.testBeqC.is <- as.vector(testval[(2 * p + 1):(3 *
            p)])
        obs.testBeq0 <- as.vector(testOBS[1:p])
        obs.testBeqC <- as.vector(testOBS[(p + 1):(2 * p)])
        obs.testBeqC.is <- as.vector(testOBS[(2 * p + 1):(3 *
            p)])
        sim.testBeq0 <- as.matrix(test[, 1:p])
        sim.testBeqC <- as.matrix(test[, (p + 1):(2 * p)])
        sim.testBeqC.is <- as.matrix(test[, (2 * p + 1):(3 *
            p)])
    }
    else {
        test <- unifCI <- Ut <- UIt <- pval.testBeq0 <- pval.testBeqC <- obs.testBeq0 <- obs.testBeqC <- sim.testBeq0 <- sim.testBeqC <- sim.testBeqC.is <- pval.testBeqC.is <- obs.testBeqC.is <- NULL
    }
    est <- matrix(out[[14]], ntimes, ps + 1)
    est[conv$convp > 0, -1] <- NA
    score <- matrix(out[[12]], ntimes, ps + 1)
    gamscore <- matrix(out[[53]], pg, 1)
    scores <- list(score = score, gamscore = gamscore)
    var <- matrix(out[[15]], ntimes, ps + 1)
    var[conv$convp > 0, -1] <- NA
    colnames(var) <- colnames(est) <- c("time", covnamesX)
    if (sim >= 1) {
        colnames(Ut) <- c("time", covnamesX)
        names(unifCI) <- names(pval.testBeq0) <- names(pval.testBeqC) <- names(pval.testBeqC.is) <- names(obs.testBeq0) <- names(obs.testBeqC) <- names(obs.testBeqC.is) <- colnames(sim.testBeq0) <- colnames(sim.testBeqC) <- colnames(sim.testBeqC.is) <- covnamesX
    }
    if (fixed == 1) {
        rownames(gamma) <- c(covnamesZ)
        colnames(var.gamma) <- rownames(var.gamma) <- c(covnamesZ)
    }
    colnames(score) <- c("time", covnamesX)
    if (is.na(sum(score)) == TRUE)
        score <- NA
    else if (sum(score[, -1]) < 1e-05)
        score <- sum(score[, -1])
    ud <- list(cum = est, var.cum = var, gamma = gamma, score = score,
        gamma2 = gamma2, var.gamma = var.gamma, robvar.gamma = var.gamma,
        pval.testBeq0 = pval.testBeq0, pval.testBeqC = pval.testBeqC,
        obs.testBeq0 = obs.testBeq0, obs.testBeqC.is = obs.testBeqC.is,
        obs.testBeqC = obs.testBeqC, pval.testBeqC.is = pval.testBeqC.is,
        conf.band = unifCI, B.iid = B.iid, gamma.iid = gamiid,
        ss = ssf, test.procBeqC = Ut, sim.test.procBeqC = UIt,
        conv = conv, weights = weights, cens.weights = cens.weights,
        scores = scores, Dscore.gamma = Dscore.gamma, step = step)
    ud$call <- call
    ud$model <- model
    ud$n <- n
    ud$clusters <- clusters
    ud$formula <- formula
    ud$response <- event.history
    ud$cause <- status
    class(ud) <- "comprisk"
    attr(ud, "Call") <- call
    attr(ud, "Formula") <- formula
    attr(ud, "time.pow") <- time.pow
    attr(ud, "causeS") <- causeS
    attr(ud, "cause") <- status
    attr(ud, "cluster.call") <- cluster.call
    attr(ud, "coarse.clust") <- coarse.clust
    attr(ud, "max.clust") <- max.clust
    attr(ud, "clusters") <- clusters
    attr(ud, "cens.code") <- cens.code
    attr(ud, "times") <- times
    return(ud)
}
environment(tcomp.risk) = environment(comp.risk)

