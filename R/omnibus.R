omnibus.test <- function(p, method=c("z", "p", "log.p"), N.sim=10000, approximate=FALSE) {
  # make sure that x is a matrix
  if(is.vector(p))
    p <- matrix(p, ncol=1)
  if(!is.matrix(p))
    stop("p has to be either a vector or matrix.")

  if(approximate && method != "p")
    stop("Beta approximation currently only available if method = 'p'")

  # determine which transformation of input p-values should be used
  trafo <- NULL
  if(is.character(method)) {
    # match with pre-defined method
    method <- match.arg(method)
    # set transform method
    trafo <- switch (method,
                    z = function(x) qnorm(x, lower.tail=FALSE),
                    p = function(x) 1-x,
                    log.p = function(x) -log(x)
    )
  } else if(is.function(method)) {
    trafo <- method
  } else {
    stop("The method parameter has not been specified properly. Has to be either a character string naming a predefined method or a function.")
  }

  # transform input p-values
  stat <- trafo(p)

  N.hyp <- nrow(p)
  # generate statistic under H0
  stat.H0 <- matrix(trafo(runif(N.sim*N.hyp)), nrow=N.hyp)

  # compute cumulative sums (sort statistic decreasing, compute cumulative sums and divide by number of combined hypotheses)
  csStat <- colCumsums(apply(stat, 2, sort, decreasing=TRUE)) / 1:N.hyp
  csStat.H0 <- colCumsums(apply(stat.H0, 2, sort, decreasing=TRUE)) / 1:N.hyp

  # transform sums using global H0
  Fh.params <- matrix(NA, nrow=N.hyp, ncol=2)
  Fh.H0 <- matrix(NA, nrow=N.hyp, ncol=N.sim)
  Fh <- matrix(NA, nrow=N.hyp, ncol=ncol(p))
  for (h in 1:N.hyp)
  {
    if(approximate) {
      # compute beta approximation under H_0 to obtain F_h
      Fh.params[h,]  <-  fitdist(csStat.H0[h,], "beta","mme")$estimate
      Fh.H0[h,] <- pbeta(csStat.H0[h,],shape1= Fh.params[h,1], shape2= Fh.params[h,2])
      Fh[h,] <- pbeta(csStat[h,],shape1= Fh.params[h,1], shape2= Fh.params[h,2])
    } else {
      # transform sums using simulations under H_0
      Fh.H0[h,] <- rank(csStat.H0[h,])/ncol(csStat.H0)
      Fh[h,] <- ecdf(csStat.H0[h,])(csStat[h,])
    }
  }

  # compute p-values based on max statistic
  if(approximate) {
    # using beta approximation
    G.params <- fitdist(colMaxs(Fh.H0), "beta","qme", probs=c(0.9, 0.999))$estimate #fitdist(colMaxs(Fh.H0), "beta","mme")$estimate
    return(1-pbeta(colMaxs(Fh), shape1=G.params[1], shape2=G.params[2]))

  } else {
    # using simulations under H_0
    return(1-ecdf(colMaxs(Fh.H0))(colMaxs(Fh)))
  }

}

hc.test <- function(p, tuning=c("half", "halfmin", "all"), N.sim=10000) {
  # make sure that p is a matrix
  if(is.vector(p))
    p <- matrix(p, ncol=1)
  if(!is.matrix(p))
    stop("p has to be either a vector or matrix.")

  tuning <- match.arg(tuning)
  N.hyp <- nrow(p)

  # make sure that there are enough hypotheses
  if(N.hyp < 2 || (tuning == "halfmin" && N.hyp < 3))
    stop("At least 2 hypotheses need to be tested. If tuning equals 'halfmin' then at least 3 hypotheses are required.")

  # generate p-values under H0
  p.H0 <- matrix(runif(N.hyp*N.sim), nrow=N.hyp)

  # sort p-values
  p <- apply(p, 2, sort)
  p.H0 <- apply(p.H0, 2, sort)

  # compute test statistic
  emp.df <- (1:N.hyp)/N.hyp
  stat  <- sqrt(N.hyp)*(emp.df-p)/sqrt(p*(1-p))
  stat.H0  <- sqrt(N.hyp)*(emp.df-p.H0)/sqrt(p.H0*(1-p.H0))
  stat <- switch(tuning,
                 half=colMaxs(stat[1:ceiling(N.hyp/2),]),
                 halfmin=colMaxs(stat[2:ceiling(N.hyp/2),]),
                 all=colMaxs(stat))
  stat.H0 <- switch(tuning,
                    half=colMaxs(stat.H0[1:ceiling(N.hyp/2),]),
                    halfmin=colMaxs(stat.H0[2:ceiling(N.hyp/2),]),
                    all=colMaxs(stat.H0))

  # compute p-values based on max statistic
  return(1-ecdf(stat.H0)(stat))
}








