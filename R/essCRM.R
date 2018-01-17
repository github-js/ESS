#' @export
#' @importFrom dfcrm titesim

essCRM <- function(PI,prior,betaSD,target,m,nsim,obswin=30,rate=2,accrual="poisson"){

  library(dfcrm)
  mRange <- 0:m
  numMC <- nsim
  getDiff <- function(d,w,y){
    denom <- d*w - 1
    term1 <- d*(log(d)^2)*w*(y-d*w)
    term2 <- d*(log(d)^2)*w
    term3 <- log(d)*(y-d*w)
    result <- term1/(denom^2) + term2/denom - term3/denom
    result
  }

  deltaBar <- rep(NA,length(mRange))

  for(q in 1:length(mRange)){
    m <- mRange[q]
    sampMC <- rep(NA,numMC)
    Dq <- 0
    if(m>0){
      for(k in 1:numMC){
        set.seed(k)
        simData <- titesim(PI, prior,
                           target, n=max(mRange),
                           x0=1, nsim=1,
                           obswin=30, rate=2,
                           accrual="poisson",
                           scale=betaSD, seed=k)

        get <- sort(sample(1:max(mRange),m))
        d <- prior[simData$level][get]
        y <- simData$tox[get]
        Tup <- 30
        u <- simData$arrival[get]
        u[u=='Inf' | u > Tup] <- Tup
        mydim <- length(simData$prior)

        w <- u/Tup
        Dq <- 0
        for(i in 1:m){
          Dq <- Dq + getDiff(d=d[i],w=w[i],y=y[i])
        }
        sampMC[k] <- Dq
      }
      deltaBar[q] <- 1/((betaSD)^2) + mean(sampMC)
    }

    if(m==0){
      deltaBar[q] <- 1/((betaSD)^2)
    }
  }

  min.n.index <- which.min(abs(deltaBar))
  min.n <- mRange[which.min(abs(deltaBar))]
  min.v <- deltaBar[which.min(abs(deltaBar))]
  interpolated <- approx(mRange, deltaBar, method = "linear")
  ESS <- interpolated$x[which.min(abs(interpolated$y))]
  ESS
}
