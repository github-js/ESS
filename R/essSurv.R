#' @export
#' @importFrom MCMCpack rinvgamma

essSurv <- function(shapeParam,scaleParam,m,nsim){
  library(MCMCpack)
  mRange <- 1:m
  numMC <- nsim
  ## Generate prior from inverse gamma distribution
  myMu <- rinvgamma(numMC, shape=shapeParam, scale = scaleParam)

  getDp <- function(alpha,beta){
    myMu <- beta/(alpha-1)
    myDp <- -(alpha+1)/(myMu^2) + 2*beta/(myMu^3)
    myDp
  }

  Dp <- getDp(alpha = shapeParam, beta = scaleParam)
  Dqm <- rep(NA,length(mRange))


  for(i in 1:length(mRange)){
    m <- mRange[i]
    myDq <- rep(NA,numMC)
    for(q in 1:numMC){
      genData <- function(m,rateParam,censtime){
        lifetime <- rexp(m, rate = rateParam)
        t0 <- pmin(lifetime, censtime)
        y <- as.numeric(censtime > lifetime)
        data <- cbind(y,t0)
        data
      }

      myData <- genData(m,rateParam=1/myMu[q],censtime=3)

      getDq <- function(myMu,y,t0){
        myDq <- (1+sum(y))*(-1/(myMu^2)) + (2/(myMu^3))*sum(t0)
        myDq
      }

      myDq[q] <- getDq(myMu=myMu[q], y = myData[,1], t0 = myData[,2])
    }

    Dqm[i] <- mean(myDq)
    #print(Dp-Dqm[i])
    rm(myDq)
  }

  deltaBar <- Dp - Dqm
  min.n.index <- which.min(abs(deltaBar))
  min.n <- mRange[which.min(abs(deltaBar))]
  min.v <- deltaBar[which.min(abs(deltaBar))]

  interpolated <- approx(mRange, deltaBar, method = "linear")

  ESS <- interpolated$x[which.min(abs(interpolated$y))]

  ESS
}
