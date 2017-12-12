

##################################################################################
##
## Main R function written by Satoshi Morita
## Downloaded 12/12/17 via https://biostatistics.mdanderson.org/softwaredownload/
##
##################################################################################

#***************************************************************************
#        R code for determining the Effective Sample Size (ESS)
#        of a logistic or linear regression model
#                 Version 1.0,  11Aug2009
#***************************************************************************


#For the example input data shown here:
#This normal linear regression example
#            should calculate the ESS  of 2 for the whole vector (theta0,theta_1,theta_2,theta_3,tau),
#                             the ESS1 of 2 for the subvector 1 (theta0,theta_1,theta_2,theta_3), and
#                             the ESS2 of 2 for the subvector 2 (tau).



## Notes on a regression model ##
## Let theta_j for j=0,...,d be a parameter for a regression model.
## Let theta_0 be a parameter for an intercept.
## Let theta_1,..., theta_d be parameters for regression coefficients of covariates X_1,..., X_d.
## Thus, a linear term of a selected regression model is given as theta_0 + X_1*theta_1 + ... + X_d*theta_d.


##### Please input the follwoing information at Steps 1 to 5 #####

## Step 1. Specify a regression model by imputing 1 or 2:
## 1 for a linear regression model, --> SEE Note (2) below
## 2 for a logistic regression model.

#  Reg_model <- 1

## Step 2. Specify the number of covariates (up to 10, 1<=d<=10):

#  Num_cov <- 3

## Step 3. Specify a prior distribution function for each theta by imputing 1 or 2,
## 1 for a normal N(mu,s2) with mean mu and variance s2,
## 2 for a gamma Ga(a,b) with mean a/b and variance a/(b*b)
## and give the numerical values of your hyperparameters:
##  For example, you assume that theta_0 follows N(0,1000), please input as "Prior_0 <- c(1, 0, 1000)".
## Note(1): If, for example, the number of covariates is set at 5, please just ignore the entries for
##              theta_6,..,theta_10.
##          Those numerical values do not affect the ESS computations.
## Note(2): Please specify a gamma prior for the precision parameter tau using the final line
##          following those for the covariates.
##          For example, the number of covariates is 5, please specify the prior using "Prior_6" as
##            "Prior_6 <- c(2, 0.001, 0.001)".

#  Prior_0 <- c(1, 0,    1)    # for theta_0
#  Prior_1 <- c(1, 0,    1)    # for theta_1
#  Prior_2 <- c(1, 0,    1)    # for theta_2
#  Prior_3 <- c(1, 0,    1)    # for theta_3
#  Prior_4 <- c(2, 1,    1)    # for theta_4
#  Prior_5 <- c(1, 0, 1000)    # for theta_5
#  Prior_6 <- c(1, 0, 1000)    # for theta_6
#  Prior_7 <- c(1, 0, 1000)    # for theta_7
#  Prior_8 <- c(1, 0, 1000)    # for theta_8
#  Prior_9 <- c(1, 0, 1000)    # for theta_9
#  Prior_10<- c(1, 0, 1000)    # for theta_10
#  Prior_11<- c(1, 0, 1000)    # for theta_11

## Step 4. Set M being a positive integer chosen so that, initially, it is reasonable to assume the prior ESS <= M.
## If M is not sufficiently large, 'NA' returns as a result of the computations.

#  M <- 10

## Step 5. Specify the number of simulations. A suggested value is 5000.
# The user can use NumSims = 10,000 to carry out the most accurate ESS computations.
# The value of NumSims as low at 1,000 may be used to reduce runtime.

# NumSims <- 5000



## Step 6. If you would like to compute ESS of a subvector of (theta_0,...,theta_11),
##         please input "1" in the corresponding elements in the following indicator vectors.
##         This program can compute ESSs of two subvectors of theta at the same time.
##         If you are interested in three or more subvectors, please repeat the computations
##         with different indicator vectors.
##  For example, if you are interested in the first three parameters, theta_0, theta_1, and theta_2,
##  please input 1's as c(    1,     1,     1,     0,     0,     0,     0,     0,     0,     0,     0,     0)

##                theta0,theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8,theta9,theta10,theta11
#  theta_sub1 <- c(    1,     1,     1,     1,     0,     0,     0,     0,     0,     0,     0,     0)
#  theta_sub2 <- c(    0,     0,     0,     0,     1,     0,     0,     0,     0,     0,     0,     0)


#########################################################################################################
###########      End of sample input information.                                               #########
#########################################################################################################




#########################################################################################################
###########      The code below performs the actual Regression Calculation.                     #########
#########################################################################################################

ESS_RegressionCalc <- function( Reg_model, Num_cov,
                                Prior_0, Prior_1, Prior_2, Prior_3, Prior_4, Prior_5,
                                Prior_6, Prior_7, Prior_8, Prior_9, Prior_10, Prior_11,
                                M, NumSims,
                                theta_sub1, theta_sub2
)
{

  ##### Start computing #####

  # Specify the prior means and Dp values under the priors and the Dq0 values under the epsilon-information priors.
  Prior  <- rbind(Prior_0,Prior_1,Prior_2,Prior_3,Prior_4,Prior_5,Prior_6,Prior_7,Prior_8,Prior_9,Prior_10,Prior_11)
  p_mn   <- numeric(1)
  Dp     <- numeric(1)
  Dq0    <- numeric(1)
  c      <- 10000
  for (j in 1:12){
    if (Prior[j,1] == 1){
      p_mn.s <-  Prior[j,2]
      Dp.s   <-  Prior[j,3]^(-1)
      Dq0.s  <- (Prior[j,3]*c)^(-1)
    }
    if (Prior[j,1] == 2){
      p_mn.s <-  Prior[j,2]/Prior[j,3]
      Dp.s   <- (Prior[j,2]  -1)/p_mn.s^2
      Dq0.s  <- (Prior[j,2]/c-1)/p_mn.s^2
    }
    p_mn <- rbind(p_mn,p_mn.s)
    Dp   <- rbind(Dp  ,Dp.s  )
    Dq0  <- rbind(Dq0 ,Dq0.s )
  }
  p_mn <- p_mn[2:13]
  Dp   <-   Dp[2:13]
  Dq0  <-  Dq0[2:13]

  th_ind <- numeric(12)
  dim_th1 <- Num_cov+2
  dim_th2 <- Num_cov+1
  if (Reg_model == 1){
    for (j in 1:dim_th1){
      th_ind[j] <- 1
    }
  }
  if (Reg_model == 2){
    for (j in 1:dim_th2){
      th_ind[j] <- 1
    }
  }
  cov_ind    <- numeric(11)
  dim_linear <- Num_cov+1
  for (j in 1:dim_linear){
    cov_ind[j] <- 1
  }

  # Compute sum_Dp, the trace of the information matrix of the prior p
  sum_Dp     <- sum(Dp*th_ind)
  sum_Dp.s1  <- sum(Dp*th_ind*theta_sub1)
  sum_Dp.s2  <- sum(Dp*th_ind*theta_sub2)

  # Compute sum_Dq0, the trace of the information matrix of the epsilon-information prior q0
  sum_Dq0    <- sum(Dq0*th_ind)
  sum_Dq0.s1 <- sum(Dq0*th_ind*theta_sub1)
  sum_Dq0.s2 <- sum(Dq0*th_ind*theta_sub2)

  # Simulate Monte Carlo samples Y from f(Y)
  DqYMrep.out    <- numeric(M+1)
  DqYMrep.out.s1 <- numeric(M+1)
  DqYMrep.out.s2 <- numeric(M+1)
  for (t in 1:NumSims)
  {
    DqYm.out    <- numeric(M)
    DqYm.out.s1 <- numeric(M)
    DqYm.out.s2 <- numeric(M)
    DqY         <- numeric(1)
    DqY.s1      <- numeric(1)
    DqY.s2      <- numeric(1)
    for (i in 1:M) {
      # Simulate Monte Carlo samples X from Unif(-1,+1)
      # If you would like, you can modify the upper and lower limits of the distributions.
      X1  <- runif(1,min=-1,max=+1)
      X2  <- runif(1,min=-1,max=+1)
      X3  <- runif(1,min=-1,max=+1)
      X4  <- runif(1,min=-1,max=+1)
      X5  <- runif(1,min=-1,max=+1)
      X6  <- runif(1,min=-1,max=+1)
      X7  <- runif(1,min=-1,max=+1)
      X8  <- runif(1,min=-1,max=+1)
      X9  <- runif(1,min=-1,max=+1)
      X10 <- runif(1,min=-1,max=+1)
      X   <- c(1,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)*cov_ind

      if (Reg_model == 1){
        Dq.lin <- X*X*p_mn[Num_cov+2]
        Dq.all <- numeric(12)
        Dq.all[1:dim_th1] <- c(Dq.lin[1:dim_linear],p_mn[Num_cov+2]^(-2)/2)
        Dq     <- sum(Dq.all)
        Dq.s1  <- sum(Dq.all*theta_sub1)
        Dq.s2  <- sum(Dq.all*theta_sub2)
      }
      if (Reg_model == 2){
        pi     <- exp(sum(p_mn[1:11]*X))/(1+exp(sum(p_mn[1:11]*X)))
        pi_pi2 <- pi - pi^2
        Dq     <- sum(X*X*(pi-pi^2))
        Dq.s1  <- sum(X*X*(pi-pi^2)*theta_sub1[1:11])
        Dq.s2  <- sum(X*X*(pi-pi^2)*theta_sub2[1:11])
      }
      DqY            <- DqY        + Dq
      DqY.s1         <- DqY.s1     + Dq.s1
      DqY.s2         <- DqY.s2     + Dq.s2
      DqYm.out[i]    <- sum_Dq0    + DqY
      DqYm.out.s1[i] <- sum_Dq0.s1 + DqY.s1
      DqYm.out.s2[i] <- sum_Dq0.s2 + DqY.s2
    }
    DqYm.out       <- c(sum_Dq0, DqYm.out)
    DqYm.out.s1    <- c(sum_Dq0.s1, DqYm.out.s1)
    DqYm.out.s2    <- c(sum_Dq0.s2, DqYm.out.s2)
    DqYMrep.out    <- rbind(DqYMrep.out,DqYm.out)
    DqYMrep.out.s1 <- rbind(DqYMrep.out.s1,DqYm.out.s1)
    DqYMrep.out.s2 <- rbind(DqYMrep.out.s2,DqYm.out.s2)
  }
  T1  <- NumSims+1
  DqYMrep.out    <- DqYMrep.out[c(2:T1),]
  DqYMrep.out.s1 <- DqYMrep.out.s1[c(2:T1),]
  DqYMrep.out.s2 <- DqYMrep.out.s2[c(2:T1),]
  Dqm.out        <- numeric(M+1)
  Dqm.out.s1     <- numeric(M+1)
  Dqm.out.s2     <- numeric(M+1)
  M1  <- M+1
  for (i in 1:M1) {
    Dqm.out[i]    <- mean(DqYMrep.out[,i])
    Dqm.out.s1[i] <- mean(DqYMrep.out.s1[,i])
    Dqm.out.s2[i] <- mean(DqYMrep.out.s2[,i])
  }

  # Compute the ESS of the whole theta.
  D.m     <- Dqm.out - sum_Dp
  D.min.n <- which(abs(D.m) == min(abs(D.m)))
  D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
  {
    if (D.min.v < 0)       {
      D.min.v.nxt <- D.m[D.min.n+1]
      ESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
    }
    else if (D.min.v > 0)  {
      D.min.v.prv <- D.m[D.min.n-1]
      ESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
    }
    else if (D.min.v == 0) {
      ESS <- D.min.n -1
    }
  }
  # Compute the ESS.1 of subvector 1 of theta.
  D.m.s1     <- Dqm.out.s1 - sum_Dp.s1
  D.min.n.s1 <- which(abs(D.m.s1) == min(abs(D.m.s1)))
  D.min.v.s1 <- D.m.s1[which(abs(D.m.s1) == min(abs(D.m.s1)))]
  {
    if (D.min.v.s1 < 0)       {
      D.min.v.nxt.s1 <- D.m.s1[D.min.n.s1+1]
      ESS.s1 <- D.min.n.s1 - 1 + (-D.min.v.s1 / (-D.min.v.s1 + D.min.v.nxt.s1))
    }
    else if (D.min.v.s1 > 0)  {
      D.min.v.prv.s1 <- D.m.s1[D.min.n.s1-1]
      ESS.s1 <- D.min.n.s1 - 1 - (D.min.v.s1 / (D.min.v.s1 - D.min.v.prv.s1))
    }
    else if (D.min.v.s1 == 0) {
      ESS.s1 <- D.min.n.s1 -1
    }
  }
  # Compute the ESS.2 of subvector 2 of theta.
  D.m.s2     <- Dqm.out.s2 - sum_Dp.s2
  D.min.n.s2 <- which(abs(D.m.s2) == min(abs(D.m.s2)))
  D.min.v.s2 <- D.m.s2[which(abs(D.m.s2) == min(abs(D.m.s2)))]
  {
    if (D.min.v.s2 < 0)       {
      D.min.v.nxt.s2 <- D.m.s2[D.min.n.s2+1]
      ESS.s2 <- D.min.n.s2 - 1 + (-D.min.v.s2 / (-D.min.v.s2 + D.min.v.nxt.s2))
    }
    else if (D.min.v.s2 > 0)  {
      D.min.v.prv.s2 <- D.m.s2[D.min.n.s2-1]
      ESS.s2 <- D.min.n.s2 - 1 - (D.min.v.s2 / (D.min.v.s2 - D.min.v.prv.s2))
    }
    else if (D.min.v.s2 == 0) {
      ESS.s2 <- D.min.n.s2 -1
    }
  }

  ### The prior ESS of the whole theta is  ESS

  ### The prior ESS of subvector 1     is  ESS.s1

  ### The prior ESS of subvector 2     is  ESS.s2
  return( list(ESSwholetheta=ESS, ESSsubvector1=ESS.s1, ESSsubvector2=ESS.s2) )


}  # end of ESS_RegressionCalc function

