#' @export

#####################################################################
##
## A wrapper function for ESS_RegressionCalc written by Jaejoon Song
##
#####################################################################


ess <- function(model,label=NULL,prior=NULL,ncov=NULL,m=NULL,nsim=NULL,svec1=NULL,svec2=NULL){

  modelStatement <- switch(EXPR = model,
                           'betaBin' = 'beta-binomial',
                           'Betabin' = 'beta-binomial',
                           'betabin' = 'beta-binomial',
                           'BetaBin' = 'beta-binomial',
                           'gammaEx' = 'gamma-exponential',
                           'Gammaex' = 'gamma-exponential',
                           'gammaex' = 'gamma-exponential',
                           'GammaEx' = 'gamma-exponential',
                           'dirMult' = 'dirichlet-multinomial',
                           'Dirmult' = 'dirichlet-multinomial',
                           'dirmult' = 'dirichlet-multinomial',
                           'DirMult' = 'dirichlet-multinomial',
  )

  if(model %in% c('betaBin','Betabin','betabin','BetaBin')){
    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- as.numeric(prior[2]) + as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the beta",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
  }


  if(model %in% c('gammaEx','Gammaex','gammaex','GammaEx')){
    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the gamma",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
  }


  if(model %in% c('dirMult','Dirmult','dirmult','DirMult')){
    label <- 'alpha1'
    for(i in 2:(length(prior)-1)){# i <- 1
      label <- c(label,paste('alpha',i,sep=""))
    }

    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- sum(as.numeric(prior[2:length(prior)]))

    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the dirichlet",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
  }

  if(model %in% c('linreg','Linreg','logistic','Logistic')){
    default_prior <- list(c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000))

    if(length(svec1)<12){svec1 <- c(svec1,rep(0,(12-length(svec1))))}
    if(length(svec2)<12){svec2 <- c(svec2,rep(0,(12-length(svec2))))}

    if(model=='linreg'){model <- 1}
    if(model=='logistic'){model <- 2}

    for(i in 1:length(prior)){
      prior[[i]][prior[[i]] %in% c('norm','normal','N','Norm','Normal')] <- 1
      prior[[i]][prior[[i]] %in% c('gamma','Gamma','Gam','gam')] <- 2
      prior[[i]] <- as.numeric(prior[[i]])
    }

    for(i in 1:length(default_prior)){
      assign(paste('Prior_',(i-1),sep=''),default_prior[[i]])
    }

    if(length(prior)>0){
      for(i in 1:length(prior)){
        assign(paste('Prior_',(i-1),sep=''),prior[[i]])
      }
    }

    if(is.null(label)){
      svec <- svec1 + svec2
      svec[1] <- 1
      label <- 'theta1'
      for(i in 2:sum(svec)){
        label <- c(label,paste('theta',i,sep=""))
      }
    }

    svec1_label <- svec1[1:length(label)]
    svec2_label <- svec2[1:length(label)]

    svec1_label <- label[which(svec1_label==1)]
    svec2_label <- label[which(svec2_label==1)]

    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    svec1_out_label <- paste("(",paste(svec1_label,collapse=","),")",sep="")
    svec2_out_label <- paste("(",paste(svec2_label,collapse=","),")",sep="")

    myESS <- ESS_RegressionCalc (Reg_model = model, Num_cov = ncov,
                                 Prior_0 = Prior_0, Prior_1 = Prior_1, Prior_2 = Prior_2,
                                 Prior_3 = Prior_3, Prior_4 = Prior_4, Prior_5 = Prior_5,
                                 Prior_6 = Prior_6, Prior_7 = Prior_7, Prior_8 = Prior_8,
                                 Prior_9 = Prior_9, Prior_10 = Prior_10, Prior_11 = Prior_11,
                                 M = m, NumSims = nsim,
                                 theta_sub1=svec1, theta_sub2=svec2 )

    ESSoutput <- list(ESSoverall = myESS$ESSwholetheta,
                      ESSsubvec1 = myESS$ESSsubvector1,
                      ESSsubvec2 = myESS$ESSsubvector2
    )

    modelStatement <- switch(EXPR = model, 'linreg' = 'linear regression',
                             'logistic' = 'logistic regression')
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESSoverall: Overall ESS for the whole vector ",overall_out_label, sep="")))
    cat("\n")
    cat(noquote(paste("ESSsubvector1: ESS for the first sub-vector ",svec1_out_label, sep="")))
    cat("\n")
    cat(noquote(paste("ESSsubvector1: ESS for the second sub-vector ",svec2_out_label, sep="")))
    cat("\n")
    cat("\n")

    ESSoutput
  }

}
