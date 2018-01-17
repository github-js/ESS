\name{essCRM}
\alias{essCRM}
\title{
Calculating ESS for time to event continual reassessment method (CRM)
}
\description{
Calculates ESS for CRM for dose-finding in phase I clinical trials, based on a binary indicator of toxicity
}
\usage{
essCRM(PI,prior,betaSD,target,mRange,numMC)
}
\arguments{
  \item{PI}{A vector of the true toxicity probabilites associated with the doses.}
  \item{prior}{A vector of initial guesses of toxicity probabilities associated with the doses. Must be of same length as PI.}
  \item{betaSD}{Standard deviation of the normal prior of the model parameter.}
  \item{target}{The target dose limiting toxicity (DLT) rate.}
  \item{m}{A positive integer specified as an maximum value in which ESS is searched.}
  \item{nsim}{umber of simulations for numerical approximation}
  \item{obswin}{The observation window with respect to which the maximum tolerated dose (MTD) is defined. Default is obswin=30.}
  \item{rate}{Patient arrival rate: Expected number of arrivals per observation window. Example: obswin=6 and rate=3 means expecting 3 patients arrive in 6 time units. Default is rate=2.}
  \item{accrual}{Patient accrual scheme. Default is accrual="poisson". Alternatively use accrual="fixed" whereby inter-patient arrival is fixed.}
}

\value{
\item{ess}{Prior effective sample size}
}

\references{
Morita, S., Thall, P. F., and Muller, P. (2010). Evaluating the impact of prior assumptions in Bayesian biostatistics. Stat Biosci, 2, 1-17.

Morita, S., Thall, P. F., and Muller, P. (2008). Determining the effective sample size of a parametric prior. Biometrics, 64, 595-602.

O'Quigley J., Pepe M., Fisher, L. (1990).Continual reassessment method: A practical design for phase I clinical trials in cancer. Biometrics, 46, 33-48.
}
\author{
Jaejoon Song <jjsong2@mdanderson.org>, Satoshi Morita <smorita@kuhp.kyoto-u.ac.jp >
}

\examples{
## Revisiting Example 4 in Morita et al. (2010, Stat Biosci).
## We assume five standardized doses to be (.02,.06,.10,.18,.30)
## For the first toxicity scenario given in Table 5 of Morita et al. (2010),
## we have (1mg/kg, 3mg/kg, 6mg/kg, 8mg/kg, 10mg/kg) = (.02,.06,.10,.18,.30)
## Assuming a normal prior for the beta with mean zero and variance of .5,
## we can compute the ESS as the following
essCRM(PI=c(.02,.06,.10,.18,.30),
        prior=c(.02,.06,.10,.18,.30),betaSD=sqrt(.5),
        target=0.2,m=7,nsim=200)
}

\keyword{ CRM }
