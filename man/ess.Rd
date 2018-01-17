\name{ess}
\alias{ess}
\title{
Determines effective sample size of a parametric prior distribution
}
\description{
Determines effective sample size of a parametric prior distribution
    in Bayesian conjugate models (beta-binomial, gamma-exponential, dirichlet-multinomial), as well as Bayesian linear and logistic regression models.
}
\usage{
ess(model,label,prior,ncov,m,nsim,svec1,svec2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Model specifications. Options are 'betaBin', 'gammaEx', 'dirMult' for conjugate beta-binomial, gamma-exponential, and dirichlet-multinomial models, respectively. For non-conjugate models, options are 'linreg' and 'logistic' for linear regression and logistic regression models, respectively.}
  \item{label}{Optional labeling for hyperparameters.}
  \item{prior}{Prior distribution specification specified as a list. Options are 'beta', 'gamma', 'dirichlet', and 'norm' (for normal prior).}
  \item{ncov}{Number of covariates}
  \item{m}{A positive integer specified as an maximum value in which ESS is searched.}
  \item{nsim}{Number of simulations for numerical approximation (specified only for model = 'linreg' or model = 'logistic').}
  \item{svec1}{Specification of first subvector for calculating ESS}
  \item{svec2}{Specification of second subvector for calculating ESS}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ESSoverall}{Overall ESS for the whole vector}
  \item{ESSsubvec1}{ESS for the first sub-vector}
  \item{ESSsubvec2}{ESS for the second sub-vector}
}
\references{
Morita, S., Thall, P. F., and Muller, P. (2008). Determining the effective sample size of a parametric prior. Biometrics, 64, 595-602.
}
\author{
Jaejoon Song <jjsong2@mdanderson.org>, Satoshi Morita <smorita@kuhp.kyoto-u.ac.jp >
}



\examples{
# Calculating ESS for a beta-binomial model with
# beta(1,2) prior
ess(model='betaBin',prior=c('beta',1,2))

# Calculating ESS for a gamma-exponential model with
# gamma(2,2) prior
ess(model='gammaEx',prior=c('gamma',2,3))

# Calculating ESS for a dirichlet-multinomial model with
# dirichlet(10,15,20) prior
ess(model='dirMult',prior=c('dirichlet',10,15,20))

# Calculating ESS for a linear regression model with
# three covariates, with priors specified as
# beta0 ~ N(0,1); beta1 ~ N(0,.1); beta2 ~ N(0,.2); beta3 ~ N(0,.3); tau ~ Gamma(1,1);
# Smaller nsim = 50 is specified for illustration purposes
# The user can use nsim = 10,000 to carry out the most accurate ESS computations.
# The value of nsim as low as 1,000 may be used to reduce runtime.
ess(model='linreg',label=c('beta0','beta1','beta2','beta3','tau'),
prior=list(c('norm',0,1),c('norm',0,.1),c('norm',0,.2),c('norm',0,.3),c('gamma',1,1)),
ncov=3,m=50,nsim=50,svec1=c(0,1,1,1,0),svec2=c(0,0,0,0,1))


# Calculating ESS for a linear regression model with
# two covariates, with priors specified as
# beta0 ~ N(0,1); beta1 ~ N(0,.1); beta2 ~ N(0,.2); tau ~ Gamma(1,1);
# Smaller nsim = 50 is specified for illustration purposes
# The user can use nsim = 10,000 to carry out the most accurate ESS computations.
# The value of nsim as low as 1,000 may be used to reduce runtime.
ess(model='linreg',label=c('beta0','beta1','beta2','tau'),
prior=list(c('norm',0,1),c('norm',0,.1),c('norm',0,.2),c('gamma',1,1)),
ncov=2,m=50,nsim=50,svec1=c(0,1,1,0),svec2=c(0,0,0,1))

# Calculating ESS for a logistic regression model with
# three covariates, with priors specified as
# beta0 ~ N(0,1); beta1 ~ N(0,1); beta2 ~ N(0,1); beta3 ~ N(0,1)
# Smaller nsim = 50 is specified for illustration purposes
# The user can use nsim = 10,000 to carry out the most accurate ESS computations.
# The value of nsim as low as 1,000 may be used to reduce runtime.
ess(model='logistic',label=c('beta0','beta1','beta2','beta3'),
prior=list(c('norm',0,1),c('norm',0,1),c('norm',0,1),c('norm',0,1)),
ncov=3,m=50,nsim=50,svec1=c(1,0,0,0),svec2=c(0,1,1,1))

}

\seealso{
\link{https://biostatistics.mdanderson.org/softwaredownload/}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ess }% use one of  RShowDoc("KEYWORDS")
