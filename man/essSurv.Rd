\name{essSurv}
\alias{essSurv}
\title{
Calculating ESS for time-to-event outcome
}
\description{
Calculates ESS for right-censored time-to-event outcome
}
\usage{
essSurv(shapeParam,scaleParam,m,nsim)
}
\arguments{
  \item{shapeParam}{Shape parameter of the inverse gamma prior}
  \item{scaleParam}{Scale parameter of the inverse gamma prior}
  \item{m}{A positive integer specified as an maximum value in which ESS is searched.}
  \item{nsim}{umber of simulations for numerical approximation}
}

\value{
\item{ess}{Prior effective sample size}
}

\references{
Morita, S., Thall, P. F., and Muller, P. (2010). Evaluating the impact of prior assumptions in Bayesian biostatistics. Stat Biosci, 2, 1-17.

Morita, S., Thall, P. F., and Muller, P. (2008). Determining the effective sample size of a parametric prior. Biometrics, 64, 595-602.

Thall, P. F., Wooten, L. H., Tannir, N. M. (2005). Monitoring event times in early phase clinical trials: some practical issues. Clinical Trials, 2, 467-478.
}
\author{
Jaejoon Song <jjsong2@mdanderson.org>, Satoshi Morita <smorita@kuhp.kyoto-u.ac.jp >
}

\examples{
## Revisiting Example 5 in Morita et al. (2010, Stat Biosci).
## This is a inverse gamma-exponential model
##  with an inverse gamma prior specified as IG(5.348,30.161)
## we can compute the ESS as the following
essSurv(shapeParam=5.348,scaleParam=30.161,m=7,nsim=1000)
}
}

\keyword{ CRM }
