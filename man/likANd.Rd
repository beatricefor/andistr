\name{likAND}
\alias{likAND}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-Likelihood function for the Asymmetric Normal Distribution
}
\description{
Log-Likelihood function for a Three-Parameter Asymmetric Normal Distribution defined in Gerlach and Chen (2015), Waldmann et al (2017) and Philipps (2022) that elicits expectiles, with location parameter equal to \code{mu}, scale parameter \code{sigma} and skewness parameter \code{p}.
}
\usage{
likAND(y, mu = 0, sigma = 1, p = 0.5, loglik = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
observation vector.
}
  \item{mu}{
location parameter \eqn{\mu}.
}
  \item{sigma}{
scale parameter \eqn{\sigma}.
}
  \item{p}{
skewness parameter \eqn{p}.
}
  \item{loglik}{
logical; if TRUE (default), the Log-likelihood is return, if not just the Likelihood.
}
}
\details{
If \code{mu}, \code{sigma} or \code{p} are not specified they assume the default values of 0, 1 and 0.5, respectively, belonging to the Symmetric Normal Distribution denoted by \eqn{AN(0,1,0.5)}.

As discussed in Waldmann et al (2017) and Philipps (2022) we say that a random variable
Y is distributed as an AN with location parameter \eqn{\mu}, scale parameter \eqn{\sigma>0} and skewness parameter \eqn{p} in \eqn{(0,1)}, if its probability density function (pdf) is given by

\deqn{
f(y|\mu,\sigma,p)= \frac{2\sqrt{p(1-p)}}{\sqrt{\pi \sigma^2}\sqrt{p} + \sqrt{1-p}} \exp{-\omega_p(\frac{y-\mu}{\sigma})}
}

where \eqn{\omega_p(.)} is the so called check (or loss) function defined by

\deqn{\omega_p(u) = u^2| p - I_{u<0} |}
with \eqn{I_{.}} denoting the usual indicator function. This  distribution is denoted by \eqn{AN(\mu,\sigma,p)} and \eqn{\mu} it's the \eqn{p}-th expectile of \eqn{y}.

The scale parameter \code{sigma} must be positive and non zero. The skew parameter \code{p} must be between zero and one (0<\code{p}<1).
Then the Log-likelihood function is given by

\deqn{\sum_{i=1}^{n}log(\frac{2\sqrt{p(1-p)}}{\sqrt{\pi \sigma^2}\sqrt p + \sqrt 1-p} \exp{-\omega_p(\frac{y_i-\mu}{\sigma})})
}

The scale parameter \code{sigma} must be positive and non zero. The skew parameter \code{p} must be between zero and one (0<\code{p}<1).

}
\value{
\code{likeALD} returns the Log-likelihood by default and just the Likelihood if \code{loglik = FALSE}.
}

\references{
Gerlach R, Chen CW (2015) Bayesian expected shortfall forecasting incorporating the intraday range. Journal of Financial Econometrics 14(1):128–158


Waldmann E, Sobotka F, Kneib T (2017) Bayesian regularisation in geoadditive expectile regression. Statistics and Computing 27(6):1539–1553

Philipps, C. (2022). The Asymmetric Normal Distribution That Elicits Expectiles. Available at SSRN 4126460.
}
\author{
Beatrice Foroni <\email{beatrice.foroni@uniroma1.it}>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{meanAND}},\code{\link{varAND}}, \code{\link{AND}}
}

\examples{
## Let's compute the log-likelihood for a given sample
library(andistr)
y = andistr:::rAND(n=1000)
loglik = andistr:::likAND(y)

#Changing the true parameters the loglik must decrease
loglik2 = andistr:::likAND(y,mu=10,sigma=2,p=0.3)

loglik;loglik2
if(loglik>loglik2){print("First parameters are Better")}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{asymmetric normal distribution}
\keyword{expectile regression}
\keyword{Log-likelihood}
\keyword{likelihood}
\keyword{AND}
