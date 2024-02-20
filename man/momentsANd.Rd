\name{momentsAND}
\alias{meanAND}
\alias{varAND}
\alias{moment_AND}
\alias{moment2_AND}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moments for the Asymmetric Normal Distribution
}
\description{
Mean, variance, first two central moments w.r.t \code{mu} for a Three-Parameter Asymmetric Normal Distribution defined in Gerlach and Chen (2015), Waldmann et al (2017) and Philipps (2022) that elicits expectiles, with location parameter equal to \code{mu}, scale parameter \code{sigma} and skewness parameter \code{p}.
}
\usage{
meanAND(mu = 0, sigma = 1, p = 0.5)
varAND(sigma = 1, p=0.5)
moment_AND(sigma = 1, p=0.5)
moment2_AND(sigma = 1, p=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
location parameter.
}
  \item{sigma}{
scale parameter.
}
  \item{p}{
skewness parameter.
}
}
\details{
If \code{mu}, \code{sigma} or \code{p} are not specified they assume the default values of 0, 1 and 0.5, respectively, belonging to the Symmetric Normal Distribution denoted by \eqn{AN(0,1,0.5)}.

As discussed in Waldmann et al (2017) and Philipps (2022) we say that a random variable
Y is distributed as an AN with location parameter \eqn{\mu}, scale parameter \eqn{\sigma>0} and skewness parameter \eqn{p} in (0,1), if its probability density function (pdf) is given by

\deqn{
f(y|\mu,\sigma,p)= \frac{2\sqrt{p(1-p)}}{\sqrt{\pi \sigma^2}\sqrt{p} + \sqrt{1-p}} \exp{-\omega_p(\frac{y-\mu}{\sigma})}
}

where \eqn{\omega_p(.)} is the so called check (or loss) function defined by

\deqn{\omega_p(u) = u^2| p - I_{u<0} |}
with \eqn{I_{.}} denoting the usual indicator function. This  distribution is denoted by \eqn{AN(\mu,\sigma,p)} and \eqn{\mu} it's the \eqn{p}-th expectile of \eqn{y}.

The scale parameter \code{sigma} must be positive and non zero. The skew parameter \code{p} must be between zero and one (0<\code{p}<1).
}
\value{
\code{meanAND} gives the mean, \code{varAND} gives the variance, \code{moment_AND} gives the \eqn{1}st central moment, i.e., \eqn{E(y-\mu)} and \code{moment2_AND} gives the \eqn{2}nd central moment, i.e., \eqn{E(y-\mu)^2}.
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
\code{\link{AND}},\code{\link{likAND}}, \code{\link{dAND}}, \code{\link{pAND}}
}

\examples{
##Let's compute mean, variance and the first two moments for the Asymmetric Normal Distribution
library(andistr)
andistr:::meanAND(mu=0,sigma=1,p=0.5)
andistr:::meanAND(mu=0,sigma=1,p=0.1)
andistr:::meanAND(mu=0,sigma=1,p=0.9)
andistr:::varAND(sigma=1,p=0.5)
andistr:::moment_AND(sigma=1,p=0.5)
andistr:::moment2_AND(sigma=1,p=0.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{moments}
\keyword{asymmetric normal distribution}
\keyword{expectile regression}
\keyword{AND}
