\name{ANd}
\alias{AND}
\alias{dAND}
\alias{pAND}
\alias{qAND}
\alias{rAND}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Expectile Normal Asymmetric Distribution
}
\description{
Density, distribution function, quantile function and random generation for a Three-Parameter Asymmetric Normal Distribution defined in Gerlach and Chen (2015), Waldmann et al (2017) and Philipps (2022) that elicits expectiles, with location parameter equal to \code{mu}, scale parameter \code{sigma} and skewness parameter \code{p}.
}
\usage{
dAND(y, mu = 0, sigma = 1, p = 0.5)
pAND(q, mu = 0, sigma = 1, p = 0.5, lower.tail = TRUE)
qAND(prob, mu = 0, sigma = 1, p = 0.5, lower.tail = TRUE)
rAND(n, mu = 0, sigma = 1, p = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y,q}{
vector of quantiles.
}
 \item{prob}{
vector of probabilities.
}
  \item{n}{
number of observations.
}
  \item{mu}{
location parameter.
}
  \item{sigma}{
scale parameter.
}
  \item{p}{
skewness parameter.
}
  \item{lower.tail}{
logical; if TRUE (default), probabilities are P[X \eqn{\le} x] otherwise, P[X > x].
}
}
\details{
If \code{mu}, \code{sigma} or \code{p} are not specified they assume the default values of 0, 1 and 0.5, respectively, belonging to the Symmetric Normal Distribution denoted by \eqn{AN(0,1,0.5)}.

As discussed in Waldmann et al (2017) and Philipps (2022) we say that a random variable
Y is distributed as an AN with location parameter \eqn{\mu}, scale parameter \eqn{\sigma>0} and skewness parameter \eqn{p} in \eqn{(0,1)}, if its probability density function (pdf) is given by

\deqn{
f(y|\mu,\sigma,p)= \frac{2\sqrt{p(1-p)}}{\sqrt{\pi \sigma^2}\sqrt{p} + \sqrt{1-p}} \exp{-\omega_p(\frac{y-\mu}{\sigma})}
}{ascii}

where \eqn{\omega_p(.)} is the so called check (or loss) function defined by

\deqn{\omega_p(u) = u^2| p - I_{u<0} |}
with \eqn{I_{.}} denoting the usual indicator function. This  distribution is denoted by \eqn{AN(\mu,\sigma,p)} and \eqn{\mu} it's the \eqn{p}-th expectile of \eqn{y}.

The scale parameter \code{sigma} must be positive and non zero. The skew parameter \code{p} must be between zero and one (0<\code{p}<1).
}

\value{
\code{dAND} gives the density, \code{pAND} gives the distribution function, \code{qAND} gives the quantile function, and \code{rAND} generates a random sample.

The length of the result is determined by n for \code{rAND}, and is the maximum of the lengths of the numerical arguments for the other functions \code{dAND}, \code{pAND} and \code{qAND}.
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
\code{\link{meanAND}},\code{\link{varAND}},\code{\link{likAND}}
}

\examples{
## Let's plot the Asymmetric Normal Distribution
##Density
library(andistr)
sseq = seq(-80,80,0.01)
dens = andistr:::dAND(y=sseq,mu=0,sigma=1,p=0.75)
plot(sseq,dens,type="l",lwd=2,col="red",xlab="x",ylab="f(x)",
      main="AN Density function")

## Distribution Function
df = andistr:::pAND(q=sseq,mu=0,sigma=1,p=0.15)
plot(sseq,df,type="l",lwd=2,col="blue",xlab="x",ylab="F(x)",
      main="AN Distribution function", xlim = c(-10,10))
abline(h=1,lty=2)

##Inverse Distribution Function
prob = seq(0,1,length.out = 1000)
idf = andistr:::qAND(prob=prob,mu=0,sigma=1,p=0.40)
plot(prob,idf,type="l",lwd=2,col="gray30",xlab="x",ylab=expression(F^{-1}~(x)))
title(main="AN Inverse Distribution function")
abline(v=c(0,1),lty=2)

#Random Sample Histogram
sample_Data = andistr:::rAND(500, p=0.75, sigma = 1, mu = 0)
hist(sample_Data,breaks = 70,freq = FALSE,ylim=c(0,max(dens)),main="")
title(main="Histogram and True density")
lines(sseq,dens,col="red",lwd=2)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{asymmetric normal distribution}
\keyword{expectile regression}
\keyword{Normal}
\keyword{AND}
