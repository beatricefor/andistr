##########################################################################################
#ASYMMETRIC NORMAL DISTRIBUTION (expectile distribution) (d,p,q and r functions)
##########################################################################################
#' @import stats

dAND = function(y,mu=0,sigma=1,p=0.5)
{
  if(length(y) == 0) stop("y must be provided.")
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  dens = ifelse(test=y<mu,yes=2*sqrt(1-p)*sqrt(p)/(sigma*sqrt(pi)*(sqrt(p)+sqrt(1-p)))*exp(-(1-p)*(y-mu)^2/sigma^2),
                no=2*sqrt(1-p)*sqrt(p)/(sigma*sqrt(pi)*(sqrt(p)+sqrt(1-p)))*exp(-p*(y-mu)^2/sigma^2))
  return(dens)
}



pAND = function(q,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(q) == 0) stop("q must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")

  a = sqrt(1-p)/(sqrt(p)+sqrt(1-p))
  b = sqrt(p)/(sqrt(p)+sqrt(1-p))

  cdf_and =   ifelse(test = q<mu, yes = 2*b*pnorm(q,mu,sigma/sqrt(2*(1-p))),
              no = (b-a) + 2*a*pnorm(q,mu,sigma/sqrt(2*p)))

  ifelse(test=lower.tail == TRUE,yes=return(cdf_and),no=return(1-(cdf_and)))
}



qAND = function(prob,mu=0,sigma=1,p=0.5,lower.tail=TRUE)
{
  if(length(prob) == 0) stop("prob must be provided.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")
  if(lower.tail != TRUE && lower.tail != FALSE) stop("lower.tail must be TRUE or FALSE.")
  if(sum(prob > 1 | prob < 0) > 0) stop("All elements of prob must be real numbers in [0,1].")


  a = sqrt(1-p)/(sqrt(p)+sqrt(1-p))
  b = sqrt(p)/(sqrt(p)+sqrt(1-p))
  cdf_in_mode = b

  q = sapply(X=prob,FUN=function(prob,mu,sigma,p){ifelse(test=lower.tail == TRUE,yes=prob,no=(1-prob));
    ifelse(test=prob<cdf_in_mode,yes= qnorm(prob/(2*b), mu, sigma/sqrt(2*(1-p))),
           no= qnorm((prob+a-b)/(2*a), mu, sigma/sqrt(2*p)) )},mu=mu,sigma=sigma,p=p)


  return(q)
}
#
#  plot(qnorm,ylim=c(-4,4))
#  prob = seq(0.01,0.99,length.out = 1000)
#  plot(qAND(prob, p=0.1),ylim=c(-4,4))
#

rAND = function(n,mu=0,sigma=1,p=0.5)
{
  if(length(n) == 0) stop("The sample size n must be provided.")
  if(n <= 0 || n%%1 !=0) stop("The sample size n must be a positive integer value.")
  if(sigma <= 0) stop("sigma must be a positive number.")
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1).")
  if(abs(mu) ==Inf) stop("mu must be a finite real number.")

  u = runif(n)
  r = sapply(X=u,FUN=qAND,mu=mu,sigma=sigma,p=p)
  return(r)
}


meanAND = function(mu=0,sigma=1,p=0.5){
  #mean = mu + (sigma/(sqrt(p)+sqrt(1-p)))*((1-2*p)/sqrt(pi*p*(1-p)))
  mean = mu + sigma*((sqrt(1-p)-sqrt(p))/(sqrt(pi*(1-p)*p)))
  return(mean)
}

varAND = function(sigma = 1, p=0.5){
  #var = (sigma^2/(sqrt(p)+sqrt(1-p)))*(1/2 * (sqrt(1-p)/p + sqrt(p)/(1-p)) - (1/(sqrt(p)+sqrt(1-p)))*((1-2*p)^2/pi*p*(1-p)))
  var = (sigma^2/(2*pi*p*(1-p))) * (pi - 2 + (4 - pi)*sqrt(p*(1-p)) )
  return(var)
}

moment_AND = function(sigma = 1, p=0.5){
  moment = sigma*((sqrt(1-p)-sqrt(p))/sqrt(pi*(1-p)*p))
  return(moment)
}

moment2_AND = function(sigma = 1, p = 0.5){
  moment = sigma^2 * (p^(3/2) + (1-p)^(3/2))/(p*(1-p)*(sqrt(p)+sqrt(1-p)))
  return(moment)
}

likAND = function(y,mu=0,sigma=1, p=0.5,loglik=TRUE){
  if(loglik != TRUE && loglik != FALSE) stop("log must be TRUE or FALSE.")
  ifelse(test=loglik == FALSE,yes=return(prod(dAND(y,mu,sigma,p=p))),no=return(sum(log(dAND(y,mu,sigma,p=p)))))
}



# charfunAND = function(y, t, mu = 0, sigma = 1, p = 0.5){
#  P1 = sqrt(pi*sigma^2)/(2*sqrt(1-p))
#  P2 = sqrt(pi*sigma^2)/(2*sqrt(p))
#  #pn1 = (-erfi(-sigma*t/(2*sqrt(p)))+1i)/(2i)
#  #pn2 = (-erfi(-sigma*t/(2*sqrt(1-p)))+1i)/(2i)
#
#  num = ifelse(test = y < mu, yes = P1*2*exp(mu*1i*t - 1/2 *t^2 *(sigma^2/(2*(1-p))))*(-erfi(-sigma*t/(2*sqrt(1-p)))+1)/(2),
#               no =  P2*2*exp(mu*1i*t)*exp(-1/2*t^2*(sigma^2/(2*p)))*(1-(-erfi(-sigma*t/(2*sqrt(p)))+1)/(2)) )
#  chfun = num/(P1+P2)
#  return(chfun)
#
# }
#
# t <- seq(-5, 5, length.out = 501)
# y = rAND(n=501)
# plot(function(y,t)
#   charfunAND(y,t), t, title = "CF of the Gaussian distribution N(0,1)")

