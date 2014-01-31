### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 8


## Example 8.5

x = seq(.45, .7, .001)
prior = dbeta(x, 330, 270)
post = dbeta(x, 950, 650)

# Plot similar to Fig. 8.1
# (minor changes for Figures 8.2, 8.3, 8.4, and 8.5)
plot(x, post, type="l", ylim=c(0, 35), lwd=2,
      xlab="Proportion in Favor", ylab="Density")
lines(x, prior)

## ----------



## Problem 8.2

alpha = 1:2000       # trial values of alpha
beta = .818*alpha    # corresponding values of beta

# Vector of probabilities for interval (.51, .59)
prob = pbeta(.59, alpha, beta) - pbeta(.51, alpha, beta)
prob.err = abs(.95 - prob)  # errors for probabilities

# Results: Target parameter values
t.al = alpha[prob.err==min(prob.err)]
t.be = round(.818*t.al)
t.al; t.be

# Checking: Achieved mean and probability
a.mean = t.al/(t.al + t.be)
a.mean
a.prob = pbeta(.59, t.al, t.be) - pbeta(.51, t.al, t.be)
a.prob



## Problem 8.3
# Plot similar to Fig. 8.6

alpha = c(.5, 1, 1.2, 2, 5);  beta = alpha
op = par(no.readonly = TRUE)  # records existing parameters
par(mfrow=c(5, 5))            # formats 5 x 5 matrix of plots
par(mar=rep(2, 4), pty="m")   # sets margins
x = seq(.001, .999, .001)

for (i in 1:5)
   {
   for (j in 1:5) {
     top = .2 + 1.2 * max(dbeta(c(.05, .2, .5, .8, .95),
        alpha[j], beta[i]))
     plot(x,dbeta(x, alpha[i], beta[j]),
        type="l", ylim=c(0, top), xlab="", ylab="",
        main=paste("BETA(",alpha[j],",", beta[i],")", sep="")) }
   }
par(op)                       # restores former parameters 



## Problem 8.5

x = 620;  n = 1000                             # data
m = 10000;  pp = seq(0, 1, length=m)           # grid points
igd = dbeta(pp, 330, 270) * dbinom(x, n, pp)   # integrand
d = mean(igd);  d                              # denominator

# Results
post.mean = mean(pp*igd)/d;  post.mean
post.pr.bigwin = (1/m)*sum(igd[pp > .6])/d;  post.pr.bigwin
1-pbeta(.6, 950, 650)                          # verify (not in text)

post.cum = cumsum((igd/d)/m)
min(pp[post.cum > .025]);  min(pp[post.cum > .975])
qbeta(c(.025, .975), 950, 650)                 # verify (not in text)



## Problem 8.6

set.seed(1234)
m = 100000
piec = numeric(m); piec[1] = 0.7                 # states of chain

for (i in 2:m) 
{
   piec[i] = piec[i-1]                           # if no jump, keep current pi
   piep = runif(1, piec[i-1]-.05, piec[i-1]+.05) # proposal
   nmtr = dnorm(piep, .55, .02)*dbinom(620, 1000, piep) %%1
   dmtr = dnorm(piec[i-1], .55, .02)*dbinom(620, 1000, piec[i-1])
   r = nmtr/dmtr; acc = (min(r,1) > runif(1))    # accept prop.?
   if(acc) {piec[i] = piep} 
}

pp = piec[(m/2+1):m]                             # after burn-in
quantile(pp, c(.025,.975)); mean(pp > .6)
qbeta(c(.025,.975), 950, 650); 1-pbeta(.6, 950, 650)

hist(pp, prob=T, col="wheat", main="")
  xx = seq(.5, .7, len=1000)
  lines(xx, dbeta(xx, 950, 650), lty="dashed", lwd=2)




## Problem 8.7

pp = seq(.001, .999, .001)     # avoid 'pi' (3.1416)
like = dbinom(620, 1000, pp)
plot(like, type="l");  pp[like==max(like)]



## Problem 8.10

alp = 5; kap = 1
p.lo = seq(.001,.05, .00001)
p.up = .95 + p.lo
q.lo = qgamma(p.lo, alp, kap)
q.up = qgamma(p.up, alp, kap)
long = q.up - q.lo          # avoid confusion with function 'length'
c(q.lo[long==min(long)], q.up[long==min(long)])



## Problem 8.11

r = 900;  s = 1100;  x = 103
nu = 6000:14000;  n = length(nu)
prior = rep(1/n, n);  like = dhyper(x, r, nu-r, s)
denom = sum(prior*like)
post = prior*like/denom;  cumpost = cumsum(post)
c(min(nu[cumpost >= .025]), max(nu[cumpost <= .975]))



















