### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 4

## Example 4.1
#  Scenario 1

# set.seed(1212)
m = 100000;  lam1 = 1/5;  lam2 = 1/4     # constants
x1 = rexp(m, lam1);  x2 = rexp(m, lam2)  # simulation
v = pmin(x1, x2)                         # m-Vector of minimums
mean(v > 5)                              # approximates P{V > 5}
hist(v, prob=T)
xx = seq(0, 15, by=.01);  lines(xx, dexp(xx, 9/20))

## --
#  Scenario 2

# set.seed(1212)
m = 100000;  lam1 = 1/5;  lam2 = 1/4
x1 = rexp(m, lam1);  x2 = rexp(m, lam2)
t = x1 + x2
mean(t > 5)
hist(t, prob=T)



## Example 4.2

# set.seed(12)
m = 100000                 # iterations
n = 5;  lam = 1/4          # constants: system parameters
x = rexp(m*n, lam)         # vector of all simulated data
DTA = matrix(x, nrow=m)    # each row a simulated system
w = apply(DTA, 1, max)     # lifetimes of m systems
mean(w)                    # approximates E(W)
mean(w > 5)                # approximates P{W > 5}
quantile(w, .9)            # approximates 90th percentile

## --

ecdf = (1:m)/m
w.sort = sort(w)      # order statistics for system simulated above
y = rexp(m, lam)      # simulation for single CPU
y.sort = sort(y)      # order statistics for single CPU
plot(w.sort, ecdf, type="l", xlim=c(0,30), xlab="Years")
lines(y.sort, ecdf, lty="dashed")



## Example 4.3

# set.seed(1237)
m = 100000;  n = 10;  mu = 100;  sg = 10
x = rnorm(m*n, mu, sg);  DTA = matrix(x, m)
x.mx = apply(DTA, 1, max);  x.mn = apply(DTA, 1, min)
x.rg = x.mx - x.mn      # vector of m sample ranges
mean(x.rg);  sd(x.rg)
quantile(x.rg, c(.025,.975))
hist(x.rg, prob=T)



## Example 4.4

# set.seed(37)
m = 100000;  n = 5;  mu = 200;  sg = 10
y = rnorm(m*n, mu, sg);  DTA = matrix(y, nrow=m)
samp.mn = rowMeans(DTA)
samp.sd = sqrt(rowSums((DTA - samp.mn)^2)/(n-1))
cor(samp.mn, samp.sd)
A = mean(samp.mn <= 195);  B = mean(samp.sd <= 5);  A;  B;  A*B
mean(samp.mn <= 195 & samp.sd <= 5)
plot(samp.mn, samp.sd, pch=".")



## Example 4.5 (assumes Example 4.4 already run)

plot(samp.mn, samp.sd, pch=".", col="darkgrey")
t.crit = qt(0.975, n-1)
t = sqrt(n)*(samp.mn - mu)/samp.sd
cover = (1:m)[abs(t) < t.crit]
points(samp.mn[cover], samp.sd[cover], pch=".")
xx = seq(min(samp.mn),max(samp.mn),length=1000)
ss = abs(sqrt(n)*(xx-mu)/t.crit)
lines(xx, ss, lwd=2)  # boundary lines



## Example 4.6

# set.seed(12)
m = 100000;  n = 5;  lam = 2
DTA = matrix(rexp(m*n, lam), nrow=m)
samp.mn = rowMeans(DTA)
samp.sd = sqrt(rowSums((DTA - samp.mn)^2)/(n-1))
plot(samp.mn, samp.sd, pch=".");  cor(samp.mn, samp.sd)



## Example 4.7

# set.seed(12)
m = 100000; n = 5; alpha = beta = 0.1
DTA = matrix(rbeta(m*n, alpha, beta), nrow=m)
samp.mn = rowMeans(DTA)
samp.sd = sqrt(rowSums((DTA - samp.mn)^2)/(n-1))
cor(samp.mn, samp.sd)
plot(samp.mn, samp.sd, pch=".")



## Example 4.8

Risk = c(38, 23, 41, 18, 37, 36, 23, 62, 31, 34, 24,
         14, 21, 17, 16, 20, 15, 10, 45, 39, 22, 35,
         49, 48, 44, 35, 43, 39, 34, 13, 73, 25, 27)
Ctrl = c(16, 18, 18, 24, 19, 11, 10, 15, 16, 18, 18,
         13, 19, 10, 16, 16, 24, 13,  9, 14, 21, 19,
          7, 18, 19, 12, 11, 22, 25, 16, 13, 11, 13)
Pair.Diff = Risk - Ctrl

#---(makes graph similar to Fig. 4.7)

plot(ecdf(Pair.Diff), pch=19, ylim=c(-.1,1))
xx = seq(-9,60, by = .01)
lines(xx, pnorm(xx, mean(Pair.Diff), sd(Pair.Diff)), lty="dashed")
stripchart(Pair.Diff, meth="stack", add=T, at=-.07, offset=1/2)

#---(makes graph similar to Fig 4.8)

# set.seed(1237)
n = length(Pair.Diff)                   # number of data pairs
d.bar = mean(Pair.Diff)                 # observed mean of diff's
B = 10000                               # number of resamples
re.x = sample(Pair.Diff, B*n, repl=T)
RDTA = matrix(re.x, nrow=B)             # B x n matrix of resamples
re.mean = rowMeans(RDTA)                # vector of B `d-bar-star's
hist(re.mean, prob=T)                   # hist. of bootstrap dist.
bci = quantile(re.mean, c(.025, .975))  # simple bootstrap CI
alt.bci = 2*d.bar - bci[2:1]            # bootstrap percentile CI
bci; alt.bci

## ----------


## Problem 4.12

# Curve for original example
m = 100000
n = 5;  lam.e = 1/4   
x = rexp(m*n, lam.e);  DTA = matrix(x, nrow=m)
w.e = apply(DTA, 1, max)
mean(w.e);  quantile(w.e, .5);  mean(w.e > 10)
ecdf = (1:m)/m;  w.e.sort = sort(w.e)
plot(w.e.sort, ecdf, type="l", xlim=c(0,40), xlab="Years")

# Overlay curve for part (a)
lam = c(1/2, 1/4, 1/4, 1/4, 1/8)   # As in part (a)
x = rexp(m*n, lam);  DTA = matrix(x, nrow=m, byrow=T)
w.a = apply(DTA, 1, max)
mean(w.a);  quantile(w.a, .5);  mean(w.a > 10)
w.a.sort = sort(w.a)
lines(w.a.sort, ecdf, lwd=2, col="darkblue", lty="dashed")



## Problem 4.24

dh = c(8.50,  9.75,  9.75,  6.00,  4.00, 10.75,  9.25, 13.25, 10.50,
      12.00, 11.25, 14.50, 12.75,  9.25, 11.00, 11.00,  8.75,  5.75,
       9.25, 11.50, 11.75,  7.75,  7.25, 10.75,  7.00,  8.00, 13.75,
       5.50,  8.25,  8.75, 10.25, 12.50,  4.50, 10.75,  6.75, 13.25,
      14.75,  9.00,  6.25, 11.75,  6.25)

##--

B = 10000;  n = length(dh)

# Parameter estimates
dh.bar = mean(dh);  sd.dh = sd(dh)

# Resampling
re.x = rnorm(B*n, dh.bar, sd.dh)
RDTA = matrix(re.x, nrow=B)

# Results
re.mean = rowMeans(RDTA)
hist(re.mean)
bci = quantile(re.mean, c(.025, .975));  bci
2*d.bar - bci[2:1]  



## Problem 4.26

m = 1000;  cover = numeric(m);  B = 1000;  n = 50
for (i in 1:m)  
{
  x = rnorm(n);  re.x = sample(x, B*n, repl=T)
  RDTA = matrix(re.x, nrow=B);  re.mean = rowMeans(RDTA)
  cover[i] = prod(quantile(re.mean, c(.025,.975)))  
}
mean(cover < 0)



## Problem 4.27

## (c)
tau = 7000:15000;  like = dhyper(95, 900, tau-900, 1100)
mle = tau[like==max(like)];  mle
plot(tau, like, type="l");  abline(v=mle, lty="dashed")

## (d) [makes new graph]
# Data
c = 900;  n = 1100;  x = 95

# Extimated population size
tau.hat = floor(c*n/x)

# Resample using estimate
B = 10000;  re.tau = floor(c*n/rhyper(B, c, tau.hat-c, n))

# Histogram and bootstrap confidence intervals
hist(re.tau)
bci = quantile(re.tau, c(.025,.975));  bci  # simple bootstrap
2*tau.hat - bci[2:1]                        # percentile method














