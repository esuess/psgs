### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 9


## Example 9.1

# set.seed(1237)
m = 50000                       # iterations
PI = numeric(m);  PI[1] = .5    # vector for results, initial value
alpha = 1;  beta = 1            # parameters of beta prior
eta = .99;  theta = .97         # sensitivity;  specificity
n = 1000;  A = 49;  B = n - A   # data

for (i in 2:m)
{
  num.x = PI[i-1]*eta;  den.x = num.x + (1-PI[i-1])*(1 - theta)
  X = rbinom(1, A, num.x/den.x)
  num.y = PI[i-1]*(1 - eta);  den.y = num.y + (1-PI[i-1])*theta
  Y = rbinom(1, B, num.y/den.y)
  PI[i] = rbeta(1, X + Y + alpha, n - X - Y + beta)
}

aft.brn = seq(m/2 + 1,m)
mean(PI[aft.brn])
quantile(PI[aft.brn], c(.025, .975))

par(mfrow=c(2,1))
   plot(aft.brn, PI[aft.brn], type="l")
   hist(PI[aft.brn], prob=T)
par(mfrow=c(1,1))


##--(new graphics page)

par(mfrow=c(1,2))
   acf(PI, ylim=c(0, .6))
   plot(1:m, cumsum(PI)/(1:m), type="l", ylim=c(.016, .024))
par(mfrow=c(1,1))



## Example 9.2

# set.seed(1237)
m = 50000                               # iterations
MU = numeric(m);  THETA = numeric(m)    # sampled values
THETA[1] = 1                            # initial value
n = 41;  x.bar = 9.6;  x.var = 2.73^2   # data
mu.0 = 0;  th.0 = 400                   # mu priors
alp.0 = 1/2;  kap.0 = 1/5               # theta priors

for (i in 2:m)
{
  th.up = 1/(n/THETA[i-1] + 1/th.0)
  mu.up = (n*x.bar/THETA[i-1] + mu.0/th.0)*th.up
  MU[i] = rnorm(1, mu.up, sqrt(th.up))

  alp.up = n/2 + alp.0
  kap.up = kap.0 + ((n-1)*x.var + n*(x.bar - MU[i])^2)/2
  THETA[i] = 1/rgamma(1, alp.up, kap.up)
}

# Bayesian point and probability interval estimates
aft.brn = (m/2 + 1):m
mean(MU[aft.brn])            # point estimate of mu
bi.MU = quantile(MU[aft.brn], c(.025,.975));  bi.MU
mean(THETA[aft.brn])         # point estimate of theta
bi.THETA = quantile(THETA[aft.brn], c(.025,.975));  bi.THETA
SIGMA = sqrt(THETA)
mean(SIGMA[aft.brn])         # point estimate of sigma
bi.SIGMA = sqrt(bi.THETA);  bi.SIGMA

par(mfrow=c(2,2))
  plot(aft.brn, MU[aft.brn], type="l")
  plot(aft.brn, SIGMA[aft.brn], type="l")
  hist(MU[aft.brn], prob=T);  abline(v=bi.MU, col="red")
  hist(SIGMA[aft.brn], prob=T);  abline(v=bi.SIGMA, col="red")
par(mfrow=c(1,1))




## Summary Data (To use in program below, comment out program lines
## in which X.bar, X.sd, g, and r are found.)

X.bar = c(91.9,  129.0,  104.1,   75.7,  108.7,  100.2,
          62.6,  107.5,   66.7,  129.1,  106.8,   93.4)
X.sd  = c( 9.96, 10.07,   4.98,  12.16,   5.06,  10.65,
           6.52, 11.05,   9.90,   8.39,   8.99,   8.14)
g = 12;  r = 10

## Note: Results on p234 of text are based on set.seed(443)
## and unrounded SDs, and so differ slightly from results
## using the rounded SDs above.


##--

#Assumes matrix X with g rows (batches), r columns (reps),
#Or provide g-vectors of batch means and SDs as the 2nd line.

set.seed(443)
X.bar = apply(X, 1, mean);  X.sd = apply(X, 1, sd)
g = dim(X)[1];  r = dim(X)[2]    # CORRECTION for compatibility with problems
m = 50000;  b = m/4              # iterations;  burn-in
MU = VAR.BAT = VAR.ERR = numeric(m)

mu.0  = 0;     th.0  = 10^10     # prior parameters for MU
alp.0 = .001;  kap.0 = .001      # prior parameters for VAR.BAT
bta.0 = .001;  lam.0 = .001      # prior parameters for VAR.ERR
MU[1] = 150;   a = X.bar         # initial values

for (k in 2:m)  
{
  alp.up = alp.0 + g/2
  kap.up = kap.0 + sum((a - MU[k-1])^2)/2
  VAR.BAT[k] = 1/rgamma(1, alp.up, kap.up)

  bta.up = bta.0 + r*g/2
  lam.up = lam.0 + (sum((r-1)*X.sd^2) + r*sum((a - X.bar)^2))/2
  VAR.ERR[k] = 1/rgamma(1, bta.up, lam.up)

  mu.up = (VAR.BAT[k]*mu.0 + th.0*sum(a))/(VAR.BAT[k] + g*th.0)  #note CORRECTION
  th.up = th.0*VAR.BAT[k]/(VAR.BAT[k] + g*th.0)                  #note CORRECTION
  MU[k] = rnorm(1, mu.up, sqrt(th.up))

  deno = r*VAR.BAT[k] + VAR.ERR[k]
  mu.a = (r*VAR.BAT[k]*X.bar + VAR.ERR[k]*MU[k])/deno
  th.a = (VAR.BAT[k]*VAR.ERR[k])/deno
  a = rnorm(g, mu.a, sqrt(th.a))  
}

mean(MU[b:m]);  sqrt(mean(VAR.BAT[b:m]));  sqrt(mean(VAR.ERR[b:m]))
bi.MU = quantile(MU[b:m], c(.025,.975))
SIGMA.BAT = sqrt(VAR.BAT);  SIGMA.ERR = sqrt(VAR.ERR)
bi.SG.B = quantile(SIGMA.BAT[b:m], c(.025,.975))
bi.SG.E = quantile(SIGMA.ERR[b:m], c(.025,.975))
ICC = VAR.BAT/(VAR.BAT+VAR.ERR);
bi.ICC = quantile(ICC[b:m], c(.025,.975))
bi.MU;  bi.SG.B;  bi.SG.E;   bi.ICC

par(mfrow=c(2,2))
  hist(MU[b:m], prob=T);        abline(v=bi.MU)
  hist(SIGMA.BAT[b:m], prob=T); abline(v=bi.SG.B)
  hist(SIGMA.ERR[b:m], prob=T); abline(v=bi.SG.E)
  hist(ICC[b:m], prob=T);       abline(v=bi.ICC)
par(mfrow=c(1,1))

## ----------



## Problem 9.8 (additional code for program of Example 9.1)


est.d = density(PI[aft.brn], from=0, to=1);  mx = max(est.d$y)
hist(PI[aft.brn], ylim=c(0, mx), prob=T, col="wheat")
lines(est.d, col="darkgreen")
median(PI[aft.brn]);  est.d$x[est.d$y==mx]



## Problem 9.12

x = c( 8.50,  9.75,  9.75,  6.00,  4.00, 10.75,  9.25, 13.25,
      10.50, 12.00, 11.25, 14.50, 12.75,  9.25, 11.00, 11.00,
       8.75,  5.75,  9.25, 11.50, 11.75,  7.75,  7.25, 10.75,
       7.00,  8.00, 13.75,  5.50,  8.25,  8.75, 10.25, 12.50,
       4.50, 10.75,  6.75, 13.25, 14.75,  9.00,  6.25, 11.75,  6.25)
mean(x)
var(x)
shapiro.test(x)

par(mfrow=c(1,2))
  boxplot(x, at=.9, notch=T, ylab="x",
    xlab = "Boxplot and Stripchart")
  stripchart(x, vert=T, method="stack", add=T, offset=.75, at = 1.2)
  qqnorm(x)
par(mfrow=c(1,1))



## Problem 9.15 (generate simulated data, use seed shown to get data shown in text)

set.seed(1212)
g = 12                                    # number of batches
r = 10                                    # replications per batch
mu = 100;  sg.a = 15;  sg.e = 9           # model parameters
a.dat = matrix(rnorm(g, 0, sg.a), nrow=g, ncol=r)
          # ith batch effect across ith row
e.dat = matrix(rnorm(g*r, 0, sg.e), nrow=g, ncol=r)
          # g x r random item variations
X = round(mu + a.dat + e.dat)             # integer data
X



## Problem 9.17 (generate simulated data, use seed shown to get data shown)

set.seed(1237)
g = 12;  r = 10
mu = 100;  sg.a = 1;  sg.e = 9
a.dat = matrix(rnorm(g, 0, sg.a), nrow=g, ncol=r)
e.dat = matrix(rnorm(g*r, 0, sg.e), nrow=g, ncol=r)
X = round(mu + a.dat + e.dat)
X.bar = apply(X, 1, mean);  X.sd = apply(X, 1, sd)
round(rbind(X.bar, X.sd), 3)



## Problem 9.20 (turnip leaf data in matrix format)

X = matrix(c(3.28, 3.09, 3.03, 3.03,
             3.52, 3.48, 3.38, 3.38,
             2.88, 2.80, 2.81, 2.76,
             3.34, 3.38, 3.23, 3.26), nrow=4, ncol=4, byrow=T)



## Problem 9.21 (dye data in matrix format)

X = matrix(c(1545, 1440, 1440, 1520, 1580,
             1540, 1555, 1490, 1560, 1495,
             1595, 1550, 1605, 1510, 1560,
             1445, 1440, 1595, 1465, 1545,
             1595, 1630, 1515, 1635, 1625,
             1520, 1455, 1450, 1480, 1445), 6, 5, byrow=T)



## Problem 9.22 (mean vector for plastic data)

X.bar = c(218,  182,  177,  174,  208,  186,
          206,  192,  187,  154,  208,  176,
          196,  179,  181,  158,  158,  198,
          160,  178,  148,  194)



## Problem 9.24 (summary data, as in Example 9.3)

X.bar =  c(124.2, 127.8, 119.4, 123.4, 110.6, 130.4, 128.4, 127.6, 122.0, 124.4)
X.sd  =  c(10.57, 14.89, 11.55, 10.14, 12.82,  9.99, 12.97, 12.82, 16.72,  8.53)

g = 10;  r = 5















