### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 2


## Example 2.1

d = 53           # modulus
a = 20           # multiplier
b = 0            # shift
s = 21           # seed
m = 60           # length of run (counting seed as #1)
r = numeric(m)   # initialize vector for random integers

r[1] = s         # set seed
for (i in 1:(m-1)) r[i+1] = (a * r[i] + b) %% d
                 # generates random integers
r                # list of random integers generated


## Example 2.2

# Initialize
a = 1093;  b = 18257;  d = 86436;  s = 7
m = 1000;  r = numeric(m);  r[1] = s

# Generate
for (i in 1:(m-1)) {r[i+1] = (a*r[i] + b) %% d}
u = (r + 1/2)/d                      # values fit in (0,1)

# Display Results
par(mfrow=c(1,2), pty="s")           # 2 panels in a plot
hist(u, breaks=10, col="wheat")      # left panel
  abline(h=m/10, lty="dashed")
u1 = u[1:(m-1)]; u2 = u[2:m]         # right panel
plot (u1, u2, pch=19)
par(mfrow=c(1,1), pty="m")           # return to default

# then change to m = 50000 and run again


## Example 2.3

a = 65539;  b = 0;  d = 2^31;  s = 10
m = 20000;  r = numeric(m);  r[1] = s
for (i in 1:(m-1)) {r[i+1] = (a*r[i] + b) %% d}
u = (r - 1/2)/(d - 1)
u1 = u[1:(m-2)]; u2 = u[2:(m-1)]; u3 = u[3:m]

par(mfrow=c(1,2), pty="s")
plot(u1, u2, pch=19, xlim=c(0,.1), ylim=c(0,.1))
plot(u1[u3 < .01], u2[u3 < .01], pch=19, xlim=c(0,1), ylim=c(0,1))
par(mfrow=c(1,1), pty="m")


## Example 2.4

# set.seed(1234)
m = 10000
u = runif(m);  x = u^2
xx = seq(0, 1, by=.001)
cut.u = (0:10)/10;  cut.x = cut.u^2

par(mfrow=c(1,2))
  hist(u, breaks=cut.u, prob=T, ylim=c(0,10))
    lines(xx, dunif(xx), col="blue")
  hist(x, breaks=cut.x, prob=T, ylim=c(0,10))
    lines(xx, .5*xx^-.5, col="blue")
par(mfrow=c(1,1))


## Example 2.6

# set.seed(1212)
m = 40000;  z1 = rnorm(m);  z2 = rnorm(m)
t = z1^2 + z2^2;  r = sqrt(t)

hist(t, breaks=30, prob=T, col="wheat")
  tt = seq(0, max(t), length=100);  dens = 0.5*exp(-0.5*tt)
  lines(tt, dens, col="blue")
mean(t);  mean(r);  mean(r < 2)

#--

quad1 = (z1 > 0) & (z2 > 0)   # logical vector
z1.q1 = z1[quad1]; z2.q1 = z2[quad1]
th = (180/pi)*atan(z2.q1/z1.q1)

hist(th, breaks=9, prob=T, col="wheat")
  aa = seq(0, 90, length = 10)
  lines(aa, rep(1/90, 10), col="blue")
sum(quad1)                    # number of hits in quadrant 1

## ----------


## Problem 2.6

a = 1093;  b = 252;  d = 86436;  s = 6
m = 50000;  n = 1:m
r = numeric(m);  r[1] = s
for (i in 1:(m-1)) {r[i+1] = (a*r[i] + b) %% d}

u = (r + 1/2)/d;  f = cumsum(u < .5)/n
plot(n, f, type="l", ylim=c(.49,.51), col="red")  # 'ell', not 1
  abline(h=.5, col="green")

set.seed(1237)            # this seed for exact plot shown in figure
g = cumsum(sample(0:1, m, repl=T))/n
lines(n, g)


## Problem 2.7

m = 10^6;  u = runif(m)      # generate one million from UNIF(0, 1)
u1 = u[1:(m-2)]; u2 = u[2:(m-1)]; u3 = u[3:m]  # 3 dimensions

par(mfrow=c(1,2), pty="s")   # 2 square panels per graph
  plot(u1, u2, pch=".", xlim=c(0,.1), ylim=c(0,.1))
  plot(u1[u3<.01], u2[u3<.01], pch=".", xlim=c(0,1), ylim=c(0,1))
par(mfrow=c(1,1), pty="m")   # restore default plotting 


## Problem 2.8

set.seed(121);  n = 100
par(mfrow=c(1,2), pty="s")     # 2 square panels per graph

# Left Panel
  s = rep(0:9, each=10)/10     # grid points
  t = rep(0:9, times=10)/10
  x = s + runif(n, .01, .09)   # jittered grid points
  y = t + runif(n, .01, .09)
  plot(x, y, pch=19, xaxs="i", yaxs="i", xlim=0:1, ylim=0:1)
    #abline(h = seq(.1, .9, by=.1), col="green")  # grid lines
    #abline(v = seq(.1, .9, by=.1), col="green")

# Right Panel
  x=runif(n);  y = runif(n)    # random points in unit square
  plot(x, y, pch=19, xaxs="i", yaxs="i", xlim=0:1, ylim=0:1)
par(mfrow=c(1,1), pty="m")     # restore default plotting 


## Problem 2.11

# set.seed(1212)
m = 10000;  lam = 1
u = runif(m);  x = -log(u)/lam

cut1 = seq(0, 1, by=.1)              # for hist of u, not plotted
cut2 = -log(cut1)/lam;  cut2[1] = max(x); cut2 = sort(cut2)
hist(x, breaks=cut2, ylim=c(0,lam), prob=T, col="wheat")
  xx = seq(0, max(x), by = .01)
  lines(xx, lam*exp(-lam*xx), col="blue")

mean(x);  1/lam                      # simulated and exact mean
median(x);  qexp(.5, lam)            # simulated and exact median
hist(u, breaks=cut1, plot=F)$counts  # interval counts for UNIF
hist(x, breaks=cut2, plot=F)$counts  # interval counts for EXP


## Problem 2.14

# set.seed(1234)
m = 2*50000; z = numeric(m)
u1 = runif(m/2);  u2 = runif(m/2)
z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)  # half of normal variates
z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)  #    other half
z[seq(1, m, by = 2)] = z1             # interleave
z[seq(2, m, by = 2)] = z2             #    two halves

cut = c(min(z)-.5, seq(-2, 2, by=.5), max(z)+.5)
hist(z, breaks=cut, ylim=c(0,.4), prob=T)
  zz = seq(min(z), max(z), by=.01)
  lines(zz, dnorm(zz), col="blue")

E = m*diff(pnorm(c(-Inf, seq(-2, 2, by=.5), Inf))); E
N = hist(z, breaks=cut, plot=F)$counts; N
Q = sum(((N-E)^2)/E);  Q;  qchisq(c(.025,.975), 9)


## Problem 2.15

# set.seed(1234)
m = 100000;  n = 12
u = runif(m*n)
UNI = matrix(u, nrow=m)
z = rowSums(UNI) - 6

cut = c(min(z)-.5, seq(-2, 2, by=.5), max(z)+.5)
hist(z, breaks=cut, ylim=c(0,.4), prob=T)
  zz = seq(min(z), max(z), by=.01)
  lines(zz, dnorm(zz), col="blue")

E = m*diff(pnorm(c(-Inf, seq(-2, 2, by=.5), Inf))); E
N = hist(z, breaks=cut, plot=F)$counts; N
Q = sum(((N-E)^2)/E);  Q;  qchisq(c(.025,.975), 9)



