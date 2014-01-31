### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 7


## Example 7.1

#Preliminaries
# set.seed(1237)
m = 100000
x = numeric(m);  x[1] = 1

#Simulation
for (i in 2:m)
{
  if (x[i-1] == 1)
     x[i] = sample(1:4, 1, prob=c(.180,.274,.426,.120))
  if (x[i-1] == 2)
     x[i] = sample(1:4, 1, prob=c(.171,.368,.274,.188))
  if (x[i-1] == 3)
     x[i] = sample(1:4, 1, prob=c(.161,.339,.375,.125))
  if (x[i-1] == 4)
     x[i] = sample(1:4, 1, prob=c(.079,.355,.384,.182))
}

#Results
summary(as.factor(x))/m          # Table of proportions
mean(x[1:(m-1)]==2 & x[2:m]==3)  # Est. Proportion of CpG

# Histogram gives same information as bar chart in Fig. 7.1 
hist(x, breaks=0:4+.5, prob=T, xlab="State", ylab="Proportion")


##--Alternative code for 7.1 
##--More 'elegant,' perhaps not as transparent for novice programmers
##--Same result as above for same seed

# set.seed(1237)
m = 100000
x = numeric(m);  x[1] = 1

P = matrix(c(.180,.274,.426,.120,
             .171,.368,.274,.188,
             .161,.339,.375,.125,
             .079,.355,.384,.182), nrow=4, byrow=T)

for (i in 2:m)
  {
  x[i] = sample(1:4, 1, prob=P[x[i-1], ])
  }

summary(as.factor(x))/m          # Table of proportions
mean(x[1:(m-1)]==2 & x[2:m]==3)  # Est. Proportion of CpG
hist(x, breaks=0:4+.5, prob=T, xlab="State", ylab="Proportion")



## Example 7.2 (No code in text; this is transition matrix only)

P = matrix(c(16,   0,   0,   0,   0,   0,
              4,   8,   4,   0,   0,   0,
              1,   4,   4,   4,   2,   1,
              0,   0,   4,   8,   0,   4,
              0,   0,  16,   0,   0,   0,
              0,   0,   0,   0,   0,  16), nrow=6, byrow=T) / 16



## Example 7.3

# set.seed(1237)
m = 1000
d = c(0, sample(c(-1,0,1), m-1, replace=T, c(1/4,1/4,1/2)))
x = cumsum(d)

# Plot similar to Fig. 7.4
plot(x, pch=".", xlab="Step", ylab="State")
lines(c(0,m), c(0,m/4), type="l")



## Example 7.4

# set.seed(1237)
m = 1000;  d = c(0, rnorm(m-1, 0.02, 1));  x = cumsum(d)

# Plot similar to Fig. 7.5 
# (for some seeds, 'fit' to line can be very bad)
plot(x, type="l", xlab="Step", ylab="State")
lines(c(0,m), c(0, 0.02*m), type="l")



## Example 7.5

# set.seed(1212)
m = 50000
d = c(0, runif(m-1, -.1, .1))
x = cumsum(d) %% 1

# Plot similar to Fig. 7.6, but with different seed
hist(x, breaks=10, prob=T, xlab="State", ylab="Proportion")

#--New plot, similar to Fig. 7.7
plot(x[1:1000], pch=".", xlab="Step", ylab="State")



## Example 7.6

# set.seed(1212)
m = 10000
x = numeric(m)
x[1] = .5         # Arbitrary starting value in (0,1)

for (i in 2:m)
{
  x[i] = rbeta(1, .001+3*x[i-1], 3.001-3*x[i-1])
}

plot(x, type="l")



## Example 7.7

# set.seed (1237)
m = 5000;  x = y = numeric(m)
x[1] = y[1] = 0

for (i in 2:m)
{
  x[i] = runif(1, 0, 1-y[i-1])
  y[i] = runif(1, 0, 1-x[i])
}

# Plot similar to Fig. 7.10
plot(x[1:100], y[1:100], type="l", xlim=c(0,1), ylim=c(0,1))

#--New plot, similar to Fig. 7.11 (but with smaller dots)
plot(x, y, pch=".", xlim=c(0,1), ylim=c(0,1))


## Example 7.8

# set.seed(1234);  m = 40000
rho = .8;  sgm = sqrt(1 - rho^2)
xc = yc = numeric(m)             # vectors of state components
xc[1] = -3; yc[1] = 3            # arbitrary starting values
jl = 1; jr = 1                   # left and right limits of proposed jumps

for (i in 2:m)
{
   xc[i] = xc[i-1]; yc[i] = yc[i-1]      # if jump rejected
   xp = runif(1, xc[i-1]-jl, xc[i-1]+jr) # proposed x coord
   yp = runif(1, yc[i-1]-jl, yc[i-1]+jr) # proposed y coord
   nmtr = dnorm(xp)*dnorm(yp, rho*xp, sgm)
   dntr = dnorm(xc[i-1])*dnorm(yc[i-1], rho*xc[i-1], sgm)
   r = nmtr/dntr                         # density ratio
   acc = (min(r, 1) > runif(1))          # jump if acc == T
   if (acc) {xc[i] = xp; yc[i] = yp}
}

x = xc[(m/2+1):m];  y = yc[(m/2+1):m]    # states after burn-in
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)), 4)
mean(diff(x)==0)                         # proportion or proposals rejected
mean(pmax(x,y) >= 1.25)                  # prop. of subj. getting certificates

# 2-panel plot, similar to upper half of Fig. 7.13
par(mfrow = c(1,2), pty="s")
  plot(xc[1:100], yc[1:100], xlim=c(-4,4), ylim=c(-4,4), type="l")
  plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")


## Example 7.9

# set.seed(1235);  m = 20000
rho = .8;  sgm = sqrt(1 - rho^2)
xc = yc = numeric(m)       # vectors of state components
xc[1] = -3;  yc[1] = 3     # arbitrary starting values

for (i in 2:m)
{
   xc[i] = rnorm(1, rho*yc[i-1], sgm)
   yc[i] = rnorm(1, rho*xc[i], sgm)
}
  
x = xc[(m/2+1):m];  y = yc[(m/2+1):m]    # states after burn-in
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)), 4)
best = pmax(x,y);  mean(best >= 1.25)    # proportion getting certif.
summary(best)

# 2-panel plot, similar to lower half of Fig. 7.13
par(mfrow = c(1,2), pty="s")
  hist(best, prob=T, col="wheat", main="")
    abline(v = 1.25, lwd=2, col="red")
  plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
    lines(c(-5, 1.25, 1.25), c(1.25, 1.25, -5), lwd=2, col="red")
par(mfrow = c(1,1), pty="m")


## ----------



## Problem 7.4

m = 10000             # number of runs
step.a = numeric(m)   # steps when absorption occurs
state.a = numeric(m)  # states where absorbed

for (j in 1:m)
{
   x <- 3  # initial state
     # in loop length of x increases to record states visited
   a = 0

   while(a==0)
   {
      i = length(x)  # current state
      # conditions below determine next state
      if (x[i]==1)  x = c(x, 1)
      if (x[i]==2)
         x = c(x, sample(1:6, 1, prob=c(4,8,4,0,0,0)))
      if (x[i]==3)
         x = c(x, sample(1:6, 1, prob=c(1,4,4,4,2,1)))
      if (x[i]==4)
         x = c(x, sample(1:6, 1, prob=c(0,0,4,8,0,4)))
      if (x[i]==5)  x = c(x, 3)
      if (x[i]==6)  x = c(x, 6)
      # condition below determines when absorbed
      if (length(x[x==1 | x==6]) > 0) a = i + 1
   }

   step.a[j] = a  # absorption step for jth run
   state.a[j] = x[length(x)]   # absorption state for jth run
}


hist(step.a)                   # simulated distribution of absorption times
mean(step.a)                   # mean time to absorption
quantile(step.a, .95)          # 95% of runs absorbed by this step
summary(as.factor(state.a))/m  # dist'n of absorption states



## Problem 7.11

# set.seed(1237)
m = 10000
d = sample(c(-1,0,1), m, replace=T, c(1/2,1/4,1/4))
x = numeric(m); x[1] = 0

for (i in 2:m)
{
  x[i] = abs(x[i-1] + d[i])
}

summary(as.factor(x))
cutp=0:(max(x)+1) - .5
hist(x, breaks=cutp, prob=T)



## Problem 7.12

# set.seed(1212)
m = 10000
x = numeric(m);  x[1] = 0

for (i in 2:m)
{
  drift = (2/8)*sign(x[i-1]);  p = c(3/8+drift, 2/8, 3/8-drift)
  x[i] = x[i-1] + sample(c(-1,0,1), 1, replace=T, prob=p)
}

summary(as.factor(x))
par(mfrow=c(2,1))  
  plot(x, type="l")
  cutp = seq(min(x), max(x)+1)-.5;  hist(x, breaks=cutp, prob=T)
par(mfrow=c(1,1))



## Problem 7.16

# set.seed(1212)
m = 5000
e = c(0, 1, 0);  f = c(0, 0, 1)
k = sample(1:3, m, replace=T)
x = numeric(m);  y = numeric(m)
x[1] = 1/2;  y[1] = 1/2

for (i in 2:m)
{
  x[i] = .5*(x[i-1] + e[k[i-1]])
  y[i] = .5*(y[i-1] + f[k[i-1]])
}

plot(x,y,pch=20) 



## Problem 7.17 
# Slightly changed from the text to make a higher-resolution
# image on screen than would be suitable for book production

m = 100000  # Text has 30K
a = c(0, .85, .2, -.15);   b = c(0, .04, -.26, .28)
c = c(0, -.04, .23, .26);  d = c(.16, .85, .22, .24)
e = c(0, 0, 0, 0);         f = c(0, 1.6, 1.6, .44)
p = c(.01, .85, .07, .07)
k = sample(1:4, m, repl=T, p)
h = numeric(m);  w = numeric(m);  h[1] = 0;  w[1] = 0

for (i in 2:m)
{
  h[i] = a[k[i]]*h[i-1] + b[k[i]]*w[i-1] + e[k[i]]
  w[i] = c[k[i]]*h[i-1] + d[k[i]]*w[i-1] + f[k[i]]
}

plot(w, h, pch=".", col="darkgreen")  # Text has 'pch=20'


## Problem 7.19

set.seed(2008)
m = 100000; xc = yc = numeric(m); xc[1] = 3; yc[1] = -3
rho = .8; sgm = sqrt(1 - rho^2); jl = 1.25; jr = .75

for (i in 2:m)
{
    xc[i] = xc[i-1]; yc[i] = yc[i-1] # if no jump
    xp = runif(1, xc[i-1]-jl, xc[i-1]+jr)
    yp = runif(1, yc[i-1]-jl, yc[i-1]+jr)
    
    nmtr.r = dnorm(xp)*dnorm(yp, rho*xp, sgm)
    dntr.r = dnorm(xc[i-1])*dnorm(yc[i-1], rho*xc[i-1], sgm)
    
    nmtr.adj = dunif(xc[i-1], xp-jl, xp+jr)*
                 dunif(yc[i-1], yp-jl, yp+jr)
    dntr.adj = dunif(xp, xc[i-1]-jl, xc[i-1]+jr)*
                 dunif(yp, yc[i-1]-jl, yc[i-1]+jr)
    r = nmtr.r/dntr.r; adj = nmtr.adj/dntr.adj
    
    acc = (min(r*adj, 1) > runif(1))
    if (acc) {xc[i] = xp; yc[i] = yp}
}

x = xc[(m/2+1):m];  y = yc[(m/2+1):m]
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)), 4)
mean(diff(xc)==0);  mean(pmax(x, y) > 1.25)

# Makes 2-panel plot, similar to upr-right & lwr-left in Fig. 7.17
par(mfrow=c(1,2), pty="s")
  jump = diff(unique(x)); hist(jump, prob=T, col="wheat")
  plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow=c(1,1), pty="m")















