### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 1

## Example 1.1

# set.seed(1237)	       # this seed is for exact result shown
m = 100000                     # number of samples to simulate
good = numeric(m)              # initialize for use in loop
for (i in 1:m)
{
  pick = sample(1:100, 5)      # vector of 5 items from ith box
  good[i] = sum(pick <= 90)    # number Good in ith box
}
mean(good == 5)                # approximates P{All Good}


## Example 1.2

n = 1:60                        # vector of room sizes
p = numeric(60)                 # initialize vector, all 0s

for (i in n)                    # index values for loop
{
  q = prod(1 - (0:(i-1))/365)   # P(No match) if i people in room
  p[i] = 1 - q                  # changes ith element of p
}

plot(n, p)                      # plot of p against n


## Example 1.3

# set.seed(1237)                  # this seed for exact result shown
m = 100000;  n = 25               # iterations; people in room
x = numeric(m)                    # vector for numbers of matches
for (i in 1:m)
{
  b = sample(1:365, n, repl=T)    # n random birthdays in ith room
  x[i] = n - length(unique(b))    # no. of matches in ith room
}
mean(x == 0);  mean(x)            # approximates P{X=0}; E(X)
cutp = (0:(max(x)+1)) - .5        # break points for histogram
hist(x, breaks=cutp, prob=T)      # relative freq. histogram


## Example 1.5

n = 30                                         # number of trials
x = 0:n;  sp = x/n                             # n+1 possible outcomes
m.err = 1.96*sqrt(sp*(1-sp)/n)                 # n+1 Margins of error
lcl = sp - m.err                               # n+1 Lower conf. limits
ucl = sp + m.err                               # n+1 Upper conf. limits

pp = .80                                       # pp = P(Success)
prob = dbinom(x, n, pp)                        # distribution vector
cover = (pp >= lcl) & (pp <= ucl)              # vector of 0s and 1s
round(cbind(x, sp, lcl, ucl, prob, cover), 4)  # 4-place printout
sum(dbinom(x[cover], n, pp))                   # total cov. prob. at pp


## Example 1.6

n = 30                                      # number of trials
alpha = .05;  k = qnorm(1-alpha/2)          # conf level = 1-alpha
adj = 0                                     # (2 for Agresti-Coull)

x = 0:n;  sp = (x + adj)/(n + 2*adj)        # vectors of
m.err = k*sqrt(sp*(1 - sp)/(n + 2*adj))     #   length
lcl = sp - m.err                            #   n + 1
ucl = sp + m.err                            #

m = 2000                                    # no. of values of pp
pp = seq(1/n, 1 - 1/n, length=m)            # vectors
p.cov = numeric(m)                          #   of length m

for (i in 1:m)                              # loop (values of pp)
{                                           # for each pp:
  cover = (pp[i] >= lcl) & (pp[i] <= ucl)   #  1 if cover, else 0
  p.rel = dbinom(x[cover], n, pp[i])        #  relevant probs.
  p.cov[i] = sum(p.rel)                     #  total coverage prob.
}
plot(pp, p.cov, type="l", ylim=c(1-4*alpha,1))
lines(c(.01,.99), c(1-alpha,1-alpha))

## ----------


## Problem 1.5

m = 100000;  n = 10;  x = numeric(m)
for (i in 1:m) {perm = sample(1:n, n);  x[i] = sum(1:n==perm)}
cutp = (-1:n) - .5;  hist(x, breaks=cutp, prob=T)
mean(x == 0);  mean(x);  sd(x)
points(0:10, dpois(0:10, 1))  # from note


## Problem 1.19

n = 30;  pp = .2                                # binomial parameters
alpha = .05;  kappa = qnorm(1-alpha/2)          # level is 1 - alpha
adj = 0                                         # 0 for traditional; 2 for Agresti-Coull
x = 0:n;  sp = (x + adj)/(n + 2*adj)
CI.len = 2*kappa*sqrt(sp*(1 - sp)/(n + 2*adj))
Prob = dbinom(x, n, pp);  Prod = CI.len*Prob
round(cbind(x, CI.len, Prob, Prod), 4)          # displays computation
sum(Prod)                                       # expected length








