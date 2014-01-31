### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 5



## Example 5.1

x = seq(40,80,1)
eta = 1 - pnorm(x, 70, 15)
theta = pnorm(x, 50, 10)
cbind(x, eta, theta)[x >= 54 & x <= 63, ]
plot(1-theta, eta, xlim=c(0,1), ylim=c(0,1), pch=19)
lines(c(0,1),c(1,0))



## Example 5.2

# set.seed(1212)
m = 100000;  n = 250
sens = .99;  spec = .97;  prev = .02
tpos <- prev*sens + (1-prev)*(1-spec)
a = rbinom(m, n, tpos);  t = a/n
lcl.t = t - 1.96*sqrt(t*(1-t)/n)
p = (t + spec - 1)/(sens + spec - 1)
lcl.p = (lcl.t + spec - 1)/(sens + spec - 1)
hist(lcl.p, nclass=5)
mean(p < 0);  mean(lcl.p < 0)



## Example 5.3

# set.seed(1237)
n = 500000
mu.s = 110;  sd.s = 5;  cut.s = 100
sd.x = 1;  cut.x = 101
s = rnorm(n, mu.s, sd.s);  x = rnorm(n, s, sd.x)
n.g  = length(s[s > cut.s])                 # number Good
n.p  = length(x[x > cut.x])                 # number Pass
n.gp = length(x[s > cut.s & x > cut.x])     # number Good & Pass
n.bf = n - n.p - n.g + n.gp                 # number Bad & Fail
pp  = (n - n.g)/n                           # prevalence pi
theta = n.gp/n.g
gamma = n.bf/(n - n.p)
pp;  theta;  gamma

## ----------



## Problem 5.6 (data only)

cre.sens = c(.939, .939, .909, .818, .758, .727, .636, .636, .545,
             .485, .485, .394, .394, .364, .333, .333, .333, .303)
cre.spec = c(.123, .203, .281, .380, .461, .535, .649, .711, .766,
             .773, .803, .811, .843, .870, .891, .894, .896, .909)
b2m.sens = c(.909, .909, .909, .909, .879, .879, .879, .879, .818,
             .818, .818, .788, .788, .697, .636, .606, .576, .576)
b2m.spec = c(.067, .074, .084, .123, .149, .172, .215, .236, .288,
             .359, .400, .429, .474, .512, .539, .596, .639, .676)



## Problem 5.15 (suggested partial code for program to be written)

mu.s = 110;  sd.s = 5;  cut.s = 100   # as in the Example
s = seq(cut.s, mu.s + 5 * sd.s, .001)
int.len = mu.s + 5 * sd.s - cut.s
integrand = dnorm(s, mu.s, sd.s) * (1 - pnorm(cut.x, s, sd.x))
pr.gp = int.len * mean(integrand);  pr.gp 



## Problem 5.19(d)

s = seq(100, 130, 0.001)
numer = 30 * mean(dnorm(s, 110, 5) * dnorm(100.5, s, 1))
denom = dnorm(100.5, 110, 5.099)
numer/denom















