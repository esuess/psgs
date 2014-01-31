### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 6


## Example 6.1

# Preliminaries
# set.seed(1237)
m = 2000;  n = 1:m;  x = numeric(m);  x[1] = 0
alpha = 1;  beta = 0.44

# Simulation
for (i in 2:m)
{
   if (x[i-1]==0)  x[i] = rbinom(1, 1, alpha)
   else            x[i] = rbinom(1, 1, 1 - beta)
}
y = cumsum(x)/n             # Running fractions of Long eruptions

# Results
y[m]       # Fraction of Long eruptions among m.  Same as: mean(x)
a = sum(x[1:(m-1)]==0 & x[2:m]==1);  a  # No. of cycles
m/a                                     # Average cycle length
plot(x[1:20], type="b", xlab="Step", ylab="State")  # Similar to Fig. 6.3

# ---- (New plot, similar to Fig. 6.4)
plot(y, type="l", ylim=c(0,1), xlab="Step", ylab="Proportion Long")

# ---- (New plot, similar to Fig. 6.5; additional output)
acf(x) # Autocorrelation plot
acf(x, plot=F) # Printed output



## Example 6.2

P = matrix(c(    0,   1,
                 .44, .56), nrow=2, ncol=2, byrow=T)
P
P2  = P %*% P;   P2
P4  = P2 %*% P2; P4
P8  = P4 %*% P4; P8
P16 = P8 %*% P8; P16



## Example 6.5

# set.seed(12)
m = 50000
eta = .99; theta = .97          # T|D
gamma = .4024;  delta = .9998   # D|T
d = numeric(m);  d[1] = 0       # vector of D's; start at 0
t = numeric (m)                 # vector of T's

for (n in 2:m)
{
  if (d[n-1]==1)  t[n-1] = rbinom(1, 1, eta)
    else          t[n-1] = rbinom(1, 1, 1 - theta)

  if (t[n-1]==1)    d[n] = rbinom(1, 1, gamma)
    else            d[n] = rbinom(1, 1, 1 - delta)
}

runprop = cumsum(d)/1:m         # vector of prevalences
mean(d[(m/2+1):m])              # Avg. prevalence after burn-in

par(mfrow=c(1,2))               # 2-panel plot, similar to Fig. 6.6
  plot(runprop, type="l", ylim=c(0,.05),
     xlab="Step", ylab="Running Proportion Infected")
  acf(d, ylim=c(-.1,.4), xlim=c(1,10))
par(mfrow=c(1,1))
acf(d, plot=F)

## ----------


## No code for problems in this chapter
