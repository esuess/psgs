### Suess & Trumbo
### Introduction to Probability 
### Simulation and Gibbs Sampling 
### with R
### Chapter 3


## Example 3.1

m = 5000;  a = 0;  b = 1            # constants
w = (b - a)/m                       # width of each rectangle
g = seq(a + w/2, b - w/2, length=m) # m  "grid" points
const = 1/sqrt(2 * pi)              # const. for density function
h = const * exp(-g^2 / 2)           # vector of m heights
sum(w * h)                          # total area (approx. prob.)


## Example 3.2

# set.seed(12)
m = 500000                   # number or random points
a = 0;  b = 1                # interval endpoints
w = (b - a)/m
u = a + (b - a) * runif(m)   # vector of m random points
h = dnorm(u)                 # hts. of density above rand. points
sum(w * h)                   # approximate probability


## Example 3.3

# set.seed(12)
m = 500000
u = runif(m, 0, 1)
h = runif(m, 0, 0.4)
frac.acc = mean(h < dnorm(u))
0.4*frac.acc


## Example 3.4

# set.seed(1212)
m = 500000                      # number of observations sampled
a = 0;  b = 1                   # interval endpoints
z = rnorm(m)                    # m random obs. from std. normal
mean(z > a & z <= b)            # proportion of obs. in (c, d)


## Example 3.5

#Initialize:
# set.seed(1212)
m = 10000              # total number of tosses
n = 1:m                # vector: n = 1, 2, ..., m; Toss number

#Simulate and Plot:
h = rbinom(m, 1, 1/2)  # vector: H = 0 or 1 each with prob. 1/2
y = cumsum(h)/n        # vector: Proportion of heads
plot (n, y, type="l", ylim=c(0,1))    # Plot trace

#Verify:
Show = cbind(n,h,y)    # matrix: 3 vectors as cols.: n, h, y
Show[1:10, ]           # print first 10 rows of Show
Show[(m-4):m, ]        # print last 5 rows of Show


## Example 3.6

#Initialize:
# set.seed(1237)
m = 10000;  n = 1:m;  alpha = 0.03;  beta = 0.06
w = numeric(m);  w[1] = 0

#Simulate:
for (i in 2:m)  
{
   if (w[i-1]==0)  w[i] = rbinom(1, 1, alpha)
   else            w[i] = rbinom(1, 1, 1 - beta)  
}
y = cumsum(w)/n

#Results:
y[m]
plot(y, type="l",  ylim=c(0,1))


## Example 3.7

# set.seed(1212)
m = 500000;  n = 12
x = runif(m*n, 0, 30);  DTA = matrix(x, m)
x.bar = rowMeans(DTA);  mean(x.bar > 18)


## Example 3.8

m = 100000
a = 0;  b = 3/2;  w = (b - a)/m
x = seq(a + w/2, b-w/2, length=m)
h = x^2
rect.areas = w*h
sum(rect.areas)                    # Riemann

# set.seed(1212)
u = runif(m, a, b)
h = u^2;  y = (b - a)*h
mean(y)                            # Monte Carlo
2*sd(y)/sqrt(m)                    # MC margin of error


## Example 3.9

m = 10000
g = round(sqrt(m))                 # no. of grid pts on each axis
x1 = rep((1:g - 1/2)/g, times=g)   # these two lines give
x2 = rep((1:g - 1/2)/g, each=g)    #   coordinates of grid points
hx = dnorm(x1)*dnorm(x2)
sum(hx[x1 + x2 < 1])/g^2           # Riemann approximation
(pnorm(sqrt(1/2)) - 0.5)^2         # exact value of J

#---

# set.seed(1237)
u1 = runif(m)                      # these two lines give a random
u2 = runif(m)                      #   point in the unit square
hu = dnorm(u1)*dnorm(u2)
hu.acc = hu[u1 + u2 < 1]           # heights above accepted points
m.prime = length(hu.acc); m.prime  # no. of points in triangle
(1/2)*mean(hu.acc)                 # Monte Carlo result
2*(1/2)*sd(hu.acc)/sqrt(m.prime)   # aprx. Marg. of Err. = 2SD(A)

## ----------


## Problem 3.2

x = seq(-2, 2, by=.25)
taylor.7 = 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720
seq.1024  = (1 + x/1024)^1024
exact = exp(x)
round(cbind(x, taylor.7, seq.1024, exact), 4)

##--(Extra: not in text)

x = seq(-2, 2, by=.01)
taylor.7 = 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120 + x^6/720
seq.1024  = (1 + x/1024)^1024
exact = exp(x)
plot(x, exact, type="l", col="yellow", lwd=5)
lines(x, taylor.7, col="red")
lines(x, seq.1024, col="blue", lty="dashed", lwd=2)


## Problem 3.10

m = 500000
u = 10*runif(m);  v = 5 + rnorm(m)
t = u + v
mean(t > 15);  mean(t);  sd(t)
hist(t)


## Problem 3.11

m = 20000
u = runif(m); y = sqrt(u)
acc = rbinom(m, 1, y); x = y[acc==T]
mean(x); sd(x); mean(x < 1/2); mean(acc)

hist(x, prob=T, ylim=c(0,3), col="wheat")
lines(c(0,1), c(0, 3), lty="dashed", lwd=2, col="darkgreen")
xx = seq(0, 1, len=1000)
lines(xx, dbeta(xx, 3, 1), lwd=2, col="blue")


## Problem 3.14

m = 5000;  n = 1:m
u = runif(m);  h = dnorm(u)
j = cumsum(h)/n
plot(n, j, type="l", ylim=c(0.32, 0.36))
abline(h=0.3413, col="blue");  j[m]



## Problem 3.18

set.seed(1212)
m = 500000 # total number of tosses
n = 1:m # vector: n = 1, 2, ..., m; Toss number

# Code from Example 3.5: vectorized (no explicit loop)
t1 = Sys.time()
h = rbinom(m, 1, 1/2) # vector: H = 0 or 1 each with prob. 1/2
y = cumsum(h)/n # vector: Proportion of heads
t2 = Sys.time(); t2 - t1

Show = cbind(n,h,y);  Show[(m-4):m, ]


 set.seed(1212)
 m = 500000 # total number of tosses
 n = 1:m # vector: n = 1, 2, ..., m; Toss number

 # First block: One operation inside loop
 t1 = Sys.time()
 h = numeric(m)
 for (i in 1:m)  {h[i] = rbinom(1, 1, 1/2)}
 y = cumsum(h)/m
 t2 = Sys.time(); t2 - t1

 Show = cbind(n,h,y);  Show[(m-4):m, ]


  set.seed(1212)
  m = 500000 # total number of tosses
  n = 1:m # vector: n = 1, 2, ..., m; Toss number

  # Second block: More operations inside loop
  t1 = Sys.time()
  y = numeric(m);  h = numeric(m)
  for (i in 1:m)  {
    if (i==1)
      {b = rbinom(1, 1, 1/2); h[i] = y[i] = b}
    else
      {b = rbinom(1, 1, 1/2); h[i] = b;
          y[i] = ((i - 1)*y[i - 1] + b)/i}  }
  t2 = Sys.time(); t2 - t1

  Show = cbind(n,h,y);  Show[(m-4):m, ]

