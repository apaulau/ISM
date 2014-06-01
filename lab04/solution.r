n <- 100000

fa <- function(x) exp(-x) * sqrt(1+x)
fc <- function(x, y) (x**2 + y**2)

fa <- Vectorize(fa)

isom <-  function(t) ( t - 0.5 )/( (1 - t) * t )
disom <- function(t) ( t^2 - t + 0.5 )/( (t-1)^2 * t^2)

fc1 <- (function(f) return(function(t,p)f( isom(t), isom(p) ) * disom(t) * disom(p)))(fc) 

ra <- mean(fa(runif(n)))

rc <- mean(
  fc1(runif(10000000), 
      runif(10000000)))

print(ra)
print(integrate(fa, 0, 1)$value)

print("")
print(rb)
print(integrate(fb, -Inf, +Inf)$value)

print(rc)  