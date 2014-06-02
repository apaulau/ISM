#' Take last element from vector
#' 
#' @param v vector
#' @param n number of elements to take
#' 
#' @return last n vector elements
take.last <- function(v, n = 1) {v[((length(v)-n + 1):length(v))]}

#' Take first element from vector
#' 
#' @param v vector
#' @param n number of elements to take
#' 
#' @return first n vector elements
take.first <- function(v, n = 1) {v[1:n]}

#' Monte-Carlo integration
#' 
#' @param f function to integrate
#' @param lowerBounds numeric or vector of lower bounds
#' @param upperBpunds numeric or vector of upper bounds
#' @param n number of tries
#' 
#' @return numeric result of integration
mc.integrate <- function(f, lowerBounds, upperBounds, n = 10000)
{
  make.transform.finite <- function(l, u) return(function(x) x * (u - l) + l)
  make.diff.transform.finite <- function(l, u) return(function(x) (u - l))
  
  v <- function(x)
  {
    t <- do.call(f, lapply(1:length(lowerBounds), function(i) {return(make.transform.finite(lowerBounds[i], upperBounds[i])(x[i]))} )) 
    t <- t* prod(sapply(1:length(lowerBounds), function(i){return(make.diff.transform.finite(lowerBounds[i], upperBounds[i])(x[i]))}))
    
    return(t)
  }
  
  mean(sapply(1:n, function(t) v(runif(length(lowerBounds)))))
}

mc.solve <- function(A, f, p, M, m = 1000, n = 500)
{
  H <- diag(length(f))
  sapply(1:length(f), function(t)
  {
    h <- H[t, ]
    sum(sapply(1:m, function(t)
    {
      chain <- make.markov(p, M)(n)
      Q <- c(ifelse(
        p[take.first(chain)] > 0, 
        h[take.first(chain)] / p[take.first(chain)], 
        0))        
      sapply(2:n, function(k)
      {         
        Q <<- c(Q, ifelse(M[chain[k-1], chain[k]] > 0, Q[k-1] * A[chain[k-1],chain[k]] / M[chain[k-1], chain[k]], 0))
      })              
      
      sum(sapply(1:n, function(k) Q[k] * f[chain[k]]))      
    }
    )) / m    
  })      
}

n <- 100000

fa <- function(x) exp(-x) * sqrt(1+x)
fc <- function(x, y) (x**2 + y**2)

fa <- Vectorize(fa)

isom <-  function(t) ( t - 0.5 )/( (1 - t) * t )
disom <- function(t) ( t^2 - t + 0.5 )/( (t-1)^2 * t^2)

fc1 <- (function(f) return(function(t,p)f( isom(t), isom(p) ) * disom(t) * disom(p)))(fc) 

mc.integrate(