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

prob.ind <- function(p, rnd = runif)
{
  find <- function(pos, prob)
  {
    ifelse(prob[pos] == 0, find(pos + 1, prob), pos)
  }
  t <-  cumsum(p)
  t[which(p == 0)]  <- NaN
  pos <- which(t < rnd(1))
  return(find(ifelse(length(pos) == 0, 1, pos + 1), p))
}

## Markov chain #################################################################################
#' Markov chain generator
#' 
#' @param p vector of aprior probabilities
#' @param m matrix of probabilities
#' 
#' @return FO witch takes only parameter n- number of values to generate
make.markov <- function(p, m)
{
  last <<- -1
  function(n)
  {            
    gen_next <- function() last <<- ifelse(last == -1, prob.ind(p), prob.ind(m[last, ]) )
    sapply(1:n, function(i) 
    {      
      gen_next()
    })
  }
}

#' Monte-Carlo integration
#' 
#' @param f function to integrate
#' @param lowerBounds numeric or vector of lower bounds
#' @param upperBpunds numeric or vector of upper bounds
#' @param n number of tries
#' 
#' @return numeric result of integration
mc.integrate <- function(f, lowerBounds, upperBounds, n = 10000) {
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
  len <- length(f)
  H <- diag(len)
  sapply(1:len, function(i)
  {
    print(paste(1, ":", i/len*100))
    h <- H[i, ]
    sum(sapply(1:m, function(j)
    {
      chain <- make.markov(h, M)(n)
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

# Integral 1
n <- 1000
sum(sapply(ksi, function(k) sqrt(1+rexp(1))))/n

mc.integrate(function(x,y) x^2 + y^2, c(0,0), c(1,2))


A <- as.matrix(read.table("lab04/A.txt", header=F, sep = " ", as.is=TRUE))
f <- scan("lab04/f.txt")

i <- 1
M <- matrix(nrow=37, ncol=37)
while (i <= 37) {
  j <- 1
  s <- sum(A[i,])
  while (j <= 37) {
    M[i, j] <- A[i, j]/s
    j <- j+1
  }
  i <- i+1
}

x <- mc.solve(A, f, seq(0.1, 1, 0.9/37), M, m=100, n=50)