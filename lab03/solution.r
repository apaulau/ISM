#' Gemoetric distribution random generator
#' 
#' @param p probability of success 
#' @param runif FO witch produces sequence of uniform values
#' 
#' @return FO witch takes only parameter n- number of values to generate
make.rand.binom <- function(p, m, rnd = runif) { function(n) replicate(n, sum(sapply(rnd(m), function(a) { ifelse(p-a <= 0, 0, 1)})))}

#' Poisson distribution random generator
#' 
#' @param lambda Poisson gistribution parameter
#' @param runif FO with produces sequence of uniform values
#' 
#' @return FO witch takes only parameter n- number of values to generate
make.rand.poisson <- function(lambda, rnd = runif) 
{
  e <- exp(-lambda)
  gen_one <- function()
  {
    a <- 1
    k <- 0
    while(a > e)
    {
      a <- a * rnd(1)
      k <- k + 1
    }
    
    return(k - 1)
  }
  
  function(n) replicate(n, gen_one())
}

#' Normal(Gaussian) distribution random generator
#' 
#' @param m mean
#' @param sigma square root of variance
#' 
#' @param N number of random values to generate normal value
#' @param rnd random generator
#' @param rmean theoretical mean of rnd sequence
#' @param rvar theoretical variance of rnd sequence 
#' 
#' @return FO witch takes only parameter n- number of values to generate
make.rand.norm <- function(m = 0, sigma = 1, N = 64, rnd = runif, rmean = 0.5, rvar = 1/12) {
  function(n) replicate(n, sigma*(sum(rnd(N)) - N*rmean) / (sqrt(rvar*N)) + m) 
}

#' Cauchy random generator
#' 
#' @param loc localisation parameter
#' @param scale scale parameter
#' 
#' @param runif uniform random values generator
#' 
#' @return FO witch takes only parameter n- number of values to generate
make.rand.weibull <- function(lambda, c, rnd = runif) {
  function(n) (-1/lambda * log(rnd(n)))**(1/c)
}

show <- function (sample) {
  print(
    data.frame(
      row.names = 1, c('Mean', 'Variance'), 
      Counted = c(mean(sample), var(sample)) 
    )
  )
  
  plot(ecdf(sample),verticals = TRUE, lty=1, pch=".",
       col.hor = "black", col.vert = "black")
  xx <- unique(sort(c(seq(-3, 2, length = 201), knots(ecdf(sample)))))
  lines(xx, ecdf(sample)(xx), col = "blue")
  
  ggplot(data.frame("sample"=sample), aes(x=sample), geom = 'blank') +
    geom_histogram(aes(y = ..density..),colour="darkgrey", fill="white", alpha = 0.5) +
    labs(color="")
}

show.binom <- function (p=.347, m=37, n=2500) { 
  show(make.rand.binom(p=p, m=m)(n))
}

show.poisson <- function (n = 10000, lambda = 5) {
  show(make.rand.poisson(lambda)(n))
}

show.normal <- function (m=-3, sigma=4, n=2500) {
  show(make.rand.norm(m, sigma)(n))
}

show.weibull <- function (lambda=12, c=23) {
  show(make.rand.weibull(lambda, c)(3700))
}
