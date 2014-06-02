library(ggplot2)
monobit.test <- function(sequence, level=0.01) {
  n <- length(sequence)
  if (n < 100) {
    warning("It is recommended that each sequence to be tested consist of a minimum of 100 bits")
  }
  Xi <- function(bit) {
    2*bit-1
  }
  conversion <- Reduce(sum, sapply(X=sequence, FUN=Xi))
  statistic <- abs(conversion)/sqrt(n)
  
  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
  
  p.value <- erfc(statistic/sqrt(2))
  
  print(paste("p-value:", p.value, " with level:", level))
  if (p.value < level) {
    print("The sequence is non-random")
  } else {
    print("The sequence is random")
  }
}

otm.test <- function(sequence, M=10, pattern=c(1,1), level=0.01) {
  n <- length(sequence)
  m <- length(pattern)
  
  if (n < 100) {
    warning("It is recommended that each sequence to be tested consist of a minimum of 100 bits")
  }
  
  lambda <- (M-m+1)/2**m
#   if (round(lambda) != 2) {
#     print(paste("lambda: ", lambda))
#     warning("lambda should be approximately 2")
#   }
  
  eta <- lambda / 2
  K <- 5
  
  P <- function(u, eta) exp(-eta)/2**u * sum(sapply(1:u, function(l) choose(u-1, l-1) * eta**l/factorial(l)))
  pi <- sapply(1:(K+1), function(i) P(i-1, eta))
  
  N <- n/M
#   if (N*min(pi) < 5) {
#     print(paste("N: ", N))
#     warning("N should be chosen so that N*(min*π_i) > 5")
#     #N <- ceiling(5/min(pi))
#   }
  
  matcher <- function(pattern, example) { 
    m <- length(pattern) 
    n <- length(example) 
    candidate <- seq.int(length=n-m+1) 
    for (i in seq.int(length=m)) { 
      candidate <- candidate[pattern[i] == example[candidate + i - 1]] 
    } 
    length(candidate)
  }
  
  i <- 1
  v <- numeric(M-m+1)
  while (i < n) {
    match <- matcher(pattern, sequence[i:(i+M-1)])
    i <- i + M
    v[match+1] <- v[match+1] + 1
  }
  
  schisq <- sapply(1:(K+1), FUN=function(i) (v[i] - N*P(i-1, eta))**2 / (N*P(i-1, eta)))
  chisq <- sum(schisq)
  
  p.value <- pgamma(chisq/2, 5/2)

  print(paste("p-value:", p.value, " with level:", level))
  if (p.value < level) {
    print("The sequence is non-random")
  } else {
    print("The sequence is random")
  }
}


rev.test <- function (sequence, level=.01) {
  n <- length(sequence)
  normalized <- sapply(sequence, function(eps) 2*eps - 1)
  
  sums <- c(0, sapply(1:n, function(i) sum(normalized[1:i])), 0)
  
  df <- data.frame(cbind(c(min(sums):-1, 1:max(sums)), sapply(c(min(sums):-1, 1:max(sums)), function(x){length(Filter(function(p) p == x, sums))})))
  J <- length(df[,1])
  
  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
  
  p.values <- sapply(1:J, function(i) erfc(abs(df[i, 2] - J) / sqrt(2*J*(4*abs(df[i, 1])-2))))
  
  print(
    data.frame(
      'State(x)'= df[,1],
      'Count'= df[,2],
      'P-value'=p.values,
      'Conclusion'=sapply(p.values, function(p) ifelse(p < level, "Non-Random", "Random"))
    )
  )
  if (any(p.values < level)) {
    print("The sequence is non-random")
  } else {
    print("The sequence is random")
  }
  
  ggplot(data.frame("x"=c(0:(n+1)), "y"=sums),aes(x=x,y=y)) + 
#     geom_point() + 
    geom_line() + 
    geom_hline(linetype="dashed", colour="#990000") +
    ggtitle("Random walk") +
    theme(plot.title = element_text(lineheight=.8, face="bold"))
}