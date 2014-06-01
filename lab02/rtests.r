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

otm.test <- function(sequence, m=2, level=0.01) {
  n <- length(sequence)
  if (n < 100) {
    warning("It is recommended that each sequence to be tested consist of a minimum of 100 bits")
  }
  
  B <- c(1, 1)
  N <- 10
  M <- n/N
  
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
  v <- numeric(n/m)
  while (i < n) {
    match <- matcher(B, sequence[i:i+M])
    i <- i + M
    v[match] <- v[match] + 1
  }
  lambda <- (M-m+1)/2**m
  eta <- lambda / 2
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PIIIII
  chisq <- sum(sapply(1:length(v), FUN=function(i) (v[i] - N*pi[i])**2 / N*pi[i]))
  print(chisq)
  p.value <- pgamma(chisq, N/2)
  print(p.value)
  print(paste("p-value:", p.value, " with level:", level))
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