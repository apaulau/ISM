# Linear congruental pseudo-random generator | continuous on [0, 1)
make.rand.lcpr <- function(seed, multiplier, increment, modulus) {
  last <- seed
  function(n) {
    sapply(
      X = numeric(n),
      FUN = function(x) {
        last <<- (multiplier * last + increment) %% modulus
        last / modulus
      }
    )
  }
}

# Linear congruental pseudo-random generator | discrete
make.rand.dlcpr <- function(seed, multiplier, increment, modulus, num) {
  lcpr <- make.rand.lcpr(seed, multiplier, increment, modulus)
  function(n) {
    sapply(
      X = lcpr(n),
      FUN = function(x) {
        for(z in 1:num) {
          if(x <= z/num) {
            return(z - 1)
          }
        }
      }
    )
  }
}

# MacLaren-Marsaglia pseudo-random generator
make.rand.mmpr <- function(rand1, rand2, tableSize = 64) {
  table <- rand1(tableSize)  
  function(n) {    
    t <- 0
    sapply(
      X = trunc(rand2(n) * tableSize + 1),
      FUN = function(s) {               
        t <<- table[s]                   
        table[s] <<- rand1(1)
        return(t)
      }
    )
  }
}

group <- function(x) {
  result <- c()
  for(i in 1:(length(x)-1)) {
    result <- rbind(result, c(x[i], x[i+1]))
  }
  result
}

dfconvert <- function (data) {
  df <- data.frame(group(data))
  colnames(df) <- c("x", "y")
  df
}

cplot <- function (data) {
  print(ggplot(dfconvert(data), aes(x=x, y=y)) + geom_point())
}

moment <- function(x, order = 1, center = FALSE, absolute = FALSE) {
  if (center)
    x <- x - mean(x)
  if (absolute)
    x <- abs(x)
  sum(x^order) / length(x)
}

moment.print <- function(data, maxOrder = 5) {
  res <- c()
  for(i in 1:maxOrder) {
    res <- rbind(res, c(moment(data, order=i), moment(data, order=i, center=T)))
  }
  res <- data.frame(res)
  colnames(res) <- c("Sample", "Central")
  print(res)
}

show.lcpr <- function(seed, multiplier, increment, modulus) {
  rand <- make.rand.lcpr(seed, multiplier, increment, modulus)
  data <- rand(1000) 
  cplot(data)
  moment.print(data)
}

show.dlcpr <- function(seed, multiplier, increment, modulus, num) {
  rand <- make.rand.dlcpr(seed, multiplier, increment, modulus, num)
  data <- rand(1000)
  moment.print(data)
}

show.mmpr <- function(seed1, multiplier1, increment1, modulus1, num1=F, seed2, multiplier2, increment2, modulus2, num2=F) {
  rand1 <- if(num1) make.rand.dlcpr(seed1, multiplier1, increment1, modulus1, num1) else make.rand.lcpr(seed1, multiplier1, increment1, modulus1)
  rand2 <- if(num2) make.rand.dlcpr(seed2, multiplier2, increment2, modulus2, num2) else make.rand.lcpr(seed2, multiplier2, increment2, modulus2)
  
  rand <- make.rand.mmpr(rand1, rand2, 2**10)
  data <- rand(1000) 
  cplot(data)
  moment.print(data)
}

show.extra.mmpr <- function() {
  rand1 <- make.rand.lcpr(17, 13, 5, 1117^12)
  rand2 <- make.rand.lcpr(3, 17, 11, 115^12)
  
  cplot(rand1(5000))
  cplot(rand2(5000))
  
  rand <- make.rand.mmpr(rand1, rand2, 2**10)
  
  cplot(rand(5000))
}

show <- function() {
  show.lcpr(43, 31, 11, 77)
  show.lcpr(17, 171, 11213, 53125)
}