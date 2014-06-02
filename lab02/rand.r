library(bitops)

make.rand.lsfr <- function () {
  s <- 1
  polynomial <- function(s) {
    bitOr(
      bitShiftL(
        bitAnd(
          bitXor(
            bitXor(
              bitXor(
                bitXor(
                  bitXor(
                    bitShiftR(s, 31), bitShiftR(s, 30)
                  ), bitShiftR(s, 29)
                ), bitShiftR(s, 27)
              ), bitShiftR(s, 25)
            ), s
          )
          , 1)
        , 31)
      , bitShiftR(s, 1))
  }  
  nextPoly <- function() {
    s <<- polynomial(s)
  }
  
  function(n) {
    sequence <- c()
    i <- 1
    while (i <= n) {
      s1 <- bitAnd(s, 1)
      nextPoly()
      
      s2 <- bitAnd(s, 1)
      if (s1 == 1) {
        sequence[i] <- s2
        i <- i + 1
      }
      nextPoly()
    }
    
    sequence
  }
}