.rwr <-
function(W,P0,gamma) {
    W <- t(W)
    PT <- P0
    k <- 0
    delta <- 1
    while  (delta > 1e-10) {
      PT1 <- (1-gamma)*W
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*P0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      PT <- PT4
      k <- k + 1
    }
    return(PT)
}

