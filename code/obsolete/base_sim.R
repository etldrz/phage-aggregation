reps <- 100000
omega <- 0.5
alpha <- 0.1
theta <- 0.2
p <- 0.3
lambda <- 0.01
nA <- 30
wG <- rep(0, nA)
wS <- rep(0, nA)

prediction_wS <- rep(NA, nA)
prediction_wG <- rep(NA, nA)



iter <- 1:nA



thetaChance <- function(fitness, nB, phageType){
  
  
  thetaP <- which(fitness %in% c(2, 3))
  
  for(t in thetaP){
    if(fitness[t] == 2){ #host A
      burst_size <- rpois(1, nA + 1)
      if(runif(1) < omega){
        outcomes_A <- sample(c(0, 1, nA), burst_size, replace=TRUE, 
                             prob=c(lambda, alpha, theta*p))
        # make this specific to generalist vs specialist
        
        fitness[t] <- sum(outcomes_A)
      } else fitness[t] <- burst_size - 1
    }
    else if(fitness[t] == 3){ #host B
      burst_size <- rpois(1, nB+1)
      if(runif(1) < omega){
        outcomes_B <- sample(c(0, 1, nA, nB), burst_size, replace=TRUE, 
                             prob=c(lambda, alpha, theta * p, theta*(1-p)))
        fitness[t] <- sum(outcomes_B)
      } else fitness[t] <- burst_size - 1
    }
  }
  
  return(fitness)
}


for(i in iter){
  nB <- i
  fit_wG <- sample(c(0, 1, 2, 3), reps, replace=TRUE, prob=c(lambda, alpha, theta*p, theta*(1-p)))
  fit_wS <- sample(c(0, 1, 2), reps, replace=TRUE, prob=c(lambda, alpha, theta*p))
  
  fit_wG <- thetaChance(fit_wG, nB, "generalist")
  fit_wS <- thetaChance(fit_wS, nB, "specialist")
  prediction_wS[i] <- alpha/(alpha + lambda + theta*p) + 
    theta*p/(alpha + lambda + theta*p)*((1 - omega)*nA + 
                                        omega*(nA + 1)*(alpha/(alpha + lambda + theta*p) + 
                                                        nA*theta*p/(alpha + lambda + theta*p)))
  
  prediction_wG[i] <- alpha/(alpha + lambda + theta) + 
    theta*p/(alpha + lambda + theta)*((1 - omega)*nA + 
                                         omega*(nA + 1)*(alpha/(alpha + lambda + theta) + 
                                                           theta/(alpha + lambda + theta) * (p * nA + (1 - p) * nB))) +
     
    theta*(1 - p)/(alpha + lambda + theta)*((1 - omega)*nB + 
                                         omega*(nB + 1)*(alpha/(alpha + lambda + theta) + 
                                                           theta/(alpha + lambda + theta) * (p * nA + (1 - p) * nB)))
  
  
  wG[i] <- mean(fit_wG)
  wS[i] <- mean(fit_wS)
}




plot(iter, wG, col='coral', pch = 16)
points(iter, wS, col='darkblue')
points(iter, prediction_wS, col='darkblue', type="l", lwd=1.5)
points(iter, prediction_wG, col='coral', lwd=1.5, type="l")
legend("topleft", legend=c("WG", "WS"), fill=c('coral', 'darkblue'))
title("WG vs WS")
