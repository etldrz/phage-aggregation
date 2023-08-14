equilibriumSimulation <- function(nA, lyse, alpha, theta, p, 
                                  lambda) {


  nBs = 1:nA
  
  treatments = length(nBs)
  
  
  burstSize <- rep(1, treatments)
  
  prediction = nA / (1 + (alpha + lambda) / (theta * p))
  
  reps = 500000
  wG = c()
  wS = c()
  
  for(i in 1:treatments)
  {
    grandchildrenG <- rep(0, reps)
    grandchildrenS <- rep(0, reps)
    
    nB = nBs[i]
    
    sOutcomes = sample(1:3, reps, replace=TRUE, prob=c(lambda, theta * p, alpha)/ (lambda + theta * p + alpha))
    for (k in 1:reps)
    {
      if (sOutcomes[k] == 2) {
        if (runif(1) < lyse) {
          burstSizeS <- rpois(1, nA)
          sOutcomesLyse <- sample(1:3, burstSizeS, replace = TRUE, prob=c(lambda, theta * p, alpha) / (lambda + theta * p + alpha))
          
          sTabLyse <- tabulate(sOutcomesLyse, 3)
          
          
          grandchildrenS[k] <- sum(sTabLyse * c(0, nA + 1, 1))
          
          sOutcomes[k] <- 4
        }
      }
      
    }
    
    sTab = tabulate(sOutcomes, 3)
    wS <- append(wS, (sum(sTab * c(0, nA+1, 1)) + sum(grandchildrenS)) / reps)
    
    gOutcomes = sample(1:4, reps, replace=TRUE, prob=c(lambda, theta * p, theta * (1-p), alpha)/ (lambda + theta + alpha))
    for (j in 1:reps)
    {
      if (gOutcomes[j] == 2 | gOutcomes[j] == 3)
      {
        if (runif(1) < lyse)
        {
          if (gOutcomes[j] == 2)
          {
            burstSizeG <- rpois(1, nA)
          }
          if (gOutcomes[j] == 3)
          {
            burstSizeG <- rpois(1, nB)
          }
          gOutcomesLyse <- sample(1:4, burstSizeG, replace=TRUE, prob=c(lambda, theta * p, theta * (1-p), alpha)/ (lambda + theta + alpha))
          
          gTabLyse <- tabulate(gOutcomesLyse, 4)
          
          grandchildrenG[j] <- sum(gTabLyse * c(0, nA+1, nB + 1, 1))

          gOutcomes[j] <- 5
        }
      }
    }

    
    gTab = tabulate(gOutcomes, 4)
    
    wG <- append(wG, (sum(gTab * c(0,(nA+1), nB+1, 1)) + sum(grandchildrenG)) / reps)
    print(c("Current iteration: " = i))
    
    #prevents the simulation from continuing past where I need
    if (wG[i] > wS[i]) {
      break
    }

  }

  rStar <- zeroHunter(wS, wG)
  return(rStar)
}