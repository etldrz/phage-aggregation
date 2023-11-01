library(mail)
source("equilibriumSimulation.R")

nA = 20
nB <- 15

lambdaChange <- 1:40/100

treatments = length(lambdaChange)


burstSize <- rep(1, treatments)

lyse = 0.001

alpha = 0.05
theta = 0.2
p = 0.3
prediction1 <- length(treatments)
rStar1 <- rep(0, length(treatments))

reps = 500000
wG = rep(0, treatments)
wS = rep(0, treatments)

bugCount <- 0

for(i in 1:treatments)
{
  lambda <- lambdaChange[i]
  
  prediction1 <- nA / (1 + (alpha + lambda) / (theta*p))
  grandchildrenG <- rep(0, reps)
  grandchildrenS <- rep(0, reps)
  
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
  
  wG[i] = (sum(gTab * c(0,(nA+1), nB+1, 1)) + sum(grandchildrenG)) / reps
  
  sOutcomes = sample(1:3, reps, replace=TRUE, prob=c(lambda, theta * p, alpha)/ (lambda + theta * p + alpha))
  for (k in 1:reps)
  {
    if (runif(1) < lyse)
    {
      if (sOutcomes[k] == 2)
      {
        burstSizeS <- rpois(1, nA)
      }
      
      sOutcomesLyse <- sample(1:3, burstSizeS, replace = TRUE, prob=c(lambda, theta * p, alpha) / (lambda + theta * p + alpha))
      
      sTabLyse <- tabulate(sOutcomesLyse, 3)
      
      grandchildrenS[k] <- sum(sTabLyse * c(0, nA + 1, 1))
      
      sOutcomes[k] <- 4
    }
  }
  
  sTab = tabulate(sOutcomes, 3)
  wS[i] = (sum(sTab * c(0, nA+1, 1)) + sum(grandchildrenS)) / reps
  bugCount <- bugCount+1
  rStar1[i] <- equilibriumSimulation(nA, lyse, alpha, theta, p, lambda)
}

plot(lambdaChange, rStar1, type="l", lwd=2, col="deepskyblue2")
points(lambdaChange, prediction1, type="p", lwd=1.5)

sendmail("evandzook@gmail.com", subject="Done?", message="Done!")