nA = 20

nBs = 1:nA

treatments = length(nBs)


burstSize <- rep(1, treatments)

lyse = 0.001


alpha = 0.05
theta = 0.2
p = 0.3
lambda = 0.01
prediction = nA / (1 + (alpha + lambda) / (theta * p))

reps = 500000
wG = rep(0, treatments)
wS = rep(0, treatments)


# every time this runs, a new scalar for the net fitness of the generalist phage 
# is used.
for(i in 1:treatments)
{
  #making a change here, by making two grandchildren vectors
  grandchildrenG <- rep(0, reps)
  grandchildrenS <- rep(0, reps)
  
  # changing fitness of phage B.
  nB = nBs[i]
  
  # logs the probabilities of different events occurring to the generalist phage.
  gOutcomes = sample(1:4, reps, replace=TRUE, prob=c(lambda, theta * p, theta * (1-p), alpha)/ (lambda + theta + alpha))
  for (j in 1:reps)
  {
    if (gOutcomes[j] == 2 | gOutcomes[j] == 3)
    {
      if (runif(1) < lyse)
      {
        # if outcome 2 then rpois would be lambda = nA if outcome 3 then rpois would nB
        if (gOutcomes[j] == 2)
        {
          burstSizeG <- rpois(1, nA)
        }
        if (gOutcomes[j] == 3)
        {
          burstSizeG <- rpois(1, nB)
        }
        # then go through gOutcomes again with a size that equals the number found above
        gOutcomesLyse <- sample(1:4, burstSizeG, replace=TRUE, prob=c(lambda, theta * p, theta * (1-p), alpha)/ (lambda + theta + alpha))
        
        # then find fitness the same way as line 64 and 67.
        gTabLyse <- tabulate(gOutcomesLyse, 4)
        
        # store the found fitness into grandchildren[j]
        grandchildrenG[j] <- sum(gTabLyse * c(0, nA+1, nB + 1, 1))
        # then set gOutcomes[j] = 5. either ignore or get rid of this 5
        
        gOutcomes[j] <- 5
      }
    }
  }
  # children of a lyse event would be treated as incoming phage
  # trying to find reproductive success of the lysing phage, this extends over two generations
  # trying to avoid double counting parents. every phage that comes out has its own set
  #   of outcomes it can run through.
  
  
  # continues the log of gOutcomes
  gTab = tabulate(gOutcomes, 4)
  
  # fitness of wG, increasing as i increases
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
}

vec <- 1:nA
for (i in 1:length(wG)) {
  vec[i] <- wS[i] - wG[i]
}
smallestAboveZero <- which(vec == min(vec[vec>0]))

for (i in 1:length(vec)) {
  
}

print(vec)
plot(nBs, wG, type="l", lwd=2, col="chocolate2")
points(nBs, wS, type="l", lwd=2, col="dodgerblue2")

# comparing the math to what the code shows
points(rep(prediction, 2), range(wG), type="l", lty="dotted", lwd=1.5)
