rStarFinder <- function(wG, wS, accuracy) {
  #R* finder
  #make a slope of a large amount of points between the neg closest to 0 and the pos closest to 0
  #you can then find the the two xes closest to 0 (pos and neg) and R* lives in here.
  
  negG <- 0
  posG <- 0
  negX <- 0
  posX <- 0

  closest <- 5
  
  for (i in 1:length(wG)) {
    if (i == 1) {
      closest <- wS[i] - wG[i]
    }
    else if (wS[i] - wG[i] < closest & wS[i] - wG[i] > 0) {
      closest <- wS[i] - wG[i]
      negG <- wG[i]
      negX <- i
    }
  }
  
  for (i in 1:length(wG)) {
    if (i == 1) {
      closest <- wS[i] - wG[i]
    }
    else if ((wG[i] - wS[i]) < closest & wG[i] - wS[i] > 0) {
      closest <- wS[i] - wG[i]
      posG <- wG[i]
      posX <- i
    }
  }
  
  slope <- (posG - negG) / (posX - negX)
  
  extensionY <- 1:accuracy
  step <- -1*negG + posG
  extensionY <- negG + step/accuracy*extensionY
  
  extensionX <- 1:accuracy
  extensionX <- negX + (posX - negX) / accuracy * extensionX 
  
  ##numeric empty
  extendedS <- 1:accuracy
  extendedS <- wS[negX] + (wS[posX] - wS[negX]) / accuracy * extendedS
  
  rStarNeg <- NaN
  rStarPos <- NaN
  closest <- 0
  for (i in 1:accuracy) {
    if (i == 1) {
      closest <- extendedS[i] - extensionY[i]
    }
    else if ((extendedS[i] - extensionY[i] < closest) & (extendedS[i] - extensionY[i]) > 0) {
      closest <- extendedS[i] - extensionY[i]
      rStarNeg <- extensionX[i]
    }
  }

  for (i in 1:accuracy) {
    if (i == 1) {
      closest <- abs(extensionY[i] - extendedS[i])
    }
    else if (abs(extensionY[i] - extendedS[i]) < closest & extensionY[i] - extendedS[i] > 0) {
      closest <- abs(extendedS[i] - extensionY[i])
      rStarPos <- extensionX[i]
    }
  }
  rStar <- (rStarNeg + rStarPos) / 2
  return(rStar)
}
