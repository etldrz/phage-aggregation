source('single_patch.R')

alpha <- 0.1; thetas <- c(rep(0.08, 3), rep(0.8, 3)); ps <- rep(c(0.1, 0.5, 0.9), 2); lambda <- 0.01; omega <- 0

par(mfrow=c(2, 3))

mapply(function(x, theta, p) {
  data <- read.table(x, header=TRUE, sep=",")
  plotBaseSimulation(data)
  pred <- baseSimPrediction(alpha, theta, p, lambda, omega)
  s <- pred[,seq(1, ncol(pred), 2)]
  g <- pred[,seq(2, ncol(pred), 2)]
  
  lines(x=burst.sizes.B, y=s, col='black')
  lines(x=burst.sizes.B, y=g, col='black')
  mtext <- paste("theta: ", theta, " p: ", p, sep="")
}, boot.files.tp, thetas, ps)
  

