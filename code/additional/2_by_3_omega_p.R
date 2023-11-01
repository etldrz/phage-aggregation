source('single_patch.R')

par(mfrow=c(2,3))

omegas <- c(0.05, 0.2, 0.35, 0.5, 0.75); ps <- c(0.1, 0.5, 0.9)
alpha <- 0.1; lambda <- 0.01; theta <- 0.35



boot.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/bootstrapped/"
boot.files.op <- c()

for(omega in c(omegas[1], omegas[5])) {
  for(p in ps){
    locat <- paste("alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
                   lambda, ", omega=", omega, ".txt", sep="")
    
    curr.boot <- paste(boot.file, locat, sep="")
    
    boot.files.op <- c(boot.files.op, curr.boot)
  }
}


omegas.curr <- c(rep(omegas[1], 3), rep(omegas[5], 3))
ps.curr <- rep(ps, 2)

mapply(function(x, omega, p){
  data <- read.table(x, header=TRUE, sep=",")
  plotBaseSimulation(data)
  pred <- baseSimPrediction(alpha, theta, p, lambda, omega)
  s <- pred[,seq(1, ncol(pred), 2)]
  g <- pred[,seq(2, ncol(pred), 2)]
  
  lines(x=burst.sizes.B, y=s, col='black', lwd=1.65)
  lines(x=burst.sizes.B, y=g, col='black', lwd=1.65)
  mtext(paste("omega:", omega, "p:", p, sep=" "))
}, boot.files.op, omegas.curr, ps.curr)


title <- paste("alpha:", alpha, "lambda:", lambda, "theta:", theta, sep=" ")
mtext(title, side=3, line=-2, outer=TRUE)
