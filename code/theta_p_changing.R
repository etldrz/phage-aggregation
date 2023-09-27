source('single_patch.R')

library(foreach)
library(ggplot2)
library(viridis)

base.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/base/"
boot.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/bootstrapped/"

thetas <- c(0.05, 0.2, 0.35, 0.5, 0.75); ps <- c(0.1, 0.5, 0.9)

alpha <- 0.1; lambda <- 0.01; omega <- 0

boot.files.tp <- c()
base.files.tp <- c()

for(theta in thetas) {
  for(p in ps){
    locat <- paste("alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
                   lambda, ", omega=", omega, ".txt", sep="")

    curr.base <- paste(base.file, locat, sep="")
    curr.boot <- paste(boot.file, locat, sep="")

    base.files.tp <- c(base.files.tp, curr.base)
    boot.files.tp <- c(boot.files.tp, curr.boot)

#     file.create(curr.base)
#     file.create(curr.boot)
#     
#     curr.base.data <- baseSimulation(alpha, theta, p, lambda, omega)
#     print(dim(curr.base.data))
#     
#     write.table(curr.base.data, curr.base, quote=FALSE, sep=",", append=FALSE)
#     
#     s <- curr.base.data[,seq(1, ncol(curr.base.data), 2)]
#     g <- curr.base.data[,seq(2, ncol(curr.base.data), 2)]
#     
#     curr.boot.data <- foreach(i = 1:ncol(g), .combine='cbind') %do% {
#       basicBootstrap(s[,i], g[,i])
#     }
#     print(dim(curr.boot.data))
#     
#     write.table(curr.boot.data, file=curr.boot, append=FALSE, quote=FALSE, sep=",")
  }
}

first.tp <- boot.files.tp[c(1, 4, 7, 10, 13)]
second.tp <- boot.files.tp[c(2, 5, 8, 11, 14)]
third.tp <- boot.files.tp[c(3, 6, 9, 12, 15)]

# fit.first.tp <- preprocessed(first.tp, thetas, "theta")
# fit.second.tp <- preprocessed(second.tp, thetas, "theta")
# fit.third.tp <- preprocessed(third.tp, thetas, "theta")

plot.p1.tp <- cbind(plotFitness(fit.first.tp, FALSE), p=ps[1])
plot.p2.tp <- cbind(plotFitness(fit.second.tp, FALSE), p=ps[2])
plot.p3.tp <- cbind(plotFitness(fit.third.tp, FALSE), p=ps[3])



theta.prediction.data <- seq(from=0.05, to=0.75, by=0.01)
prediction.p1 <- complexSimPrediction(alpha, theta.prediction.data, ps[1], lambda)
prediction.p2 <- complexSimPrediction(alpha, theta.prediction.data, ps[2], lambda)
prediction.p3 <- complexSimPrediction(alpha, theta.prediction.data, ps[3], lambda)



cols=c('darkorchid4', 'darkorange3', 'darkblue')

#simulated points
plot(y=plot.p1.tp[,1], x=plot.p1.tp[,4], type='p', ylab="R*", xlab="\U03B8", cex.lab=1.3,
     ylim=c(0, 1), xlim=c(0.05, 0.75), pch=16, cex=1.15, col=cols[1])
points(y=plot.p2.tp[,1], x=plot.p2.tp[,4], pch=16, cex=1.15, col=cols[2])
points(y=plot.p3.tp[,1], x=plot.p3.tp[,4], pch=16, cex=1.15, col=cols[3])

#error bars
for(i in 1:nrow(plot.p1.tp)){
  arrows(y0=plot.p1.tp$lower.quantile[i], x0=plot.p1.tp$theta[i], 
         y1=plot.p1.tp$upper.quantile[i], x1=plot.p1.tp$theta[i], length=0, code=3, angle=90,
         col=cols[1], lwd=1.7)
}
for(i in 1:nrow(plot.p2.tp)){
  arrows(y0=plot.p2.tp$lower.quantile[i], x0=plot.p2.tp$theta[i], 
         y1=plot.p2.tp$upper.quantile[i], x1=plot.p2.tp$theta[i], length=0, code=3, angle=90,
         col=cols[2], lwd=1.7)
}
for(i in 1:nrow(plot.p3.tp)){
  arrows(y0=plot.p3.tp$lower.quantile[i], x0=plot.p3.tp$theta[i], 
         y1=plot.p3.tp$upper.quantile[i], x1=plot.p3.tp$theta[i], length=0, code=3, angle=90,
         col=cols[3], lwd=1.7)
}

#prediction lines
lines(y=prediction.p1, x=theta.prediction.data, col=cols[1], lty='dashed', lwd=1.7)
lines(y=prediction.p2, x=theta.prediction.data, col=cols[2], lty='dashed', lwd=1.7)
lines(y=prediction.p3, x=theta.prediction.data, col=cols[3], lty='dashed', lwd=1.7)
legend('topleft', legend=c("p = 0.1", "p = 0.5", "p = 0.9"), 
       fill=c(cols[1], cols[2], cols[3]), border='white', bty='n')


line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight') - .14
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {
  
  text(x = line2user(xl, 2), y = line2user(yl, 3), 
       lab, xpd = NA, font = 2, cex = 1.5, adj = c(0, 1))
  
}

addfiglab("A")
