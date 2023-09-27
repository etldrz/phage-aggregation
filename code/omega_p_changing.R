source('single_patch.R')

library(foreach)
library(ggplot2)
library(viridis)

base.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/base/"
boot.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/bootstrapped/"

omegas <- c(0.05, 0.2, 0.35, 0.5, 0.75); ps <- c(0.1, 0.5, 0.9)
alpha <- 0.1; lambda <- 0.01; theta <- 0.35

boot.files.op <- c()
base.files.op <- c()

for(omega in omegas) {
  for(p in ps){
    locat <- paste("alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
                   lambda, ", omega=", omega, ".txt", sep="")
    
    curr.base <- paste(base.file, locat, sep="")
    curr.boot <- paste(boot.file, locat, sep="")
    
    base.files.op <- c(base.files.op, curr.base)
    boot.files.op <- c(boot.files.op, curr.boot)
    
    # file.create(curr.base)
    # file.create(curr.boot)
    # 
    # curr.base.data <- baseSimulation(alpha, theta, p, lambda, omega)
    # print(dim(curr.base.data))
    # 
    # write.table(curr.base.data, curr.base, quote=FALSE, sep=",", append=FALSE)
    # 
    # s <- curr.base.data[,seq(1, ncol(curr.base.data), 2)]
    # g <- curr.base.data[,seq(2, ncol(curr.base.data), 2)]
    # 
    # curr.boot.data <- foreach(i = 1:ncol(g), .combine='cbind') %do% {
    #   basicBootstrap(s[,i], g[,i])
    # }
    # print(dim(curr.boot.data))
    # 
    # write.table(curr.boot.data, file=curr.boot, append=FALSE, quote=FALSE, sep=",")
  }
}

first.op <- boot.files.op[c(1, 4, 7, 10, 13)]
second.op <- boot.files.op[c(2, 5, 8, 11, 14)]
third.op <- boot.files.op[c(3, 6, 9, 12, 15)]

# fit.first.op <- preprocessed(first.op, omegas, "omega")
# fit.second.op <- preprocessed(second.op, omegas, "omega")
# fit.third.op <- preprocessed(third.op, omegas, "omega")

plot.p1.op <- cbind(plotFitness(fit.first.op, FALSE), p=p[1])
plot.p2.op <- cbind(plotFitness(fit.second.op, FALSE), p=p[2])
plot.p3.op <- cbind(plotFitness(fit.third.op, FALSE), p=p[3])

prediction.p1 <- complexSimPrediction(alpha, theta, ps[1], lambda)
prediction.p2 <- complexSimPrediction(alpha, theta, ps[2], lambda)
prediction.p3 <- complexSimPrediction(alpha, theta, ps[3], lambda)

title.txt <- paste("alpha:", alpha,  "theta:", theta, "lambda:", lambda)

cols=c('darkorchid4', 'darkorange3', 'darkblue')

#simulated points
plot(y=plot.p1.op[,1], x=plot.p1.op[,4], type='p', ylab="R*", xlab="\U03C9",
     ylim=c(0, 1), xlim=c(0.05, 0.75), pch=16, cex=1.15, col=cols[1])
points(y=plot.p2.op[,1], x=plot.p2.op[,4], pch=16, cex=1.15, col=cols[2])
points(y=plot.p3.op[,1], x=plot.p3.op[,4], pch=16, cex=1.15, col=cols[3])

#error bars
for(i in 1:nrow(plot.p1.op)){
  arrows(y0=plot.p1.op$lower.quantile[i], x0=plot.p1.op$omega[i], 
         y1=plot.p1.op$upper.quantile[i], x1=plot.p1.op$omega[i], length=0, code=3, angle=90,
         col=cols[1], lwd=1.7)
}
for(i in 1:nrow(plot.p2.op)){
  arrows(y0=plot.p2.op$lower.quantile[i], x0=plot.p2.op$omega[i], 
         y1=plot.p2.op$upper.quantile[i], x1=plot.p2.op$omega[i], length=0, code=3, angle=90,
         col=cols[2], lwd=1.7)
}
for(i in 1:nrow(plot.p3.op)){
  arrows(y0=plot.p3.op$lower.quantile[i], x0=plot.p3.op$omega[i], 
         y1=plot.p3.op$upper.quantile[i], x1=plot.p3.op$omega[i], length=0, code=3, angle=90,
         col=cols[3], lwd=1.7)
}

#prediction lines
lines(y=rep(prediction.p1, 5), x=plot.p1.op$omega, col=cols[1], lty='dashed', lwd=1.7)
lines(y=rep(prediction.p2, 5), x=plot.p1.op$omega, col=cols[2], lty='dashed', lwd=1.7)
lines(y=rep(prediction.p3, 5), x=plot.p1.op$omega, col=cols[3], lty='dashed', lwd=1.7)
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

addfiglab("B")