source('single_patch.R')

library(foreach)
library(ggplot2)
library(viridis)

base.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/base/"
boot.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/bootstrapped/"

lambdas <- c(0.05, 0.2, 0.35, 0.5, 0.75); p <- 0.5
alphas <- c(0.1, 0.5, 0.9); omega <- 0; theta <- 0.35

boot.files.la <- c()
base.files.la <- c()

for(lambda in lambdas) {
  for(alpha in alphas){
    locat <- paste("alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
                   lambda, ", omega=", omega, ".txt", sep="")
    
    curr.base <- paste(base.file, locat, sep="")
    curr.boot <- paste(boot.file, locat, sep="")
    
    base.files.la <- c(base.files.la, curr.base)
    boot.files.la <- c(boot.files.la, curr.boot)
    
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

# first.la <- boot.files.la[c(1, 4, 7, 10, 13)]
# second.la <- boot.files.la[c(2, 5, 8, 11, 14)]
# third.la <- boot.files.la[c(3, 6, 9, 12, 15)]
# 
# fit.first.la <- preprocessed(first.la, lambdas, "lambda")
# fit.second.la <- preprocessed(second.la, lambdas, "lambda")
# fit.third.la <- preprocessed(third.la, lambdas, "lambda")
# 
plot.alpha1.la <- cbind(plotFitness(fit.first.la, FALSE), p=p[1])
plot.alpha2.la <- cbind(plotFitness(fit.second.la, FALSE), p=p[2])
plot.alpha3.la <- cbind(plotFitness(fit.third.la, FALSE), p=p[3])

lambda.prediction.data <- seq(from=0.05, to=0.75, by=0.01)
prediction.alpha1 <- complexSimPrediction(alphas[1], theta, p, lambda.prediction.data)
prediction.alpha2 <- complexSimPrediction(alphas[2], theta, p, lambda.prediction.data)
prediction.alpha3 <- complexSimPrediction(alphas[3], theta, p, lambda.prediction.data)


cols=c('darkred', 'darkorange3', 'darkblue')

#simulated points
plot(y=plot.alpha1.la[,1], x=plot.alpha1.la[,4], type='p', ylab="R*", xlab="\U03BB", cex.lab=1.3,
     ylim=c(0, 1), xlim=c(0.05, 0.75), pch=16, cex=1.15, col=cols[1])
points(y=plot.alpha2.la[,1], x=plot.alpha2.la[,4], pch=16, cex=1.15, col=cols[2])
points(y=plot.alpha3.la[,1], x=plot.alpha3.la[,4], pch=16, cex=1.15, col=cols[3])

#error bars
for(i in 1:nrow(plot.alpha1.la)){
  arrows(y0=plot.alpha1.la$lower.quantile[i], x0=plot.alpha1.la$lambda[i], 
         y1=plot.alpha1.la$upper.quantile[i], x1=plot.alpha1.la$lambda[i], length=0, code=3, angle=90,
         col=cols[1], lwd=1.7)
}
for(i in 1:nrow(plot.alpha2.la)){
  arrows(y0=plot.alpha2.la$lower.quantile[i], x0=plot.alpha2.la$lambda[i], 
         y1=plot.alpha2.la$upper.quantile[i], x1=plot.alpha2.la$lambda[i], length=0, code=3, angle=90,
         col=cols[2], lwd=1.7)
}
for(i in 1:nrow(plot.alpha3.la)){
  arrows(y0=plot.alpha3.la$lower.quantile[i], x0=plot.alpha3.la$lambda[i], 
         y1=plot.alpha3.la$upper.quantile[i], x1=plot.alpha3.la$lambda[i], length=0, code=3, angle=90,
         col=cols[3], lwd=1.7)
}

#prediction lines
lines(y=prediction.alpha1, x=lambda.prediction.data, col=cols[1], lty='dashed', lwd=1.7)
lines(y=prediction.alpha2, x=lambda.prediction.data, col=cols[2], lty='dashed', lwd=1.7)
lines(y=prediction.alpha3, x=lambda.prediction.data, col=cols[3], lty='dashed', lwd=1.7)
legend('topleft', legend=c("\U003B1 = 0.1", "\U003B1 = 0.5", "\U003B1 = 0.9"), 
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

