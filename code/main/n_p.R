setwd(dirname(rstudioapi::documentPath()))
source('single_patch.R')
setwd("../../")

library(ggplot2)
library(foreach)

base.file <- "~/Desktop/p.a_revisions/varying_p/base/"
boot.file <- "~/Desktop/p.a_revisions/varying_p/boot/"

ps <- c(.1, .5, .9); n <- c(5, 10, 20, 30, 50, 500, 5000)
omega <- 0.05; lambda <- 0.01; alpha <- 0.1; theta <- 0.35
boot.files.psd <- base.files.psd <- c()

for(p in ps) {
 for(popsize in n) {
     locat <- paste("n=", popsize, ", alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
                    lambda, ", omega=", omega, ".txt", sep="")
  
     curr.base <- paste(base.file, locat, sep="")
     curr.boot <- paste(boot.file, locat, sep="")
     

     base.files.psd <- c(base.files.psd, curr.base)
     boot.files.psd <- c(boot.files.psd, curr.boot)

    # file.create(curr.base)
    # file.create(curr.boot)

    # curr.base.data <- baseSimulation(alpha, theta, p, lambda, omega, n=popsize)
    # print(dim(curr.base.data))
    # 
    # # curr.base.data <- read.csv(curr.base)
    # write.table(curr.base.data, curr.base, quote=FALSE, sep=",", append=TRUE,
    #             col.names=FALSE, row.names=FALSE)
    # # rm(curr.base.data)
    # 
    # s <- curr.base.data[,seq(1, ncol(curr.base.data), 2)]
    # g <- curr.base.data[,seq(2, ncol(curr.base.data), 2)]
    # 
    # rm(curr.base.data)
    # curr.boot.data <- foreach(i = 1:ncol(g), .combine='cbind') %do% {
    #   basicBootstrap(s[,i], g[,i])
    # }
    # print(dim(curr.boot.data))
    # 
    # rm(s); rm(g)
    # write.table(curr.boot.data, file=curr.boot, append=F, quote=FALSE, sep=",")
 }
}

.psd <- boot.files.psd[1:7]
second.psd <- boot.files.psd[8:14]
third.psd <- boot.files.psd[15:21]

load("./plots/simple_model_outputs/popsize_p.RData")
# fit.first.psd <- preprocessed(first.psd, n, "n")
# fit.second.psd <- preprocessed(second.psd, n, "n")
# fit.third.psd <- preprocessed(third.psd, n, "n")

plot.1.psd <- cbind(plotFitness(fit.first.psd, F), p=ps[1])
plot.2.psd <- cbind(plotFitness(fit.second.psd, F), p=ps[2])
plot.3.psd <- cbind(plotFitness(fit.third.psd, F), p=ps[3])

plot.1.psd$lgn <- log(plot.1.psd$n, 10)
plot.2.psd$lgn <- log(plot.2.psd$n, 10)
plot.3.psd$lgn <- log(plot.3.psd$n, 10)

prediction.1.psd <- complexSimPrediction(alpha, theta, ps[1], lambda)
prediction.2.psd <- complexSimPrediction(alpha, theta, ps[2], lambda)
prediction.3.psd <- complexSimPrediction(alpha, theta, ps[3], lambda)

title.txt <- paste("alpha:", alpha,  "theta:", theta, "lambda:", lambda)

cols <- c('darkred', 'darkorange3', 'darkblue')

#simulated points
plot(y=plot.1.psd$r.star.mean, x=plot.1.psd$lgn, type='p', ylab="R*", 
     xlab=expression("log"[10]*"(population size)"), cex.lab=1.5,
     ylim=c(0, 1), pch=16, cex=1.15, col=cols[1])
points(y=plot.2.psd$r.star.mean, x=plot.2.psd$lgn, pch=16, cex=1.15, col=cols[2])
points(y=plot.3.psd$r.star.mean, x=plot.3.psd$lgn, pch=16, cex=1.15, col=cols[3])

#error bars
for(i in 1:nrow(plot.1.psd)){
  arrows(y0=plot.1.psd$lower.quantile[i], x0=plot.1.psd$lgn[i], 
         y1=plot.1.psd$upper.quantile[i], x1=plot.1.psd$lgn[i], length=0, code=3, angle=90,
         col=cols[1], lwd=1.7)
}
for(i in 1:nrow(plot.2.psd)){
  arrows(y0=plot.2.psd$lower.quantile[i], x0=plot.2.psd$lgn[i], 
         y1=plot.2.psd$upper.quantile[i], x1=plot.2.psd$lgn[i], length=0, code=3, angle=90,
         col=cols[2], lwd=1.7)
}
for(i in 1:nrow(plot.3.psd)){
  arrows(y0=plot.3.psd$lower.quantile[i], x0=plot.3.psd$lgn[i], 
         y1=plot.3.psd$upper.quantile[i], x1=plot.3.psd$lgn[i], length=0, code=3, angle=90,
         col=cols[3], lwd=1.7)
}

#prediction lines
lines(y=rep(prediction.1.psd, 7), x=plot.1.psd$lgn, col=cols[1], lty='dashed', lwd=1.7)
lines(y=rep(prediction.2.psd, 7), x=plot.2.psd$lgn, col=cols[2], lty='dashed', lwd=1.7)
lines(y=rep(prediction.3.psd, 7), x=plot.3.psd$lgn, col=cols[3], lty='dashed', lwd=1.7)
legend('topleft', legend=c("p = 0.1", "p = 0.5", "p = 0.9"), 
       fill=c(cols[1], cols[2], cols[3]), border='white', bty='n')

#adds figure label
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
