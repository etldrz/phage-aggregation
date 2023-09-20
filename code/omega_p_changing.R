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

fit.first.op <- preprocessed(first.op, omegas, "omega")
fit.second.op <- preprocessed(second.op, omegas, "omega")
fit.third.op <- preprocessed(third.op, omegas, "omega")

plot.1.op <- plotFitness(fit.first.op, FALSE)
plot.2.op <- plotFitness(fit.second.op, FALSE)
plot.3.op <- plotFitness(fit.third.op, FALSE)

total.op <- rbind(plot.1.op, plot.2.op, plot.3.op)
total.op <- cbind(total.op, p=rep(c(.1, .5, .9), each=5))
total.op <- cbind(total.op, prediction =
                    complexSimPrediction(alpha=alpha, theta=theta,
                                         p=rep(c(.1, .5, .9), each=5), lambda=lambda))

title.txt <- paste("alpha:", alpha,  "theta:", theta, "lambda:", lambda)

# 
# 
# plot.op <- ggplot(total.op, aes(x=omega, y=r.star.mean, group=factor(p),
#                                 color=factor(p))) +
#   geom_line(aes(x=omega, y=prediction), linetype='solid', alpha=0.9) +
#   geom_line(linetype='longdash') +
#   geom_point(size=1) +
#   geom_errorbar(aes(ymin=lower.quantile, ymax=upper.quantile), width=0.002,
#                 alpha=0.75) +
#   scale_fill_viridis(option="A") +
#   theme_classic() +
#   labs(color = "p") + 
#   ylim(0, 1) +
#   ggtitle(title.txt)
