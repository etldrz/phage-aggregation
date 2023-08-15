source('single_patch.R')

library(foreach)
library(ggplot2)
library(wesanderson)

base.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/base/"
boot.file <- "C:/Users/Evan/Desktop/repos/phage-aggregation/data/bootstrapped/"

alpha <- 0.1; thetas <- c(0.08, 0.8); ps <- c(0.1, 0.5, 0.9); lambda <- 0.01; omega <- 0

# boot.files.tp <- c()
# base.files.tp <- c()
# 
# for(theta in thetas) {
#   for(p in ps){
#     locat <- paste("alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
#                    lambda, ", omega=", omega, ".txt", sep="")
# 
#     curr.base <- paste(base.file, locat, sep="")
#     curr.boot <- paste(boot.file, locat, sep="")
# 
#     base.files.tp <- c(base.files.tp, curr.base)
#     boot.files.tp <- c(boot.files.tp, curr.boot)
#     
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
#   }
# }

first.tp <- boot.files.tp[c(1, 4)]
second.tp <- boot.files.tp[c(2, 5)]
third.tp <- boot.files.tp[c(3, 6)]

fit.first.tp <- preprocessed(first.tp, thetas, "theta")
fit.second.tp <- preprocessed(second.tp, thetas, "theta")
fit.third.tp <- preprocessed(third.tp, thetas, "theta")

plot.1.tp <- plotFitness(fit.first.tp, FALSE)
plot.2.tp <- plotFitness(fit.second.tp, FALSE)
plot.3.tp <- plotFitness(fit.third.tp, FALSE)

total.tp <- rbind(plot.1.tp, plot.2.tp, plot.3.tp)
total.tp <- cbind(total.tp, p=c(.1, .1, .5, .5, .9, .9))
total.tp <- cbind(total.tp, prediction =
                    complexSimPrediction(alpha, rep(c(0.08, 0.8), 3),
                                         c(.1, .1, .5, .5, .9, .9), lambda))
total.tp <- cbind(total.tp, predicted.nB = 
                    total.tp[,6] * (burst.size.A - 1))

data <- c()
data <- c(data, (total.tp[1,7] - mean(fit.first.tp[[1]][[14]])))
data <- c(data, (total.tp[2,7] - mean(fit.first.tp[[2]][[14]])))
data <- c(data, (total.tp[3,7] - mean(fit.second.tp[[1]][[14]])))
data <- c(data, (total.tp[4,7] - mean(fit.second.tp[[2]][[14]])))
data <- c(data, (total.tp[5,7] - mean(fit.third.tp[[1]][[14]])))
data <- c(data, (total.tp[6,7] - mean(fit.third.tp[[2]][[14]])))

total.tp <- cbind(total.tp, difference = data)

title.txt <- paste("alpha:", alpha, "lambda:", lambda, "omega:", omega)

plot.tp <- ggplot(total.tp, aes(x=theta, y=r.star.mean, group=factor(p),
                                color=factor(p))) +
  geom_line(aes(x=theta, y=prediction), linetype='longdash', alpha=0.9) +
  geom_line() +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=lower.quantile, ymax=upper.quantile), width=0.002,
                alpha=0.75) +
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1")) +
  theme_classic() +
  labs(color = "p") +
  ggtitle(title.txt)

total.tp
