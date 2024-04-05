# library(ggplot2)
# library(foreach)
# setwd("~/Desktop/repos/phage-aggregation/code/main")
# source("single_patch.R")
# 
# base.file <- "~/Desktop/p.a_revisions/varying_p/base/"
# boot.file <- "~/Desktop/p.a_revisions/varying_p/boot/"
# 
# omega <- 0.05; lambda <- 0.01; alpha <- 0.1; theta <- 0.35; ps <- c(.1, .5, .9)
# n <- c(5, 10, 20, 30, 50, 500, 5000)
# boot.files.psd <- base.files.psd <- c()
# 
# for(p in ps) {
#   for(popsize in n) {
#     locat <- paste("n=", popsize, ", alpha=", alpha, ", theta=", theta, ", p=", p, ", lambda=",
#                    lambda, ", omega=", omega, ".txt", sep="")
# 
#     curr.base <- paste(base.file, locat, sep="")
#     curr.boot <- paste(boot.file, locat, sep="")
# 
#     base.files.psd <- c(base.files.psd, curr.base)
#     boot.files.psd <- c(boot.files.psd, curr.boot)
# # 
# #     file.create(curr.base)
# #     file.create(curr.boot)
# # 
# #     curr.base.data <- baseSimulation(alpha, theta, p, lambda, omega, n=popsize)
# #     print(dim(curr.base.data))
# # 
# #     write.table(curr.base.data, curr.base, quote=FALSE, sep=",", append=FALSE)
# # 
# #     s <- curr.base.data[,seq(1, ncol(curr.base.data), 2)]
# #     g <- curr.base.data[,seq(2, ncol(curr.base.data), 2)]
# # 
# #     curr.boot.data <- foreach(i = 1:ncol(g), .combine='cbind') %do% {
# #       basicBootstrap(s[,i], g[,i])
# #     }
# #     print(dim(curr.boot.data))
# # 
# #     write.table(curr.boot.data, file=curr.boot, append=FALSE, quote=FALSE, sep=",")
#   }
# }

# p.1 <- boot.files.psd[1:7]
# p.5 <- boot.files.psd[8:14]
# p.9 <- boot.files.psd[15:21]
# 
# fit.1 <- preprocessed(p.1, n, "n")
# fit.5 <- preprocessed(p.5, n, "n")
# fit.9 <- preprocessed(p.9, n, "n")
#
plt.1 <- cbind(plotFitness(fit.1, F), p=ps[1])
plt.5 <- cbind(plotFitness(fit.5, F), p=ps[2])
plt.9 <- cbind(plotFitness(fit.9, F), p=ps[3])

sdp1 <- c()
sdp5 <- c()
sdp9 <- c()
for(i in 1:nrow(plt.1)) {
  sdp1 <- c(sdp1, sd(rbinom(reps, n[i], plt.1$p[i]) / n[i]))
  sdp5 <- c(sdp5, sd(rbinom(reps, n[i], plt.5$p[i]) / n[i]))
  sdp9 <- c(sdp9, sd(rbinom(reps, n[i], plt.9$p[i]) / n[i]))
}
plt.1$sdp <- sdp1
plt.5$sdp <- sdp5
plt.9$sdp <- sdp9


plt.1$sdp <- rbinom(reps, n, plt.1$p) / n
plt.5$n <- exp(plt.5$n)
plt.9$n <- exp(plt.9$n)
cols <- c("0.1"="firebrick", "0.5"="skyblue",
          "0.9"="forestgreen")

plt.1$lgn <- log(plt.1$n)
plt.5$lgn <- log(plt.5$n)
plt.9$lgn <- log(plt.9$n)

ggplot() +
  geom_point(data=plt.1, mapping=aes(x=lgn, y=plt.1[,1], color=cols[1])) +
  geom_line(data=plt.1, mapping=aes(x=lgn, y=plt.1[,1], color=cols[1])) +
  geom_errorbar(data=plt.1, mapping=aes(x=lgn, y=plt.1[,1],
                                        ymin=plt.1[,2], ymax=plt.1[,3], color=cols[1]),
                width=0.001) +
  geom_point(data=plt.5, mapping=aes(x=lgn, y=plt.5[,1], color=cols[2])) +
  geom_line(data=plt.5, mapping=aes(x=lgn, y=plt.5[,1], color=cols[2])) +
  geom_errorbar(data=plt.5, mapping=aes(x=lgn, y=plt.5[,1],
                                        ymin=plt.5[,2], ymax=plt.5[,3], color=cols[2]),
                width=0.001) +
  geom_point(data=plt.9, mapping=aes(x=lgn, y=plt.9[,1], color=cols[3])) +
  geom_line(data=plt.9, mapping=aes(x=lgn, y=plt.9[,1], color=cols[3])) +
  geom_errorbar(data=plt.9, mapping=aes(x=lgn, y=plt.9[,1],
                                        ymin=plt.9[,2], ymax=plt.9[,3], color=cols[3]),
                width=0.001) +
  theme_classic() +
  labs(x = "log(population size)",
       y = "R*") +
  scale_color_discrete(name="p=",
                      breaks=c("firebrick", "skyblue",
                               "forestgreen"),
                      labels=c("0.1", "0.5", "0.9"))

ggplot() +
  geom_point(data=plt.1, mapping=aes(x=sdp, y=r.star.mean, color=cols[1])) +
  geom_line(data=plt.1, mapping=aes(x=sdp, y=r.star.mean, color=cols[1])) +
  geom_errorbar(data=plt.1, mapping=aes(x=sdp, y=r.star.mean,
                                        ymin=plt.1[,2], ymax=plt.1[,3], color=cols[1]),
                width=0.001) +
  geom_point(data=plt.5, mapping=aes(x=sdp, y=r.star.mean, color=cols[2])) +
  geom_line(data=plt.5, mapping=aes(x=sdp, y=r.star.mean, color=cols[2])) +
  geom_errorbar(data=plt.5, mapping=aes(x=sdp, y=r.star.mean,
                                        ymin=plt.5[,2], ymax=plt.5[,3], color=cols[2]),
                width=0.001) +
  geom_point(data=plt.9, mapping=aes(x=sdp, y=r.star.mean, color=cols[3])) +
  geom_line(data=plt.9, mapping=aes(x=sdp, y=r.star.mean, color=cols[3])) +
  geom_errorbar(data=plt.9, mapping=aes(x=sdp, y=r.star.mean,
                                        ymin=plt.9[,2], ymax=plt.9[,3], color=cols[3]),
                width=0.001) +
  theme_classic() +
  labs(x = "sd(p)",
       y = "R*") +
  scale_color_discrete(name="p=",
                       breaks=c("firebrick", "skyblue",
                                "forestgreen"),
                       labels=c("0.1", "0.5", "0.9"))



# 
# sdp <- function(x) {
#   sqrt(x * (1-x) / n)
# }
# 
# sd.df <- data.frame(n=rep(n, 3))
# sd.df$p.1 <- rep(ps, 3)
# sd.df$sd.p <- sqrt(sd.df$p*(1-sd.df$p)/sd.df$n)
# 
# plot(x=n, y=sdp(0.1), type="l", ylab="sd(p)", col="gold", lwd=3)
# lines(x=n, y=sdp(0.5), col="purple", lwd=3)
# legend("topright", legend=c("p=0.1", "p=0.5"), fill=c("gold", "purple"))
