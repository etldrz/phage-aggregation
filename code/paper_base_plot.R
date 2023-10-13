source('single_patch.R')

burst.size.A <- 50
burst.sizes.B <- 0:(burst.size.A + as.integer(burst.size.A * 0.05))


# paper.base.plot.base <- "C:\\Users\\Evan\\Desktop\\repos\\phage-aggregation\\data\\base\\alpha=0.1, theta=0.2, p=0.5, lambda=0.01, omega=0.txt"
# 
# data <- read.table(paper.base.plot.base, header=TRUE, sep=",")

#r.star <- rStar(paper.base.plot.base, 0.2)

nB <- complexSimPrediction(alpha=0.1, theta=0.2, p=0.5, lambda=0.01)

burst.B <- (nB * (burst.size.A - 1)) + 1

s <- colMeans(data[,seq(1, ncol(data), 2)])
g <- colMeans(data[,seq(2, ncol(data), 2)])

min <- min(g)
max <- max(g)


plot(x=burst.sizes.B, y=s, col='firebrick',  type='p', cex=.75, ylim=c(min, max), pch=16, 
     xlab="Burst size of host B", ylab="Fitness", cex.lab=1.35)
points(x=burst.sizes.B, y=g, col='darkorchid3', cex=.75, pch=16)
abline(v=burst.B, lty='dotted', col='black', lwd=1.5)

legend('topleft', legend=c("Specialist", "Generalist"), 
       col=c('firebrick', 'darkorchid3'), border='white', pch=16, bty='n', cex=1.25)



