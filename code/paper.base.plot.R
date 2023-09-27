burst.size.A <- 50
burst.sizes.B <- 0:(burst.size.A + as.integer(burst.size.A * 0.05))


# paper.base.plot <-  
#"C:\\Users\\Evan\\Desktop\\repos\\phage-aggregation\\data\\
#bootstrapped\\alpha=0.1, theta=0.2, p=0.5, lambda=0.01, omega=0.txt"
# 
# data <- read.table(paper.base.plot, header=TRUE, sep=",")



s <- colMeans(data[,seq(1, ncol(data), 2)])
g <- colMeans(data[,seq(2, ncol(data), 2)])

min <- min(g)
max <- max(g)


plot(x=burst.sizes.B, y=s, col='firebrick',  type='p', cex=.75, ylim=c(min, max), pch=16, 
     xlab="Burst size of host B", ylab="Fitness")
points(x=burst.sizes.B, y=g, col='darkorchid3', cex=.75, pch=16)
legend('topleft', legend=c("Specialist", "Generalist"), 
       col=c('firebrick', 'darkorchid3'), pch=16)