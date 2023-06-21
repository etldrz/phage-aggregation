
# file <- "C:\\Users\\Evan\\Desktop\\repos\\phageAggregation\\plottingData\\alpha and theta change. p=.3, lambda=.01, omega=.01\\base\\ alpha is 0.9 theta is 0.2 .txt"
# data <- read.table(file, header=T, sep=",")

library(foreach)
library(doParallel)

g <- data[,seq(2, ncol(data), 2)]
s <- data[,seq(1, ncol(data), 2)]



basicBootstrap <- function(g, s){
  bt_size <- 1e4

  size <- length(g)
  
  dfG <- as.data.frame(table(g))
  uniqueG <- as.numeric(levels(dfG[,1]))
  freqG <- dfG[,2] / size
  
  dfS <- as.data.frame(table(s))
  uniqueS <- as.numeric(levels(dfS[,1]))
  freqS <- dfS[,2] / size
  
  g_means <- replicate(bt_size, mean(sample(uniqueG, size, replace=T, prob=freqG)))
  s_means <- replicate(bt_size, mean(sample(uniqueS, size, replace=T, prob=freqS)))
  
  
  return(cbind(g_means, s_means))
}

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


start <- Sys.time()
currBoot <- foreach(i=1:2, .combine='cbind') %dopar% {
  basicBootstrap(g[,i], s[,i])
}
end <- Sys.time()

print(end - start)
stopCluster(cl)



