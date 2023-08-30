diver <- readRDS('diver.rds')

#diver covers alpha: 0.1, lambda: 0, omega: 0, p:.1
  #for rows 1:1e4 theta: 0.08. for rows 1e4+1:2e4 theta: 0.8

foo <- function(choice){
  lim <- c(1, 1e4)
  miny <- 4
  maxy <- 7
  if(choice=="upper"){
    miny <- 20
    maxy <- 30
    lim <- c(1e4 + 1, 2e4)
  }
  
  x <- seq(1, 52, by=50)
  plot(NA, xlim=c(mean(diver[lim[1]:lim[2],5] + 1) - 3, mean(diver[lim[1]:lim[2],5] + 1) + 3), 
       ylim=c(miny, maxy))
  
  for(i in lim[1]:lim[2]){
    c <- diver[i,]
    
    lines(x=x, y=(x*c[1] + c[2]), col='black', lwd=1)
    lines(x=x, y=(x*c[3] + c[4]), col='black', lwd=1)
  }
  
  tot.loc <- 1
  if(lim[1] > 1e4) tot.loc <- 2
  
  abline(v=c(min(diver[lim[1]:lim[2],5]), max(diver[lim[1]:lim[2],5])), col='red')
  abline(v=total.tp[tot.loc,7] + 1, col='black')
  
}

foo('upper')