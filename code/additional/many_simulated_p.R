source('single_patch.R')

alpha <- 0.1
lambda <- 0.01
theta <- seq(0, 1, .05)
p <- seq(0, 1, 0.03)

container <- matrix(nrow=length(theta), ncol=length(p))

for(i in 1:nrow(container)){
  for(s in 1:ncol(container)){
    container[i,s] <- complexSimPrediction(alpha=alpha, theta=theta[i], p=p[s],
                                        lambda=lambda)
  }
}


plot.new()
by=seq(0, 1, 0.1)
axis(4, at=by); axis(1, at=by); axis(2, at=c(0, 0.5, 1))
for(i in 1:ncol(container)){
  lines(x=theta,y=container[,i])
}