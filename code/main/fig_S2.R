library("minpack.lm")

dr = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_003.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr2 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_020.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr = rbind(dr, dr2)
dAlt = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_008b.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr30 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_021.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr15 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_022.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
hits = print(tapply(dr$rep[dr$t == 1], dr$id[dr$t == 1], length))
altHits = print(tapply(dAlt$rep[dAlt$t == 1], dAlt$id[dAlt$t == 1], length))
weaker = c(13, 23, 50, 65, 72, 74, 80, 86, 94, 95, 97)
picks = sort(c(98, 92, 89, 79, 69, 54, 40, 35, 34, 15, 6, 11))
totalHits = hits
totalHits[weaker] = altHits
ts = sort(unique(dr$t))
threshold = ceiling(length(ts) * 1/10)
ids = sort(unique(dr$id))
nId = length(ids)
bees = sort(unique(dr$burstB))
midpoints = rep(0, nId)
transparentRed = rgb(255,0,0,150, maxColorValue = 255)
transparentPurple = rgb(191,62,255,150, maxColorValue = 255)
par(mfrow=c(4,3), mar=c(2.5,2.5,1,1), mgp=c(1.25, 0.4, 0), cex.axis=0.95)
for(id in picks)
{
  if(id %in% weaker)
  {
    d = dAlt[dAlt$id == id,]
    col = "dodgerblue"
  } else {
    d = dr[dr$id == id,]
    col = "black"
  }
  
  xs = NULL
  ys = NULL
  endB = NULL
  starts = which(d$t == 1)
  ends = which(d$t == max(d$t))
  for(i in 1:length(starts))
  {
    thisRep = d[(starts[i]):(ends[i]),]
    meanG = mean(thisRep$G[thisRep$t >= threshold])
    meanS = mean(thisRep$S[thisRep$t >= threshold])
    if((meanS > 0 | meanG > 0) & thisRep$A[5214] > 0 & thisRep$B[5214] > 0)
    {
      xs = c(xs, thisRep$burstB[thisRep$t == threshold])
      ys = c(ys, meanS / (meanS + meanG))
    }
  }
  plot(c(0, 50), c(0,1), type="n", xlab=expression(F[B]), ylab="Fraction specialists")
  points(rep(1,2), c(0,1), type="l", lty="dotted")
  points(c(0, 55), rep(0.5, 2), type="l", lty="dotted")
  points(xs, ys, col=col)
  bs = seq(from = min(xs), to=max(xs), by = 0.01)
  if(any(ys > 0.5) & any(ys < 0.5))
  {
    estimates = rep(-10, 10)
    for(i in 1:10)
    {
      try(
        {
          mod = nlsLM(ys~d + f/(1 + exp(b * (xs-c))), start=list(b = 1, c=runif(1, 2, 20), d = 0.01, f = 1))
          slope = summary(mod)$coefficients[1]
          intercept = summary(mod)$coefficients[3]
          height = summary(mod)$coefficients[4]
          estimates[i] = summary(mod)$coefficients[2] + log(height / (0.5 - intercept) - 1) / slope
          points(bs, intercept + height / (1 + exp(slope * (bs - summary(mod)$coefficients[2]))), type="l")
          points(estimates[i], 0.5, pch=16)
        }, silent = TRUE)
    }
  }
  midpoints[id] = mean(estimates[which(estimates != -10)])
  
  d = dr30[dr30$id == id,]
  if(dim(d)[1] > 0)
  {
    xs = NULL
    ys = NULL
    endB = NULL
    starts = which(d$t == 1)
    ends = which(d$t == max(d$t))
    for(i in 1:length(starts))
    {
      thisRep = d[(starts[i]):(ends[i]),]
      meanG = mean(thisRep$G[thisRep$t >= threshold])
      meanS = mean(thisRep$S[thisRep$t >= threshold])
      if((meanS > 0 | meanG > 0) & thisRep$A[5214] > 0 & thisRep$B[5214] > 0)
      {
        xs = c(xs, thisRep$burstB[thisRep$t == threshold])
        ys = c(ys, meanS / (meanS + meanG))
      }
    }
    points(xs, ys, col=transparentRed, pch=3, cex=1.1)
  }
  
  d = dr15[dr15$id == id,]
  if(dim(d)[1] > 0)
  {
    xs = NULL
    ys = NULL
    endB = NULL
    starts = which(d$t == 1)
    ends = which(d$t == max(d$t))
    for(i in 1:length(starts))
    {
      thisRep = d[(starts[i]):(ends[i]),]
      meanG = mean(thisRep$G[thisRep$t >= threshold])
      meanS = mean(thisRep$S[thisRep$t >= threshold])
      if((meanS > 0 | meanG > 0) & thisRep$A[5214] > 0 & thisRep$B[5214] > 0)
      {
        xs = c(xs, thisRep$burstB[thisRep$t == threshold])
        ys = c(ys, meanS / (meanS + meanG))
      }
    }
    points(xs, ys, col=transparentPurple, pch=4, cex=1.1)
  }
}

