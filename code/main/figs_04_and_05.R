library("minpack.lm")

# Figure 4

dr = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_003.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr2 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_020.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr = rbind(dr, dr2)
dAlt = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_008b.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
hits = print(tapply(dr$rep[dr$t == 1], dr$id[dr$t == 1], length))
altHits = print(tapply(dAlt$rep[dAlt$t == 1], dAlt$id[dAlt$t == 1], length))
weaker = c(13, 23, 50, 65, 72, 74, 80, 86, 94, 95, 97)
totalHits = hits
totalHits[weaker] = altHits
ts = sort(unique(dr$t))
threshold = ceiling(length(ts) * 1/10)
ids = sort(unique(dr$id))
nId = length(ids)
bees = sort(unique(dr$burstB))
midpoints = rep(0, nId)
par(mfrow=c(1,1))
for(id in ids)
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
  plot(c(0, 50), c(0,1), type="n", main = paste0(id, " ", length(xs)), xlab=expression(F[B]), ylab="Fraction specialists")
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
}

exclude = c(20, 68, 80, 86, 94, 95, 97)
goodIds = ids[-exclude]
midpoints = midpoints[goodIds]

dr = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_010.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dZero = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_011.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dExtra = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_012.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dr = rbind(dr, dZero, dExtra)
ts = sort(unique(dr$t))
threshold = ceiling(length(ts) * 1/10)
ids = sort(unique(dr$id))
nId = length(ids)
bees = sort(unique(dr$burstB))

dWeak = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_014.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeak2 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_015.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeak3 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_016.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeak4 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_018.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeak5 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_019.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeak6 = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_023.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dBatch = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_017.txt", header=FALSE, col.names=c("id", "rep", "burstB", "t", "G", "S", "A", "B", "opCost", "tries"))
dWeakNew = rbind(dWeak, dWeak2, dWeak3, dBatch, dWeak4, dWeak5, dWeak6)

controlWeaker = c(3, 8, 10, 11, 12, 14, 15, 18, 22, 23, 24, 26, 27, 29, 30, 32, 34, 40, 43, 45, 48, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 67, 68, 69, 70, 72, 74, 75, 77, 79, 81, 86, 89, 90, 91)
controlRegular = c(1,  2,  4,  5,  6,  7, 13, 16, 17, 19, 21, 25, 31, 33, 35, 36, 37, 38, 39, 44, 46, 47, 49, 64, 65, 66, 78, 82, 85, 93)
controlMidpoints = NULL
hits = NULL
for(id in ids)
{
  if(id %in% controlWeaker)
  {
    d = dWeakNew[dWeakNew$id == id,]
    col = "dodgerblue"
  } else {
    d = dr[dr$id == id,]
    col = "black"
  }
  if(id %in% controlWeaker) col = "dodgerblue"
  xs = NULL
  ys = NULL
  starts = which(d$t == 1)
  ends = which(d$t == max(d$t))
  for(i in 1:length(starts))
  {
    thisRep = d[(starts[i]):(ends[i]),]
    meanG = mean(thisRep$G[thisRep$t >= threshold])
    meanS = mean(thisRep$S[thisRep$t >= threshold])
    if(meanS > 0 | meanG > 0 & thisRep$A[5214] > 0 & thisRep$B[5214] > 0)
    {
      xs = c(xs, thisRep$burstB[thisRep$t == threshold])
      ys = c(ys, meanS / (meanS + meanG))
    }
  }
  hits = c(hits, length(ys))
  par(mfrow=c(1,1))
  plot(c(0, max(xs)), c(0,1), type="n", main=paste0(id, ", ", length(ys)))
  points(rep(1,2), c(0,1), type="l", lty="dotted")
  points(c(0, 12), rep(0.5, 2), type="l", lty="dotted")
  points(xs, ys, col=col)
  bs = seq(from = min(xs), to=max(xs), by = 0.01)
  if(any(ys > 0.5) & any(ys < 0.5))
  {
    estimates = rep(-10, 10)
    for(i in 1:100)
    {
      try(
        {
          mod = nlsLM(ys~d + f/(1 + exp(b * (xs-c))), start=list(b = 1, c=runif(1, 0, 4), d = 0.01, f = 1))
          slope = summary(mod)$coefficients[1]
          intercept = summary(mod)$coefficients[3]
          height = summary(mod)$coefficients[4]
          estimates[i] = summary(mod)$coefficients[2] + log(height / (0.5 - intercept) - 1) / slope
          points(bs, intercept + height / (1 + exp(slope * (bs - summary(mod)$coefficients[2]))), type="l")
          points(estimates[i], 0.5, pch=16)
        }, silent = TRUE)
    }
  }
  print(c(id, mean(estimates[which(estimates != -10)])))
  controlMidpoints = c(controlMidpoints, mean(estimates[which(estimates != -10)]))
}

# exclude those parameter sets for which the inferred midpoint is NaN or strongly negative
exclusions = c(8, 30, 40, 51, 54, 60, 61, 67, 77, 86)

controlMidpoints = controlMidpoints[-which(ids %in% exclusions)]

togetherRs = (midpoints - 1) / 49
apartRs = (controlMidpoints - 1) / 49
brks = seq(from=-0.075, to = 0.7, by=0.05)
localBlue = rgb(0.1, 0.52,0.96, 0.5)
localRed = rgb(0.9, 0.1, 0.1, 0.5)
histA = hist(apartRs, breaks=brks, col=localRed, xlab="R*", ylab="Counts", main="", cex.lab=1.2, cex.axis=0.9)
histB = hist(togetherRs, breaks=brks, add=TRUE, col=localBlue)
legend("topright", c("together", "apart"), col=c(localBlue, localRed), pch=15)


# Figure 5

do = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/op_cost_001.txt", header=FALSE, col.names=c("id", "rep", "tracked", "hits"))
do = do[-which(do$tracked < 100),]
ops = tapply(do$hits / do$tracked, do$id, mean)
print(tapply(do$hits / do$tracked, do$id, length))
xs = ops[goodIds]
ys = togetherRs

par(mfrow=c(1,1), mgp=c(1.9, 0.8, 0), mar=c(4,4,1,1))
plot(xs, ys, xlab="Opportunity cost P(A|B)", ylab="R*")
mod = lm(ys~xs)
b = summary(mod)$coefficients[1]
a = summary(mod)$coefficients[2]
deets = seq(from = min(xs), to = max(xs), by = 0.005)
points(deets, a * deets + b, type="l", lwd=2, col="darkgrey")


