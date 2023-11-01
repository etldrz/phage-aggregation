library(viridis)
expNumber = 4

resFile = paste0("/Users/draghi/My Drive/phage_metapopulations/complex_model_outputs/res_", formatC(expNumber, width=3, flag="0"), ".txt")
dt = read.table(resFile, header=FALSE, col.names=c("burst", "arrival", "b", "WG", "WS", "reps"))

bursts = sort(unique(dt$burst))
arrivals = sort(unique(dt$arrival))
res = matrix(0, nrow=length(arrivals), ncol=length(bursts))
par(mfrow=rep(ceiling(sqrt(length(bursts) * length(arrivals))), 2), mar=c(1,2,1,1), mgp=c(1.9, 0.4, 0))
for(i in 1:length(bursts))
{
  for(j in 1:length(arrivals))
  {
    burst = bursts[i]
    arrival = arrivals[j]
    d = dt[which(dt$burst == burst & dt$arrival == arrival),]
    totalReps = tapply(d$reps, d$b, sum)
    print(totalReps)
    wg = tapply(d$WG * d$reps, d$b, sum) / totalReps
    ws = tapply(d$WS * d$reps, d$b, sum) / totalReps
    bs = sort(unique(d$b))
    plot(bs, ws, col="dodgerblue", ylim=range(c(ws, wg)))
    points(bs, wg, col="orange")
    modS = lm(ws ~ bs)
    bsSq = bs^2
    modG = lm(wg ~ bs + bsSq)
    sC = modS$coefficients
    gC = modG$coefficients
    
    intercept = (sC[2] - gC[2] + sqrt((gC[2] - sC[2])^2  - 4 * gC[3] * (gC[1] - sC[1]))) / (2 * gC[3])
    points(bs, sC[1] + sC[2] * bs, type="l")
    points(bs, gC[1] + gC[2] * bs + gC[3] * bs^2, type="l")
    points(intercept, sC[1] + sC[2] * intercept, col="red", pch=16)
    res[j,i] = intercept
  }
}

cols = plasma(length(bursts)+1)[-3]
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.1, 0.7, 0), cex.lab=1.4, cex.axis=1.1)
plot(range(arrivals), c(0,1), type="n", xlab=expression(paste("", phi)), ylab=expression(paste("R"^"*")), log="x")
for(i in 1:length(bursts))
{
  points(arrivals, (res[,i]-1) / 49, type="b", col=cols[i])
}
legend("topright", legend=c(expression(paste(beta, " = ", 0.0005)), expression(paste(beta, " = ", 0.0025))), lwd=1.5, col=cols)


