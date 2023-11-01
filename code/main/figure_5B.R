library(viridis)
expNumber = 3

resFile = paste0("/Users/draghi/My Drive/phage_metapopulations/complex_model_outputs/res_", formatC(expNumber, width=3, flag="0"), ".txt")
dt = read.table(resFile, header=FALSE, col.names=c("theta", "alpha", "b", "WG", "WS", "reps"))

alphas = sort(unique(dt$alpha))
thetas = sort(unique(dt$theta))
res = matrix(0, nrow=length(alphas), ncol=length(thetas))
par(mfrow=rep(ceiling(sqrt(length(alphas) * length(thetas))), 2), mar=c(1,2,1,1), mgp=c(1.9, 0.4, 0))
for(i in 1:length(thetas))
{
  for(j in 1:length(alphas))
  {
    theta = thetas[i]
    alpha = alphas[j]
    d = dt[which(dt$theta == theta & dt$alpha == alpha),]
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

cols = plasma(length(thetas)+1)[-4]
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.1, 0.7, 0), cex.lab=1.4, cex.axis=1.1)
plot(range(alphas), c(0,1), type="n", xlab=expression(paste("", alpha)), ylab=expression(paste("R"^"*")), log="x")
for(i in 1:length(thetas))
{
  points(alphas, (res[,i]-1) / 49, type="b", col=cols[i])
}
legend("topright", legend=c(expression(paste(theta, " = ", 1e-4)), expression(paste(theta, " = ", 1e-3)), expression(paste(theta, " = ", 1e-2))), lwd=1.5, col=cols)


