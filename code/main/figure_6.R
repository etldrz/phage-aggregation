exps = 2
bLength = 31
wG = rep(0, bLength)
wS = wG

resFile = paste0("/Users/draghi/My Drive/phage_metapopulations/complex_model_outputs/res_", formatC(exps[i], width=3, flag="0"), ".txt")
dt = read.table(resFile, header=FALSE, col.names=c("b", "WG", "WS", "reps"))
bs = sort(unique(dt$b))
totalReps = tapply(dt$reps, dt$b, sum)
print(totalReps)
wG = tapply(dt$WG * dt$reps, dt$b, sum) / totalReps
wS = tapply(dt$WS * dt$reps, dt$b, sum) / totalReps

colors = c("darkorange3", "darkolivegreen3")
par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.1, 0.7, 0), cex.lab=1.4, cex.axis=1.1)
plot(range(bs), range(wG), type="n", xlab=expression(F[B]), ylab="Mean reproductive success")
points(bs, wS, type="l", lwd=2, col=colors[1])
points(bs, wG, type="l", lwd=2, col=colors[2])

legend("bottomright", legend=c("specialist", "generalist"), lwd=2, col=colors, cex=0.8)

