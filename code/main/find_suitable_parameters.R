# Note: the "good_parameters_together" file is comprised of the results of four batches.
# res_001 and p_001 are used for the first batch of 11.
# res_004 and p_004 are used fro the second batch of 35
# res_006 and p_006 are used for the third batch of 35.
# res_007 and p_007 are used for the fourth batch of 22, of which 19 were used.

d = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/res_004.txt", header=FALSE, col.names=c("id", "rep", "meanBurstB", "t", "G", "S", "A", "B", "opCost", "reps"))
dp = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/p_004.txt", header=TRUE)

maxT = max(d$t)
minT = maxT / 2
ids = sort(unique(d$id))
minA = tapply(d$A[d$t >= minT], d$id[d$t >= minT], min)
minB = tapply(d$B[d$t >= minT], d$id[d$t >= minT], min)
minG = tapply(d$G[d$t >= minT], d$id[d$t >= minT], min)
meanA = tapply(d$A[d$t >= minT], d$id[d$t >= minT], mean)
meanB = tapply(d$B[d$t >= minT], d$id[d$t >= minT], mean)
meanG = tapply(d$G[d$t >= minT], d$id[d$t >= minT], mean)
completed = tapply(d$rep, d$id, max)
deltaA = rep(0, length(ids))
deltaB = rep(0, length(ids))
deltaG = rep(0, length(ids))

for(i in 1:length(ids))
{
  As = d$A[d$t >= minT & d$id == ids[i]]
  Bs = d$B[d$t >= minT & d$id == ids[i]]
  Gs = d$G[d$t >= minT & d$id == ids[i]]
  ts = d$t[d$t >= minT & d$id == ids[i]]
  modA = lm(As ~ ts)
  modB = lm(Bs ~ ts)
  modG = lm(Gs ~ ts)
  deltaA[i] = abs(modA$coefficients[2] / mean(As))
  deltaB[i] = abs(modB$coefficients[2] / mean(Bs))
  deltaG[i] = abs(modG$coefficients[2] / mean(Gs))
}

good = which(minA > 5000 & minB > 5000 & minG > 5000 & meanA > 10000 & meanB > 10000 & meanG > 10000 & deltaA < 5e-3 & deltaB < 5e-3 & deltaG < 5e-3)

write.table(cbind(1:length(good), dp$theta[good], dp$phageDrift[good], dp$cellDeathRate[good], dp$phageDeathRate[good], dp$cellMigrationCoefficient[good], dp$phageMigrationCoefficient[good], dp$patchPeriod[good], dp$maxRateA[good], dp$maxRateB[good]), "/Users/draghi/My Drive/phage_metapopulations/pool_returns/good_parameters_together_D.txt", quote=FALSE, col.names = c("id", "theta", "phageDrift", "cellDeathRate", "phageDeathRate", "cellMigrationCoefficient", "phageMigrationCoefficient", "patchPeriod", "maxRateA", "maxRateB"), row.names=FALSE)
