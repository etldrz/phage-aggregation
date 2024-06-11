library(doParallel)
library("flock")
library("minpack.lm")

# Organism types are generalist phage (G), specialist phage (S), good host (A), bad host (B)
# States are: G, S, A, B, G+A, G+B, S+A, RA, RB
#             1  2  3  4   5    6    7    8   9
stateNames = c("G", "S", "A", "B", "G+A", "G+B", "S+A", "RA", "RB")
source("/Users/jeremydraghi/My Drive/phage_metapopulations/change_matrices.R")

# Mean of Poisson-distributed number of virions produced (not decremented for loss of parent virion)
meanBurstA = 50

# Cell division consumes one unit of resource; these define the maximum yield per aggregations.
YA = 200
YB =  200

# Per-capita division rate is maxRate * R / (R + K), where R is the relevant resource
# and K is the relevant half-saturation constant.
maxRateA = 0.03
maxRateB = 0.03
KA = 40
KB = 40

# Resource decay rate--used to ensure that the simulation terminates by resource
# exhaustion even if all cells have been eliminated.
rd = 0.001

# The simulation ends when both resources are below this fraction of their initial abundances.
resourceThreshold = 0.05

# Per-capita burst rate of infected cells.
burstRate = 0.005

# theta here is a rate constant for the infection-generating interaction 
# rate between phage and cells
theta = 1e-3

# Rate at which phage drift out of the aggregation
phageDrift = 1e-2

# Cells depart the aggregate actively based on a sigmoidal rate steepness * (1 - 1 / (1 + mid / R))
midA = 10
midB = 10
steepnessA = 0.1
steepnessB = 0.1

# Death rates
cellDeathRate = 2e-6
phageDeathRate = 2e-5

# Rates of entering a patch
cellMigrationCoefficient = 1e-5
phageMigrationCoefficient = 1.5e-6

# Update interval--controls how often patches are synchronized and pool death processed.
interval = 30

skips = as.integer(60/interval) * 24 * 7

initialA = 1e5
initialB = 1e5
initialPhage = 1e5

# Rate at which new patches occur. 
patchPeriod = 18

# Replicates to perform
reps = 25

# Migration strength modulated here
migrants = 2e-5
#migrants = 0.4e-5

expNumber = 21

# together == G and S with A+B patches
# apart == G and S with A or B patches
# pTester == G only, with variation in nuisance parameters
mode = "together"
pTesting = FALSE

if(pTesting == TRUE) maxT = 60 * 24 * 365 * 5
if(pTesting == FALSE) maxT = 60 * 24 * 365 * 100
burnIn = as.integer(60/interval) * 24 * 365 * 5

paraFile = paste0("/Users/jeremydraghi/My Drive/phage_metapopulations/pool_returns/p_", formatC(expNumber, width=3, flag="0"), ".txt")
resFile = paste0("/Users/jeremydraghi/My Drive/phage_metapopulations/pool_returns/res_", formatC(expNumber, width=3, flag="0"), ".txt")

if(pTesting == FALSE)
{
  paras = read.table("/Users/jeremydraghi/My Drive/phage_metapopulations/pool_returns/good_parameters_together.txt", header=TRUE)
  pSets = dim(paras)[1]
  nTreatments = pSets
}

if(pTesting == TRUE)
{
  nTreatments = 800
  cellDeathRates = exp(runif(nTreatments, min = log(1e-8), max=log(1e-5)))
  phageDeathRates = exp(runif(nTreatments, min = log(1e-5), max=log(1e-2)))
  cellMigrationCoefficients = exp(runif(nTreatments, min = log(1e-7), max=log(1e-4)))
  phageMigrationCoefficients = exp(runif(nTreatments, min = log(2e-8), max=log(2e-5)))
  phageDrifts = exp(runif(nTreatments, min = log(1e-3), max=log(1e-1)))
  thetas = exp(runif(nTreatments, min = log(1e-4), max=log(1e-2)))
  patchPeriods = sample(20:40, nTreatments, replace=TRUE)
  maxRatesA = exp(runif(nTreatments, min = log(0.005), max = log(0.1)))
  maxRatesB = exp(runif(nTreatments, min = log(0.005), max = log(0.1)))
  
  write.table(cbind("rep", "meanBurstA", "YA", "YB", "maxRateA", "maxRateB", "KA", "KB", "rd", "burstRate", "theta", "phageDrift", "midA", "midB", "steepnessA", "steepnessB", "resourceThreshold", "cellDeathRate", "phageDeathRate", "cellMigrationCoefficient", "phageMigrationCoefficient", "patchPeriod"), paraFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(cbind(1:nTreatments, rep(meanBurstA, nTreatments), rep(YA, nTreatments), rep(YB, nTreatments), rep(maxRatesA, nTreatments), rep(maxRatesB, nTreatments), rep(KA, nTreatments), rep(KB, nTreatments), rep(rd, nTreatments), rep(burstRate, nTreatments), thetas, phageDrifts, rep(midA, nTreatments), rep(midB, nTreatments), rep(steepnessA, nTreatments), rep(steepnessB, nTreatments), rep(resourceThreshold, nTreatments), cellDeathRates, phageDeathRates, cellMigrationCoefficients, phageMigrationCoefficients, patchPeriods), paraFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
}

seeds = sample((1+3e7):4e7, nTreatments, replace=FALSE)

hubO = makeForkCluster(6)
registerDoParallel(hubO)
temp = foreach(trial=rev(c(98, 92, 89, 79, 69, 54, 40, 35, 34, 15, 6, 11))) %dopar%
  {
    set.seed(seeds[trial])
    if(pTesting == FALSE)
    {
      cellDeathRate = paras$cellDeathRate[trial]
      phageDeathRate = paras$phageDeathRate[trial]
      cellMigrationCoefficient = paras$cellMigrationCoefficient[trial]
      phageMigrationCoefficient = paras$phageMigrationCoefficient[trial]
      phageDrift = paras$phageDrift[trial]
      theta = paras$theta[trial]
      patchPeriod = paras$patchPeriod[trial]
      maxRateA = paras$maxRateA[trial]
      maxRateB = paras$maxRateB[trial]
    }
    if(pTesting == TRUE)
    {
      cellDeathRate = cellDeathRates[trial]
      phageDeathRate = phageDeathRates[trial]
      cellMigrationCoefficient = cellMigrationCoefficients[trial]
      phageMigrationCoefficient = phageMigrationCoefficients[trial]
      phageDrift = phageDrifts[trial]
      theta = thetas[trial]
      patchPeriod = patchPeriods[trial]
      maxRateA = maxRatesA[trial]
      maxRateB = maxRatesB[trial]
    }
    
    for(r in 1:reps)
    {
      # G, S, A, B 
      pool = c(initialPhage, 0, initialA, initialB)
      patches = list()
      if(mode == "together")
      {
        patches[[1]] = c(0, 0, 0, 0, 0, 0, 0, YA, YB)
      } 
      if(mode == "apart")
      {
        patches[[1]] = c(0, 0, 0, 0, 0, 0, 0, YA, 0)
        patches[[2]] = c(0, 0, 0, 0, 0, 0, 0, 0, YB)
      }
      if(pTesting == TRUE)
      {
        meanBurstB = sample(0:25,1)
      } else {
        if(mode == "together")
        {
          meanBurstB = sample(0:40, 1)
        }
        if(mode == "apart")
        {
          meanBurstB = sample(seq(from = 0, to = 4, by = 0.5), 1)
        }
      }
      
      intervals = seq(from = 0, to = maxT, by = interval)
      cycles = as.integer(maxT / interval)
      res = NULL
      opCost = 0
      opCostReps = 0
      nextPatch = patchPeriod
      for(cycle in 1:cycles)
      {
        if(length(patches) > 0)
        {
          # Process patches in random order to minimize order effects
          aliases = sample(1:length(patches), length(patches))
          for(i in 1:length(patches))
          {
            d = patches[[aliases[i]]]
            localT = 0
            while(localT < interval)
            {
              rates = rep(0, 28)
              # 1: A reproduces (RA decreases by one)
              # 2: B reproduces (RB decreases by one)
              rates[1:2] = c(maxRateA / (KA / d[8] + 1) , maxRateB / (KB/d[9] + 1)) * d[3:4]
              # 3: G+A bursts
              # 4: G+B bursts
              # 5: S+A bursts
              rates[3:5] = burstRate * d[5:7]
              # 6: New G+A infection
              # 7: New G+B infection
              # 8: New S+A infection
              rates[6:8] = theta * d[c(1, 1, 2)] * d[c(3, 4, 3)]
              # 9: G leaves
              # 10: S leaves
              rates[9:10] = phageDrift * d[1:2]
              # 11: A leaves
              # 12: B leaves
              # 13: G+A leaves
              # 14: G+B leaves
              # 15: S+A leaves
              rates[c(11, 13, 15)] = steepnessA * (1 - 1 / (1 + midA / d[8])) * d[c(3, 5, 7)]
              rates[c(12, 14)] = steepnessB * (1 - 1 / (1 + midB / d[9])) * d[c(4, 6)]
              # 16: G arrives
              # 17: S arrives
              rates[16:17] = phageMigrationCoefficient * pool[1:2]
              # 18: A arrives
              # 19: B arrives
              rates[18:19] = cellMigrationCoefficient * pool[3:4] * (d[8:9] > 0)
              # 20: RA decays by one
              # 21: RB decays by one
              rates[20:21] = rd * d[8:9]
              # 22: A dies
              # 23: B dies
              # 24: G+A dies
              # 25: G+B dies
              # 26: S+A dies
              rates[22:26] = cellDeathRate * d[3:7]
              # 27: G dies
              # 28: S dies
              rates[27:28] = phageDeathRate * d[1:2]
              
              # Choose which event happens at that time step based on sampling each event
              # in proportion to its rate.
              summedRate = sum(rates)
              if(summedRate > 0)
              {
                event = sample(1:28, 1, prob=rates)
                deltaT = rexp(1, summedRate)
                localT = localT + deltaT
              } else {
                localT = interval
              }
              if(localT < interval)
              {
                # For simple changes (not involving Poisson-distributed outcomes), use this matrix to make changes.
                d = d + poolChanges[event,]
                stopifnot(all(d >= 0))
                # If an infected cell bursts, add its progeny to the simulation
                if(event == 3) d[1] = d[1] + rpois(1, meanBurstA)
                if(event == 4) d[1] = d[1] + rpois(1, meanBurstB)
                if(event == 5) d[2] = d[2] + rpois(1, meanBurstA)
                
                if(event == 7)
                {
                  opCostReps = opCostReps + 1
                  opCost = opCost + rates[6] / (rates[6] + rates[9] + rates[27])
                }
                
                if(event == 9) pool[1] = pool[1] + 1
                if(event == 10) pool[2] = pool[2] + 1
                if(event == 11) pool[3] = pool[3] + 1
                if(event == 12) pool[4] = pool[4] + 1
                
                # If an infected cell leaves, add its progeny to total reproduction
                if(event == 13) pool[1] = pool[1] + rpois(1, meanBurstA)
                if(event == 14) pool[1] = pool[1] + rpois(1, meanBurstB)
                if(event == 15) pool[2] = pool[2] + rpois(1, meanBurstA)
                if(event == 16) pool[1] = pool[1] - 1
                if(event == 17) pool[2] = pool[2] - 1
                if(event == 18) pool[3] = pool[3] - 1
                if(event == 19) pool[4] = pool[4] - 1
              }
            }
            patches[[aliases[i]]] = d
          }
        }
        # Delete exhausted patches. Use a while loop with dynamic limit to account for behavior of
        # deleting list elements, which slides elements past the deleted elements to the left. Hence,
        # after deleting an element we continue with the current index $i.
        i = 1
        while(i <= length(patches))
        {
          d = patches[[i]]
          if(d[8] < resourceThreshold * YA & d[9] < resourceThreshold * YB & sum(d[3:4]) <= 5)
          {
            # Before deleting a patch, process infections and empty all contents into pool
            pool[1] = pool[1] + d[1] + rpois(1, d[5] * meanBurstA) + rpois(1, d[6] * meanBurstB)
            pool[2] = pool[2] + d[2] + rpois(1, d[7] * meanBurstA)
            pool[3] = pool[3] + d[3]
            pool[4] = pool[4] + d[4]
            patches[[i]] = NULL
          } else {
            i = i + 1
          }
        }
        
        # Process pool deaths
        pool = rbinom(4, pool, exp(-c(phageDeathRate, phageDeathRate, cellDeathRate, cellDeathRate) * interval))
        if(pTesting == FALSE & cycle > burnIn)
        {
          phage = sum(pool[1:2])
          pool[1] = pool[1] + rpois(1,migrants * phage)
          pool[2] = pool[2] + rpois(1,migrants * phage)
        } 
        
        # Make new patches
        if(nextPatch == 0)
        {
          nextPatch = patchPeriod
          if(mode == "together")
          {
            patches[[length(patches)+1]] = c(0, 0, 0, 0, 0, 0, 0, YA, YB)
          } 
          if(mode == "apart")
          {
            patches[[length(patches)+1]] = c(0, 0, 0, 0, 0, 0, 0, YA, 0)
            patches[[length(patches)+1]] = c(0, 0, 0, 0, 0, 0, 0, 0, YB)
          }
        }
        
        if(cycle %% skips == 0)
        {
          # Census
          counts = rep(0, 4)
          if(length(patches) > 0)
          {
            for(i in 1:length(patches))
            {
              d = patches[[i]]
              counts = counts + d[1:4]
              counts[1] = counts[1] + d[5] * meanBurstA + d[6] * meanBurstB
              counts[2] = counts[2] + d[7] * meanBurstA
            }
          }
          counts = pool + counts
          res = rbind(res, c(counts, opCost, opCostReps))
          opCost = 0
          opCostReps = 0
        }
        nextPatch = nextPatch - 1
      }
      finalRes = res[dim(res)[1], 1:4]
      fl = lock(resFile)  
      write.table(cbind(rep(trial, cycles %/% skips), rep(r, cycles %/% skips), rep(meanBurstB, cycles %/% skips), 1:(cycles %/% skips), res), resFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
      unlock(fl)
      if(pTesting == TRUE & (finalRes[1] == 0 | finalRes[3] == 0 | finalRes[4] == 0)) break
    }
  } 
stopCluster(hubO)
