library(doParallel)
library("flock")

source("/Users/draghi/My Drive/phage_metapopulations/change_matrices.R")

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
interval = 20

skips = 3 * 24 * 7

initialA = 1e5
initialB = 1e5
initialPhage = 1e5

# Rate at which new patches occur. 
patchPeriod = 18

# Replicates to perform
reps = 11

migrants = 2e-5

expNumber = 1
maxT = interval * 3 * 24 * 365 * 20
burnIn = 3 * 24 * 365 * 5
limit = 3 * 24 * 365 * 15

resFile = paste0("/Users/draghi/My Drive/phage_metapopulations/pool_returns/op_cost_", formatC(expNumber, width=3, flag="0"), ".txt")

paras = read.table("/Users/draghi/My Drive/phage_metapopulations/pool_returns/good_parameters_together.txt", header=TRUE)
pSets = dim(paras)[1]
nTreatments = pSets

seeds = sample((1+5e8):6e8, nTreatments, replace=FALSE)

hubW = makeForkCluster(20)
registerDoParallel(hubW)

temp = foreach(trial=c(1:100)) %dopar%
{
  set.seed(seeds[trial])

  cellDeathRate = paras$cellDeathRate[trial]
  phageDeathRate = paras$phageDeathRate[trial]
  cellMigrationCoefficient = paras$cellMigrationCoefficient[trial]
  phageMigrationCoefficient = paras$phageMigrationCoefficient[trial]
  phageDrift = paras$phageDrift[trial]
  theta = paras$theta[trial]
  patchPeriod = paras$patchPeriod[trial]
  maxRateA = paras$maxRateA[trial]
  maxRateB = paras$maxRateB[trial]

  for(r in 1:reps)
  {
    # S, S_ghost, A, B 
    pool = c(initialPhage, 0, initialA, initialB)
    patches = list()
    # S, S_ghost, A, B, S + A, YA, YB
    patches[[1]] = c(0, 0, 0, 0, 0, YA, YB)
   
    intervals = seq(from = 0, to = maxT, by = interval)
    cycles = as.integer(maxT / interval)
    tracked = 0
    fits = 0
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
            rates = rep(0, 22)
            # 1: A reproduces (RA decreases by one)
            # 2: B reproduces (RB decreases by one)
            rates[1:2] = c(maxRateA / (KA / d[6] + 1) , maxRateB / (KB/d[7] + 1)) * d[3:4]
            # 3: S+A bursts
            rates[3] = burstRate * d[5]
            # 4: New S+A infection
            # 5: New S+B ghost infection
            # 6: New S_ghost+A infection
            rates[4:6] = theta * d[c(1, 1, 2)] * d[c(3, 4, 3)]
            if(cycle < burnIn | cycle >= limit) rates[5] = 0
            # 7: S leaves
            # 8: S_ghost leaves
            rates[7:8] = phageDrift * d[1:2]
            # 9: A leaves
            # 10: B leaves
            # 11: S+A leaves
            rates[c(9, 11)] = steepnessA * (1 - 1 / (1 + midA / d[6])) * d[c(3, 5)]
            rates[10] = steepnessB * (1 - 1 / (1 + midB / d[7])) * d[4]
            # 12: S arrives
            # 13: S_ghost arrives
            rates[12:13] = phageMigrationCoefficient * pool[1:2]
            # 14: A arrives
            # 15: B arrives
            rates[14:15] = cellMigrationCoefficient * pool[3:4] * (d[6:7] > 0)
            # 16: RA decays by one
            # 17: RB decays by one
            rates[16:17] = rd * d[6:7]
            # 18: A dies
            # 19: B dies
            # 20: S+A dies
            rates[18:20] = cellDeathRate * d[3:5]
            # 21: S dies
            # 22: S_ghost dies
            rates[21:22] = phageDeathRate * d[1:2]
            
            # Choose which event happens at that time step based on sampling each event
            # in proportion to its rate.
            summedRate = sum(rates)
            if(summedRate > 0)
            {
            	event = sample(1:22, 1, prob=rates)
            	deltaT = rexp(1, summedRate)
            	localT = localT + deltaT
            } else {
            	localT = interval
            }
            if(localT < interval)
            {
              # For simple changes (not involving Poisson-distributed outcomes), use this matrix to make changes.
              d = d + specialistChanges[event,]
              stopifnot(all(d >= 0))
              # If an infected cell bursts, add its progeny to the simulation
              if(event == 3) d[1] = d[1] + rpois(1, meanBurstA)
              if(event == 5) tracked = tracked + 1
              if(event == 6) fits = fits + 1
              
              if(event == 7) pool[1] = pool[1] + 1
              if(event == 8) pool[2] = pool[2] + 1
              if(event == 9) pool[3] = pool[3] + 1
              if(event == 10) pool[4] = pool[4] + 1
              
              # If an infected cell leaves, add its progeny to total reproduction
              if(event == 11) pool[1] = pool[1] + rpois(1, meanBurstA)
              
              if(event == 12) pool[1] = pool[1] - 1
              if(event == 13) pool[2] = pool[2] - 1
              if(event == 14) pool[3] = pool[3] - 1
              if(event == 15) pool[4] = pool[4] - 1
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
        if(d[6] < resourceThreshold * YA & d[7] < resourceThreshold * YB & sum(d[3:4]) <= 5)
        {
          # Before deleting a patch, process infections and empty all contents into pool
          pool[1] = pool[1] + d[1] + rpois(1, d[5] * meanBurstA)
          pool[2] = pool[2] + d[2]
          pool[3] = pool[3] + d[3]
          pool[4] = pool[4] + d[4]
          patches[[i]] = NULL
        } else {
          i = i + 1
        }
      }

      # Process pool deaths
      pool = rbinom(4, pool, exp(-c(phageDeathRate, phageDeathRate, cellDeathRate, cellDeathRate) * interval))
            
      # Make new patches
      if(nextPatch == 0)
      {
      	nextPatch = patchPeriod
    		patches[[length(patches)+1]] = c(0, 0, 0, 0, 0, YA, YB)
      }
      nextPatch = nextPatch - 1
    }
    
    fl = lock(resFile)  
    write.table(cbind(trial, r, tracked, fits), resFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    unlock(fl)
  }
} 
stopCluster(hubW)
