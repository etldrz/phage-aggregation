library(doParallel)
library("flock")

# Organism types are generalist phage (G), specialist phage (S), good host (A), bad host (B)
# States are: G, S, A, B, G+A, G+B, S+A, RA, RB
#             1  2  3  4   5    6    7    8   9
stateNames = c("G", "S", "A", "B", "G+A", "G+B", "S+A", "RA", "RB")
source("/Users/draghi/My Drive/phage_metapopulations/change_matrices.R")

# Mean of Poisson-distributed number of virions produced (not decremented for loss of parent virion)
meanBurstA = 50

# Cell division consumes one unit of resource; these define the maximum yield per aggregations.
YA = 200
YB =  200

# Per-capita division rate is maxRate * R / (R + K)
maxRateA = 0.03
maxRateB = 0.03
KA = 40
KB = 40

# Resource decay rate--used to ensure that the simulation terminates by resource
# exhaustion even if all cells have been eliminated.
rd = 0.0005

# Per-capita burst rate of infected cells.
burstRate = 0.0005

# theta here is a rate constant for the interaction rate between phage and cells
theta = 1e-3

# Rate at which phage drift out of the aggregation
phageDrift = 1e-2

# Cells depart the aggregate actively based on a sigmoidal rate steepness * (1 - 1 / (1 + mid / R))
midA = 10
midB = 10
steepnessA = 0.1
steepnessB = 0.1

# The simulation ends when both resources are below this fraction of their initial abundances.
resourceThreshold = 0.05

# The rate at which a new phage particles arrives.
phageSettlingRate = 0.001

# The rate at which a bacterium arrives--this is independent of resource abundances.
hostASettlingRate = 0.2
hostBSettlingRate = 0.2

# Replicates to perform
reps = 5e5

expNumber = 8

mode = "together"
experiment = "burst_vs_settling"

if(experiment == "regular")
{
	alpha = 1e-2
	theta = 1e-3
	bs = seq(from = 0, to = 50, by = 5)
	bees = bs
}

if(experiment == "alpha_vs_theta")
{
	# Values of meanBurstB to use as treatments
	alphas = c(0.125e-2, 0.25e-2, 0.5e-2, 1e-2, 2e-2, 4e-2, 8e-2)
	thetas = c(1e-4, 1e-3, 1e-2)
	bs = seq(from = 0, to = 50, by = 5)
	alphaTreatments = rep(rep(alphas, length(thetas)), each = length(bs))
	thetaTreatments = rep(thetas, each = length(alphas) * length(bs))
	bees = rep(bs, length(thetas) * length(alphas))
}
if(experiment == "burst_vs_settling")
{
	# Values of meanBurstB to use as treatments
	settlings = c(0.000125, 0.00025, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.016, 0.032)
	bursts = c(0.0005, 0.0025)
	bs = seq(from = 0, to = 50, by = 5)
	settlingTreatments = rep(rep(settlings, length(bursts)), each = length(bs))
	burstTreatments = rep(bursts, each = length(settlings) * length(bs))
	bees = rep(bs, length(bursts) * length(settlings))
}

nTreatments = length(bees)

paraFile = paste0("/Users/draghi/My Drive/phage_metapopulations/no_pool/p_", formatC(expNumber, width=3, flag="0"), ".txt")
resFile = paste0("/Users/draghi/My Drive/phage_metapopulations/no_pool/res_", formatC(expNumber, width=3, flag="0"), ".txt")

write.table(cbind(c("meanBurstA", "YA", "YB", "maxRateA", "maxRateB", "KA", "KB", "rd", "burstRate", "theta", "phageDrift", "midA", "midB", "steepnessA", "steepnessB", "resourceThreshold", "phageSettlingRate", "hostASettlingRate", "hostBSettlingRate"), c(meanBurstA, YA, YB, maxRateA, maxRateB, KA, KB, rd, burstRate, theta, phageDrift, midA, midB, steepnessA, steepnessB, resourceThreshold, phageSettlingRate, hostASettlingRate, hostBSettlingRate)), paraFile, quote=FALSE, row.names=FALSE, col.names=FALSE)

seeds = sample((1+12e8):13e8, nTreatments, replace=FALSE)

hubQ = makeForkCluster(18)
registerDoParallel(hubQ)

temp = foreach(trial=1:nTreatments) %dopar%
{
  set.seed(seeds[trial])
  meanBurstB = bees[trial]
  if(experiment == "alpha_vs_theta")
  {
     theta = thetaTreatments[trial]
     phageDrift = alphaTreatments[trial]
  }
  if(experiment == "burst_vs_settling")
  {
  	phageSettlingRate = settlingTreatments[trial]
  	burstRate = burstTreatments[trial]
  }
  inputs = matrix(0, ncol=2, nrow=reps)
  outputs = matrix(0, ncol=2, nrow=reps)

  for(r in 1:reps)
  {
    # Initial state has full resources and no cells or phage
    if(mode == "together")
    {
      phageRates = c(phageSettlingRate, phageSettlingRate)
      d = c(0, 0, 0, 0, 0, 0, 0, YA, YB)
    } else if(mode == "apart") {

      phageRates = 2 * c(phageSettlingRate, phageSettlingRate)
      if(rbinom(1, 1, 0.5) == 1)
      {
        d = c(0, 0, 0, 0, 0, 0, 0, YA, 0)
      } else {
        d = c(0, 0, 0, 0, 0, 0, 0, 0, YB)
      }
    }
    insG = 0
    outsG = 0
    insS = 0
    outsS = 0
    
    while(d[8] >= resourceThreshold * YA | d[9] >= resourceThreshold * YB | sum(d[3:7]) > 5)
    {
      rates = rep(0, 21)
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
      rates[16:17] = phageRates
      # 18: A arrives
      # 19: B arrives
      hostRates = c(hostASettlingRate, hostBSettlingRate) * (d[8:9] > 0)
      rates[18:19] = hostRates
      # 20: RA decays by one
      # 21: RB decays by one
      rates[20:21] = rd * d[8:9]
      
      # Choose which event happens at that time step based on sampling each event
      # in proportion to its rate.
      event = sample(1:21, 1, prob=rates)
      
      # For simple changes (not involving Poisson-distributed outcomes), use this matrix to make changes.
      d = d + noPoolChanges[event,]
      
      # If an infected cell bursts, add its progeny to the simulation
      if(event == 3) d[1] = d[1] + rpois(1, meanBurstA)
      if(event == 4) d[1] = d[1] + rpois(1, meanBurstB)
      if(event == 5) d[2] = d[2] + rpois(1, meanBurstA)
      
      # Phage that leave the aggregate are counted as one "offspring" each
      if(event == 9) outsG = outsG + 1
      if(event == 10) outsS = outsS + 1
      
      # If an infected cell leaves, add its progeny to total reproduction
      if(event == 13) outsG = outsG + rpois(1, meanBurstA)
      if(event == 14) outsG = outsG + rpois(1, meanBurstB)
      if(event == 15) outsS = outsS + rpois(1, meanBurstA)
      
      # If a phage arrives, increment this counter for the denominator of fitness calculations
      if(event == 16) insG = insG + 1
      if(event == 17) insS = insS + 1
    }
    
    inputs[r,1] = insG
    inputs[r,2] = insS
    outputs[r,1] = outsG + d[1] + rpois(1, d[5] * meanBurstA + d[6] * meanBurstB)
    outputs[r,2] = outsS + d[2] + rpois(1, d[7] * meanBurstA)
  }
  
  x = colSums(outputs) / colSums(inputs)
  WG = x[1]
  WS = x[2]
  fl = lock(resFile)  
  if(experiment == "regular") write.table(cbind(meanBurstB, WG, WS, reps), resFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  if(experiment == "alpha_vs_theta") write.table(cbind(theta, phageDrift, meanBurstB, WG, WS, reps), resFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  if(experiment == "burst_vs_settling") write.table(cbind(burstRate, phageSettlingRate, meanBurstB, WG, WS, reps), resFile, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  unlock(fl)
} 
stopCluster(hubQ)
