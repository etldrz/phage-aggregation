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
theta = 1e-2

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
phageSettlingRate = 0
phageRates = c(phageSettlingRate,phageSettlingRate)

# The rate at which a bacterium arrives--this is independent of resource abundances.
hostASettlingRate = 0.2
hostBSettlingRate = 0.2

# Replicates to perform
reps = 4

increment = 1

# Values of meanBurstB to use as treatments
meanBurstB = 20

par(mfrow=c(4,1), mar=c(1.5,3,1,1), mgp=c(2.5, 0.6, 0))
for(r in 1:reps)
{
  # Initial state has full resources and no cells or phage
  d = c(0, 0, 0, 0, 0, 0, 0, YA, YB)
  record = NULL
  cumulativePhageG = NULL
  cumulativePhageS = NULL
  ts = NULL
  insG = 0
  outsG = 0
  insS = 0
  outsS = 0
  t = 0
  nextT = increment
  
  while(d[8] >= resourceThreshold * YA | d[9] >= resourceThreshold * YB | sum(d[3:7]) > 5)
  {
    if(t >= nextT)
    {
      ts = c(ts, t)
      record = rbind(record, d)
      cumulativePhageG = c(cumulativePhageG, outsG)
      cumulativePhageS = c(cumulativePhageS, outsS)
      nextT = nextT + increment
    }
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
    rates[18:19] = c(hostASettlingRate, hostBSettlingRate)
    # 20: RA decays by one
    # 21: RB decays by one
    rates[20:21] = rd * d[8:9]
    
    # Choose time step based on the exponential distribution for summed rates
    deltaT = rexp(1, sum(rates))
    t = t + deltaT
    
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
  # States are: G, S, A, B, G+A, G+B, S+A, RA, RB
  #             1  2  3  4   5    6    7    8   9
  plot(c(0, 5), c(0, 210), type="n", ylab="", xlab="", cex.axis=1.2)
  points(ts / 60, record[,3] + record[,5] + record[,7], type="l", col="deepskyblue4", lwd=1.5)
  points(ts / 60, record[,4] + record[,6], type="l", col="darkorange4", lwd=1.5)
  #points(ts, cumulativePhageG, type="l", col="firebrick2")
  #points(ts, cumulativePhageS, type="l", col="firebrick2", lty="dotted")
  points(ts / 60, record[,8], type="l", lwd=2.5, col="deepskyblue")
  points(ts / 60, record[,9], type="l", lwd=2.5, col="darkorange")
}
 