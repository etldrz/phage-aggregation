
prop_allowed <- 0.01


#' Helper function used throughout
#' For verifying floating numbers
floatMatch <- function(x, y) {
  return(abs(x - y) < 1e-6)
}



#' This attempts to swap overlapping locations with non-overlapping locations 
#' from the matrix so that there is no overlap at all. The first column is 
#' expected to be larger than the second; locations where the second is larger
#' than the first is considered to be overlap.
swap <- function(means){
  
  # Vector of locations inside means where the second column is > than the first
  to_fix <- which(means[,1] <= means[,2])
  left_to_fix <- length(to_fix)

  # This is true when perc == 0
  prop.allowed <- FALSE
  
  for(index in to_fix){
    # The current means location trying to be fixed.
    overlap <- means[index,]
    
    # Locations that could potentially swap with the current overlap row
    viable <- which(!1:nrow(means) %in% to_fix)
    
    # Records how many times a redraw happens (how many items from viable were 
    # used) Not currently being used anywhere.
    count <- 1 
    
    while(length(viable) > 0){
      
      rand <- sample(viable, 1)
      
      curr_viable <- means[rand,]
      
      if(curr_viable[1] > overlap[2] & overlap[1] > curr_viable[2]){
        
        # Swapping the first column of the two rows
        temp <- curr_viable[1]
        means[rand,1] <- overlap[1]
        means[to_fix[index],1] <- temp
        
        left_to_fix <- left_to_fix - 1
        prop.allowed <- left_to_fix == 0
        
        if(prop.allowed){
          attributes(means)$swap_count <- count
          return(means)
        }
        else
          break
      }
      
      # If the above block does not trigger, then viable will be updated
      # to exclude the just-used random selection
      viable <- viable[which(!viable %in% rand)]
      count <- count + 1
    }
  }
  # If the algorithm does not generate a viable matrix
  return(NULL)
}


#' Helper function used by baseSimulation
#' This will be called to deal with phage interactions with hosts, in order to
#' calculate the fitnesses of children and grandchildren.
findLyseFitness <- function(outcomes, nA, nB, lambda, alpha, omega, theta, 
                            p, is_specialist) {
  # found_theta represents a phage interacting with a host
  found_theta <- which(outcomes %in% c(2, 3))
  
  # The loop covers both specialist and generalist by dealing with
  # outcomes having 2 and 3 via an if else statement. The fitnesses of
  # children and grandchildren are then found and put into outcomes, which is 
  # then returned.
  for(k in found_theta){
    if(outcomes[k] == 2){
      burst_size <- rpois(1, nA + 1)
      if(runif(1) < omega){
        outcomes_lyse <- NULL
        
        if(is_specialist){
          outcomes_lyse <- sample(c(0, 1, nA), burst_size, replace=TRUE,
                                  prob=c(lambda, alpha, theta*p))
        } else {
          outcomes_lyse <- sample(c(0, 1, nA, nB), burst_size, replace=TRUE,
                                  prob=c(lambda, alpha, theta*p, theta*(1 - p)))
        }
        outcomes[k] <- sum(outcomes_lyse)
      } else {
        outcomes[k] <- burst_size - 1
      }
    } else if(outcomes[k] == 3){
      burst_size <- rpois(1, nB + 1)
      
      if(runif(1) < omega){
        outcomes_lyse <- sample(c(0, 1, nA, nB), burst_size, replace=TRUE, 
                                prob=c(lambda, alpha, theta * p, theta*(1-p)))
        outcomes[k] <- sum(outcomes_lyse)
      } else {
        outcomes[k] <- burst_size - 1
      }
    }
  }
  return(outcomes)
}


#' Helper function used by zeroHunter
#' Finds the point where wG = wS and return that point by solving the 
#' linear slope equations for both.
slopeEqual <- function(min_x, min_wg, min_ws, max_x, max_wg, max_ws) {
  
  g_slope <- (max_wg - min_wg) / (max_x - min_x)
  s_slope <- (max_ws - min_ws) / (max_x - min_x)
  
  intercept_g <- min_wg - (g_slope * min_x)
  intercept_s <- min_ws - (s_slope * min_x)
  
  
  # Now solving g_slope*x + intercept_g = s_slope*x + intercept_s
  x <- (intercept_s - intercept_g) / (g_slope - s_slope)
  
  return(x / nA) # Dividing nB by nA will return R*
}




#' Helper function used in runChanging AND standalone function
#' Generates two vectors representing W_S and W_G for the given parameters.
#' Used by equilibrium simulation.
#' RETURN: a matrix where the first column is the fitness vector for W_S and the
#'  second is a fitness_vector for W_G
baseSimulation <- function(nA, alpha, lambda, omega, theta, p, inc_past=0.15, 
                           reps=5e5) {
  # Simulated fitnesses
  base <- c()
  
  runs <- 0:(nA + as.integer(nA*inc_past))
  
  for(i in runs){
    nB <- i
    
    fit_wS <- sample(c(0, 1, 2), reps, replace=TRUE, prob=c(lambda, alpha, theta*p))
    fit_wG <- sample(c(0, 1, 2, 3), reps, replace=TRUE, prob=c(lambda, alpha, theta*p, theta*(1-p)))
    
    fit_wS <- findLyseFitness(fit_wS, nA=nA, nB=nB, lambda=lambda, alpha=alpha, 
                              omega=omega, theta=theta, p=p, is_specialist=TRUE)
    fit_wG <- findLyseFitness(fit_wG, nA=nA, nB=nB, lambda=lambda, alpha=alpha, 
                              omega=omega, theta=theta, p=p, is_specialist=FALSE)
    
    base <- cbind(base, fit_wS, fit_wG)
  }
  
  return(base)   
}


#' Standalone function
#' Generates a nA by 2 matrix which contains the predictions for W_s and W_g of 
#' the base simulation.
baseSimPrediction <- function(nA, alpha, theta, p, lambda, omega, 
                              inc_past=0.15) {
  
  nB <- 0:(nA + as.integer(nA*inc_past))

  prediction_wS <- alpha/(alpha + lambda + theta*p) + 
    theta*p/(alpha + lambda + theta*p)*((1 - omega)*nA + 
                                          omega*(nA + 1)*(alpha/(alpha + lambda + theta*p) + 
                                                            nA*theta*p/(alpha + lambda + theta*p)))
  
  prediction_wG <- alpha/(alpha + lambda + theta) + 
    theta*p/(alpha + lambda + theta)*((1 - omega)*nA + 
                                        omega*(nA + 1)*(alpha/(alpha + lambda + theta) + 
                                                          theta/(alpha + lambda + theta) * (p * nA + (1 - p) * nB))) +
    theta*(1 - p)/(alpha + lambda + theta)*((1 - omega)*nB + 
                                              omega*(nB + 1)*(alpha/(alpha + lambda + theta) + 
                                                                theta/(alpha + lambda + theta) * (p * nA + (1 - p) * nB)))
  
  
  return(cbind(prediction_wS, prediction_wG))
}


#' Standalone function that plots the results generated by baseSimulation
#' and baseSimPrediction
plotBaseSimulation <- function(base, prediction=NA, with.prediction=FALSE, 
                               colors=c('coral4', 'darkorchid4')) {
  
  curr_wS <- colMeans(base[,seq(from=1, to=ncol(base), by=2)])
  curr_wG <- colMeans(base[,seq(from=2, to=ncol(base), by=2)])

  y_min <- min(curr_wG)
  y_max <- max(curr_wG)
  
  if(max(curr_wS) > max(curr_wG))
    y_max <- max(curr_wS)
  
  plot(y=curr_wS, x=1:length(curr_wS), type='p', col=colors[1], 
       xlim=c(1, length(curr_wS)), ylim=c(y_min, y_max))
  
  points(y=curr_wG, x=1:length(curr_wG), col=colors[2])
  
  if(with.prediction){
    if(is.na(prediction))
      stop("Needs a prediction matrix")
    lines(y=prediction[1],x=1:length(curr_wS),col=colors[1], lwd=1.5)
    lines(y=prediction[2],x=1:length(curr_wG), col=colors[2], lwd=1.5)
  }
}





basicBootstrap <- function(g, s){
  bt_size <- 1e4
  
  size <- length(g)
  
  dfG <- as.data.frame(table(g))
  uniqueG <- as.numeric(levels(dfG[,1]))
  freqG <- dfG[,2] / size
  
  dfS <- as.data.frame(table(s))
  uniqueS <- as.numeric(levels(dfS[,1]))
  freqS <- dfS[,2] / size
  
  g_means <- replicate(bt_size, mean(sample(uniqueG, size, replace=T, prob=freqG)))
  s_means <- replicate(bt_size, mean(sample(uniqueS, size, replace=T, prob=freqS)))
  
  
  return(cbind(s_means, g_means))
}





preprocessed <- function(files, changing_name, changing) {
  
  exist <- which(sapply(files, file.exists))
  if(length(exist) != 0) stop(exist)
  
  
  
}

rStar <- function(file, current_changing) {
  boot <- read.table(file, header=TRUE, sep=",")
  
  g <- boot[,seq(2, ncol(boot), 2)]
  s <- boot[,seq(1, ncol(boot), 2)]
  
  viable_lower <- c()
  viable_upper <- c()
  for(i in 1:ncol(g)){
    s_greater <- which(s[,i] > g[,i])
    g_greater <- which(g[,i] > s[,i])
    
    if(length(s_greater) / nrow(boot) <= prop_allowed)
      viable_lower <- append(viable_lower, i)
    else if(length(g_greater) / nrow(boot) <= prop_allowed)
      viable_upper <- append(viable_upper, i)
  }
  
  
  lower_boot <- NULL
  upper_boot <- NULL
  r_star <- NULL
  done <- FALSE
  
  while(!done){
    if(length(viable_upper) == 0){
      message(paste("No upper found for ", file, "\n",
                    "Setting R* equal to 1", sep=""))
      r_star <- 1
    } else if(length(viable_lower) == 0){
      message(paste("No lower found for ", file, "\n",
                    "Setting R* equal to 0", sep=""))
      r_star <- 0
    }
    
    curr_upper <- min(viable_upper)
    curr_lower <- max(viable_lower)
    
    lower_
    
  }
  
  
}

#' Technical overview
#' A single fitness vector is has a length of 500,000 and is created using sample().
#' Each entry on the vector represents the fitness of a unique phage. If the phage dies, 
#' it has a fitness of 0, if it leaves the aggregation it has a fitness of 1, 
#' and if it experiences a lyse event, it has a fitness determined by the function findLyseFitness,
#' which will calculate the fitness of the children and grandchildren of the phage, should it 
#' burst inside the aggregation. If the phage does not lyse inside the aggregation, then the fitness
#' of the phage is simply equal to the burst size of the infected host. To test our hypotheses, 
#' the net burst size of host B, nB, ranges from 0 to the net burst size of host A, nA, plus nA * .15.
#' Once a suitable range of fitness-vectors have been found for various parameters, the fitness vectors are
#' bootstrapped by taking the mean of a fitness vector that has had sample used on it once again. This is done 10,000
#' times for each fitness vector.
#' 
#' We consider a grouping of fitness vectors for both specialists and generalists where nB is the only changing value
#' to be a set. The ratio nB/nA is found for each set by calculating the numerical point where W_S = W_G and dividing 
#' this by nA. The fitness lines are assumed to be linear between close enough points; the upper point is the first
#' bootstrapped specialist/generalist fitness vector pair where each value for the generalist is greater than
#' its opposing specialist value. The opposite is true for the lower point. If the amount of unwanted overlap is less than 
#' or equal to 1%, then an algorithm is called that randomly cycles through the vectors and tries to swap points
#' so that the found overlap is 0%. If this is a success, then the swapped-up pair of vectors will be used.
#' 
#' Since linearity is assumed, the linear slope equation will be used to find the value of nB for when the specialist
#' fitness equals the generalist fitness, and the mean of the resulting vector is divided by nA to yield R*. 95% quantiles
#' are generated before the mean is taken.
