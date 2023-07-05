burst_A <- 30 # burst size of host A
inc_past <- 0.15 # how far past nA nB goes (limites failures)
reps <- 5e5 # how large each fitness vector is
bt_size <- 1e4 # how many times each fitness vector is bootstrapped
prop_allowed <- 0.01 # proportion of allowed overlap between bootstrapped vectors



#' Helper function used throughout
#' Used to verify floating numbers
floatMatch <- function(x, y) {
  return(abs(x - y) < 1e-6)
}


#' This attempts to swap overlapping locations with non-overlapping locations 
#' from the matrix so that there is no overlap at all. The first column is 
#' expected to be larger than the second; locations where the second is larger
#' than the first is considered to be overlap.
swap <- function(data){
  
  # Vector of locations inside data where the second column is > than the first
  to_fix <- which(data[,1] <= data[,2])
  left_to_fix <- length(to_fix)

  if(left_to_fix == 0) return(data)
  
  
  for(index in to_fix){
    # The current data location trying to be fixed.
    overlap <- data[index,]
    
    # Locations that could potentially swap with the current overlap row
    viable <- which(!1:nrow(data) %in% to_fix)

    # Records how many times a redraw happens
    count <- 1 
    
    while(length(viable) > 0){
      
      rand <- sample(viable, 1)
      
      curr_viable <- data[rand,]
      
      if(curr_viable[1] > overlap[2] & overlap[1] > curr_viable[2]){
        
        # Swapping the first column of the two rows
        temp <- curr_viable[1]
        data[rand,1] <- overlap[1]
        data[to_fix[index],1] <- temp
        
        left_to_fix <- left_to_fix - 1

        if(left_to_fix > 0){
          break
        }
        attributes(data)$swap_count <- count
        print(attributes(data))
        return(data)
      }
      # If the above block does not trigger, then viable will be updated
      # to exclude the just-used random selection
      viable <- viable[which(!viable %in% rand)]
      count <- count + 1
    }
  }
  return(NULL) # If the algorithm does not generate a viable matrix
}


#' Helper function used by baseSimulation
#' This will be called to deal with phage interactions with hosts, in order to
#' calculate the fitnesses of children and grandchildren.
findLyseFitness <- function(outcomes, burst_B, alpha, theta, p, lambda, omega, 
                            is.specialist) {
  # found_theta represents a phage interacting with a host
  found_theta <- which(outcomes %in% c(2, 3))
  if(length(found_theta) == 0)
    return(outcomes)
  
  burst_chance <- rbinom(found_theta, 1, omega)
  
  # The loop covers both specialist and generalist by dealing with
  # outcomes having 2 and 3 via an if else statement. The fitnesses of
  # children and grandchildren are then found and put into outcomes, which is 
  # then returned.
  for(k in 1:length(found_theta)){
    loc <- found_theta[k]
    
    if(outcomes[loc] == 2){
      burst_size <- rpois(1, burst_A)
      if(burst_chance[k] == TRUE){
        outcomes_lyse <- NULL
        
        if(is.specialist){
          outcomes_lyse <- sample(c(0, 1, burst_A), burst_size, replace=TRUE,
                                  prob=c(lambda, alpha, theta*p))
        }else {
          outcomes_lyse <- sample(c(0, 1, burst_A, burst_B), burst_size, replace=TRUE,
                                  prob=c(lambda, alpha, theta*p, theta*(1 - p)))
        }
        outcomes[loc] <- sum(outcomes_lyse)
      }else {
        outcomes[loc] <- burst_size
      }
    }else if(outcomes[loc] == 3){
      burst_size <- rpois(1, burst_B)
      if(burst_chance[k] == FALSE){
        outcomes_lyse <- sample(c(0, 1, burst_A, burst_B), burst_size, replace=TRUE, 
                                prob=c(lambda, alpha, theta * p, theta*(1-p)))
        outcomes[loc] <- sum(outcomes_lyse)
      }else {
        outcomes[loc] <- burst_size
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
  
  return(x / burst_A) # Dividing burst_B by burst_A will return R*
}


#' Helper function used in runChanging AND standalone function
#' Generates two vectors representing W_S and W_G for the given parameters.
#' Used by equilibrium simulation.
#' RETURN: a matrix where the first column is the fitness vector for W_S and the
#'  second is a fitness_vector for W_G
baseSimulation <- function(alpha, theta, p, lambda, omega) {
  
  #Effective fitness: E(W_S) = (theta*p*W_S + alpha)/alpha + lambda + theta*p
  
  # Simulated fitnesses
  base <- c()
  
  runs <- 0:(burst_A + as.integer(burst_A*inc_past))
  
  for(i in runs){
    burst_B <- i
    
    fit_wS <- sample(c(0, 1, 2), reps, replace=TRUE, prob=c(lambda, alpha, theta*p))
    fit_wG <- sample(c(0, 1, 2, 3), reps, replace=TRUE, prob=c(lambda, alpha, theta*p, theta*(1-p)))
    
    fit_wS <- findLyseFitness(fit_wS, burst_B=burst_B, lambda=lambda, alpha=alpha, 
                              omega=omega, theta=theta, p=p, is.specialist=TRUE)
    fit_wG <- findLyseFitness(fit_wG, burst_B=burst_B, lambda=lambda, alpha=alpha, 
                              omega=omega, theta=theta, p=p, is.specialist=FALSE)
    
    base <- cbind(base, fit_wS, fit_wG)
  }
  
  return(base)   
}


#' Standalone function
#' Generates a nA by 2 matrix which contains the predictions for W_s and W_g of 
#' the base simulation.
baseSimPrediction <- function(alpha, theta, p, lambda, omega) {
  
  burst_B <- 0:(burst_A + as.integer(burst_A*inc_past))
  
  prediction_wS <- alpha/(alpha + lambda + theta*p) + 
    theta*p/(alpha + lambda + theta*p)*((1 - omega)*(burst_A - 1) + 
                                          omega*((burst_A - 1) + 1)*(alpha/(alpha + lambda + theta*p) + 
                                                            (burst_A - 1)*theta*p/(alpha + lambda + theta*p)))
  
  prediction_wG <- alpha/(alpha + lambda + theta) + 
    theta*p/(alpha + lambda + theta)*((1 - omega)*(burst_A - 1) + 
                                        omega*((burst_A - 1) + 1)*(alpha/(alpha + lambda + theta) + 
                                                          theta/(alpha + lambda + theta) * (p * (burst_A - 1) + (1 - p) * (burst_B- 1)))) +
    theta*(1 - p)/(alpha + lambda + theta)*((1 - omega)*(burst_B- 1) + 
                                              omega*((burst_B- 1) + 1)*(alpha/(alpha + lambda + theta) + 
                                                                theta/(alpha + lambda + theta) * (p * (burst_A - 1) + (1 - p) * (burst_B- 1))))
  data <- cbind(prediction_wS, prediction_wG)
  
  plot(x=burst_B, y=data[,2], type='l', col='firebrick', lwd=1.5, ylab="fitness")
  lines(x=burst_B, y=data[,1], col='darkblue', lwd=1.5)
  legend("topleft", legend=c("fitness.s", "fitness.g"), lty=1, 
         col=c('darkblue', 'firebrick'), lwd=1.5)
  return(data)
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
    lines(y=prediction[1],x=1:length(curr_wS),col=colors[1], lwd=1.5)
    lines(y=prediction[2],x=1:length(curr_wG), col=colors[2], lwd=1.5)
  }
}

#' Generates a pair of bootstrapped fitness vectors
basicBootstrap <- function(s, g){

  size <- length(g)
  
  dfG <- as.data.frame(table(g))
  uniqueG <- as.numeric(levels(dfG[,1]))
  freqG <- dfG[,2] / size
  
  dfS <- as.data.frame(table(s))
  uniqueS <- as.numeric(levels(dfS[,1]))
  freqS <- dfS[,2] / size
  
  g_means <- replicate(bt_size, sum(rmultinom(1, size, freqG) * uniqueG) / size)
  s_means <- replicate(bt_size, sum(rmultinom(1, size, freqS) * uniqueS) / size)

  # g_means <- replicate(bt_size, mean(sample(uniqueG, size, replace=T, prob=freqG)))
  # s_means <- replicate(bt_size, mean(sample(uniqueS, size, replace=T, prob=freqS)))
  
  return(cbind(s_means, g_means))
}


# c("lower.bound", "lower.boot.s", "lower.boot.g", "upper.bound", 
#   "upper.boot.s", "upper.boot.g", changing_name, "lower.quantile", 
#   "upper.quantile")

#' Returns a matrix with 9 columns and number of files * bt_size rows
#' The columns are:
#'  lower.bound: the nB where the lower boundry was found
#'  lower.boot.s: the bootstrapped fitness vector for S at the lower bound
#'  lower.boot.g: same as above except for G
#'  upper.bound: the nB where the upper boundry was found.
#'  upper.boot.s: the bootstrapped fitness vector for S at the upper bound
#'  upper.boot.g: same as above except for G
#'  r.star: the value of nB where WG = WS. Solved for using the linear slope equation
#'  r.star.mean: the mean of r.star
#'  changing: the current value of the changing variable
#'  lower.quantile: lower bound of a 95% confidence interval for the r.star vector
#'  upper.quantile: upper bound of a 95% confidence interval for the r.star vector
preprocessed <- function(files, changing, changing_name) {

  exist <- which(sapply(files, file.exists))
  if(length(exist) != length(files)) 
    stop(files[1:length(files)[!1:length(files) %in% exist]])
  if(length(files) != length(changing))
    stop("mismatch")
  

  data <- list()
  for(f in 1:length(files)){
    data[[f]] <- rStar(files[f], changing[f], changing_name)
  }

  return(data)
}
  
#' HEADER=FALSE
#' Finds R* and all acompanying data and returns it as a matrix
rStar <- function(file, current_changing, changing_name) {
  boot <- read.table(file, header=TRUE, sep=",")
  
  # Files are organized by repeating the sequence W.S W.G
  s <- boot[,seq(1, ncol(boot), 2)]
  g <- boot[,seq(2, ncol(boot), 2)]
  
  # Containers for possibly useful upper/lower bounds in the file from which
  # R* can be found
  viable_lower <- c()
  viable_upper <- c()
  
  for(i in 1:ncol(g)){
    s_greater <- which(s[,i] > g[,i])
    g_greater <- which(g[,i] > s[,i])
    
    if(1 - (length(s_greater) / nrow(boot)) <= prop_allowed)
      viable_lower <- append(viable_lower, i)
    else if(1 - (length(g_greater) / nrow(boot)) <= prop_allowed)
      viable_upper <- append(viable_upper, i)
  }
  
  # Upper and lower points used to calculate R* via the linear slope equation
  lower_pair <- NULL
  upper_pair <- NULL
  r_stars <- NULL
  lower_burst_B <- NULL
  upper_burst_B <- NULL
  
  while(is.null(lower_pair)){
    if(length(viable_lower) == 0){
      message(paste("No lower found for ", file, "\n",
                    "Setting R* equal to 0", sep=""))
      r_stars <- 0
      lower_burst_B <- -1
      break
    }
    lower_burst_B <- max(viable_lower)
    lower_pair <- swap(cbind(s[,lower_burst_B], g[,lower_burst_B]))
    viable_lower <- viable_lower[!viable_lower %in% lower_burst_B]
  }
  
  while(is.null(upper_pair)){
    if(length(viable_upper) == 0){
      message(paste("No upper found for ", file, "\n",
                    "Setting R* equal to 1", sep=""))
      r_stars <- 1
      upper_burst_B <- burst_A + as.integer(burst_A*inc_past) + 1
      break
    }
    upper_burst_B <- min(viable_upper)
    upper_pair <- swap(cbind(g[,upper_burst_B], s[,upper_burst_B]))
    viable_upper <- viable_upper[!viable_upper %in% upper_burst_B]
  }
  
  if(is.null(r_stars)){
    r_stars <- mapply(slopeEqual, lower_burst_B, lower_pair[,2], lower_pair[,1],
                      upper_burst_B, upper_pair[,1], upper_pair[,2])

  }
  
  quant <- quantile(r_stars, probs=c(0.025, 0.975))
  lower_quantile <- quant[1]
  upper_quantile <- quant[2]
  
  data <- list(lower.bound = lower_burst_B, lower.boot.s = lower_pair[,1], 
               lower.boot.g = lower_pair[,2], upper.bound = upper_burst_B,
               upper.boot.s = upper_pair[,2], upper.boot.g = upper_pair[,1], 
               r.star = r_stars, r.star.mean = mean(r_stars), 
               changing.name = current_changing, lower.quantile = lower_quantile, 
               upper.quantile = upper_quantile, 
               swap.count.lower = attributes(lower_pair)$swap_count, 
               swap.count.upper = attributes(upper_pair)$swap_count)
  return(data)
}


checkOverlap <- function(s, g) {
  out <- mapply(function(x, y) {
    out <- matrix(ncol=2, nrow=1)
    colnames(out) <- c("S.greater", "G.greater")
    
    out[1,1] <- length(which(x > y))
    out[1,2] <- bt_size - out[1,1]
    if(out[1,1] + out[1,2] != bt_size)
      message("at least one pair of values are equal")
    return(out)
  }, s, g, SIMPLIFY=FALSE)
  
  return(do.call(rbind, out))
}

plotFitness <- function(fitness) {
  change <- unique(fitness[,9])
  print(change)
  plotting <- c()
  
  for(cur in change){
    plotting <- rbind(plotting, fitness[which(floatMatch(fitness[,9], cur))[1],])
  }
  
  plotting <- as.data.frame(plotting)

  plot <- ggplot2::ggplot(plotting, mapping=aes(x=plotting[,9], y=plotting[,8])) +
    geom_point() +
    geom_line() +
    theme_classic() +
    geom_errorbar(aes(ymin=plotting[,10], ymax=plotting[,11]), width=0.0002)
  
  return(plot)
}

#' Technical overview
#' A single fitness vector is has a length of 500,000 and is created using 
#' sample(). Each entry on the vector represents the fitness of a unique phage. 
#' If the phage dies, it has a fitness of 0, if it leaves the aggregation it has
#' a fitness of 1, and if it experiences a lyse event, it has a fitness 
#' determined by the function findLyseFitness, which will calculate the fitness 
#' of the children and grandchildren of the phage, should it burst inside the 
#' aggregation. If the phage does not lyse inside the aggregation, then the 
#' fitness of the phage is simply equal to the burst size of the infected host. 
#' To test our hypotheses, the net burst size of host B, nB, ranges from 0 to
#' the net burst size of host A, nA, plus nA * .15. Once a suitable range of 
#' fitness-vectors have been found for various parameters, the fitness vectors
#' are bootstrapped by taking the mean of a fitness vector that has had sample 
#' used on it once again. This is done 10,000 times for each fitness vector.
#' 
#' We consider a grouping of fitness vectors for both specialists and 
#' generalists where nB is the only changing value to be a set. The ratio nB/nA
#' is found for each set by calculating the numerical point where W_S = W_G and 
#' dividing this by nA. The fitness lines are assumed to be linear between close
#' enough points; the upper point is the first bootstrapped specialist/generalist 
#' fitness vector pair where each value for the generalist is greater than
#' its opposing specialist value. The opposite is true for the lower point. If 
#' the amount of unwanted overlap is less than or equal to 1%, then an algorithm
#' is called that randomly cycles through the vectors and tries to swap points
#' so that the found overlap is 0%. If this is a success, then the swapped-up 
#' pair of vectors will be used.
#' 
#' Since linearity is assumed, the linear slope equation will be used to find 
#' the value of nB for when the specialist fitness equals the generalist 
#' fitness, and the mean of the resulting vector is divided by nA to yield R*. 
#' 95% quantiles are generated before the mean is taken.