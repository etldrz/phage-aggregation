
#' Helper function used throughout
#' For verifying floating numbers
floatMatch <- function(x, y) {
  return(abs(x - y) < 1e-6)
}


#' Helper function used by equilibriumSimulation
#' When reverse is false, then this will check to see if each of the generated g_means
#' is less than the generated s_means, and will return false if this is the case.
#' A return of false indicates to the calling fucntion, equilibriumSimulation,
#' that it needs to keep searching for viable vectors. If all(g_means > s_means), 
#' then a matrix of the two vectors will be returned. A return of a matrix
#' indicates to equilibriumSimulation that the returned bootstrapped vectors are good,
#' and it will start to go backwards, to find where all(g_means < s_means). reverse
#' will be true in this case.
bootstrap <- function(fitness_g, fitness_s, bt_size=1e4, reverse=FALSE,
                      allowed_overlap=0.1, gather.overlap=FALSE) {
  
  g_means <- rep(0, bt_size)
  s_means <- rep(0, bt_size)
  
  # indicates the amount of overlapping that takes place; used when gather.overlap = T
  g_overlap <- 0
  s_overlap <- 0
  
  
  for (i in 1:bt_size) {
    g_means[i] <- mean(sample(fitness_g, length(fitness_g), replace=TRUE))
    s_means[i] <- mean(sample(fitness_s, length(fitness_s), replace=TRUE))
    
    if(gather.overlap){
      if(g_means[i] < s_means[i]){
        s_overlap <- s_overlap + 1
      }
      else if(s_means[i] < g_means[i]){
        g_overlap <- g_overlap + 1
      }
    }
  }
  
  if(gather.overlap)
    return(cbind(s_overlap, g_overlap))
  
  # reverse determines the order of the comb_means matrix, it it is true
  # then s_means is first and g_means is first if it is true.
  curr_overlap <- (length(which(s_means > g_means)) / bt_size)*100
  comb_means <- cbind(g_means, s_means)
  if(reverse){
    curr_overlap <- (length(which(g_means > s_means)) / bt_size)*100
    comb_means <- cbind(s_means, g_means)
  }
  
  
  if(curr_overlap <= allowed_overlap){
    return(swap(comb_means))
  }
  else return(NULL)
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
  perc.allowed <- FALSE
  
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
        perc.allowed <- left_to_fix == 0
        
        if(perc.allowed)
          return(means)
        else
          break
      }
      
      # If the above block does not trigger, then viable will be updated
      # to exclude the just-used random selection
      viable <- which(!viable %in% rand)
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


#' Helper function used in runChanging
#' This will find the equilibrium point for each fitness$changing_value by using
#' the helper function slopeEqual. It will then find the upper quantile, the 
#' lower quantile, and the mean of the found R* and return all four of these 
#' vectors as a matrix.
zeroHunter <- function(fitness) {
  r_star <- rep(NA, nrow(fitness))
  
  for(i in 1:nrow(fitness)){
    r_star[i] <- slopeEqual(min_x = fitness$x_lower[i], min_wg = fitness$bt_wg_lower[i], 
                            min_ws = fitness$bt_ws_lower[i], max_x = fitness$x_upper[i], 
                            max_wg = fitness$bt_wg_upper[i], max_ws = fitness$bt_ws_upper[i])
  }
  
  quant <- quantile(r_star, probs=c(0.025, 0.975))
  
  lower_quantile <- rep(quant[1], nrow(fitness))
  upper_quantile <- rep(quant[2], nrow(fitness))
  r_star_mean <- rep(mean(r_star), nrow(fitness))
  combined <- cbind(lower_quantile, upper_quantile, r_star, r_star_mean)
  return(combined)
}


#' TODO: make sure it deal with changing parameters other than alpha and theta
#' Helper function used in runChanging
#' This runs the simulation and finds all the necessary data to find R* and the 
#' bootstrapped confidence intervals. base_outcomes are two matrices which hold nA columns,
#' one for each point in the base plot.
#' changing: a string of the current changing variable name
#' bt_size: length of the bootstrap vector
equilibriumSimulation <- function(plot_file, nA, omega, alpha, theta, p, lambda, 
                                  inc_past=0.25, fitness_s, fitness_g, changing, 
                                  allowed_overlap, bt_size) {
  
  # The names of the dataframe generated by equilibriumSimulaiton
  legend <- c('x_lower', 'bt_ws_lower','bt_wg_lower','x_upper', 'bt_ws_upper', 
              'bt_wg_upper', 'lower_quantile', 'upper_quantile', 'r_star',
              'r_star_mean', 'changing_value', 'lambda', 'omega', 'p', 'prediction')
  
  
  # These hold the fitnesses during the length of total_treatments; the fitness vectors
  # that are not needed by the end are discarded.
  g_outcomes_dataframe <- c()
  s_outcomes_dataframe <- c()
  
  current_value <- theta
  if(changing == "alpha") 
    current_value <- alpha
  
  nB_vec <- 0:(nA + as.integer(nA*inc_past))
  
  total_treatments <- length(nB_vec)
  
  
  # When this is true, the parent loop will know to break. It will be true once
  # all four fitness vectors are found.
  halt <- FALSE
  
  for(i in nB_vec){
    if(halt)
      break
    
    # The nB value to be used in this loop
    nB <- i
    
    # The mathematical prediction for the current parameters.
    prediction <- 1 / (1 + (alpha + lambda) / (theta * p))
    
    # This holds all of the relevant data that will be saved and used later.
    bound <- cbind(NA, NA, NA, 
                   NA, NA, NA, 
                   NA, NA, NA, NA,
                   rep(current_value, bt_size), rep(lambda, bt_size), 
                   rep(omega, bt_size), rep(p, bt_size), rep(prediction, bt_size))
    
    # Bound legend:
    #             'x_lower', 'bt_ws_lower', 'bt_wg_lower', 
    #             'x_upper', 'bt_ws_upper', 'bt_wg_upper', 
    #             'lower_quantile', 'upper_quantile', 'r_star', 'r_star_mean',
    #             'changing_value', 'lambda',
    #             'omega', 'p', 'prediction'
    
    
    # Will return a matrix if a viable one is found, null if not. The order of 
    # the columns for the matrix will be g,s for !reverse and s,g for reverse.
    # Goes until a non-null matrix is found, or the upper bound is reached. Same
    # for finding the lower bound.
    bootstrapped <- bootstrap(fitness_g[,i], fitness_s[,i], bt_size=bt_size, 
                              allowed_overlap=allowed_overlap, reverse=FALSE)
    
    
    if(is.null(bootstrapped)){
      
      # If fitness_g and fitness_s are not viable, then they will be saved
      # so that they can be accessed when reverse=TRUE
      g_outcomes_dataframe <- cbind(g_outcomes_dataframe, fitness_g[,i])
      s_outcomes_dataframe <- cbind(s_outcomes_dataframe, fitness_s[,i])
      
      # This check just makes sure that the user is notified when the simulation
      # doesn't find a value when it should.
      if(i < total_treatments){
        next
      }
      else{
        stop("Reached the upper bound without finding a viable vector", 
                ", ", current_value)
      }
    } 
    else if(is.matrix(bootstrapped)){
      print(paste("Uppermost found at", i, "out of", total_treatments))
      
      # Using bound to record the bootstrapped viable upper-limit fitness vectors,
      # as well as the nB they were found at. Dividing by nA so R* will be found.
      bound[,4:6] <- cbind(rep(nB, bt_size), 
                           bootstrapped[,2], # bootstrapped specialist
                           bootstrapped[,1]) # bootstrapped generalist
      
      for(y in (nB - 1):1) {
        
        bootstrapped <- bootstrap(g_outcomes_dataframe[,y],
                                  s_outcomes_dataframe[,y], bt_size=bt_size,
                                  allowed_overlap=allowed_overlap, reverse=TRUE)
        
        # Will continue reversing if bootstrap() says the vectors aren't viable by
        # returning a null value
        if(is.null(bootstrapped) & y == min(nB_vec)){
          stop("Reached the lower bound without finding a viable vector", ", ", current_value)
        }
        else if(is.null(bootstrapped) & y > min(nB_vec)){
          next
        }
        
        # Place the lower entries into the bound vector
        bound[,1:3] <- cbind(rep(y, bt_size), #y is the found lower value of nB
                             bootstrapped[,1], bootstrapped[,2])
        
        bound <- as.data.frame(bound)
        
        names(bound) <- legend

        # These ifs deal with moving the values into a text file by checking to see
        # if the text file is empty or not
        if(file.info(plot_file)$size != 0){
          
          # If the text file isn't empty, then the data will be pulled from it,
          # and bound will be appended. This is to ensure smooth transition of 
          # data legibility.
          fitness_bootstrap <- read.table(file=plot_file, header=TRUE, sep=",")
          names(fitness_bootstrap) <- legend
          
          fitness_bootstrap <- rbind(fitness_bootstrap, bound)
          write.table(fitness_bootstrap, file=plot_file, append=FALSE, sep=",",
                      quote=FALSE)
          
          # When halt is true, the parent loop of equilibriumSimulation() will break.
          halt <- TRUE
          break
        } else {
          
          fitness_bootstrap <- bound
          names(fitness_bootstrap) <- legend
          write.table(fitness_bootstrap, file=plot_file, append=FALSE, sep=",", quote=FALSE)
          halt <- TRUE
          break
        }
      }
    }
  }
} #end of equilibriumSimulation


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
                              inc_past=0.25) {
  
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


#' Standalone function
#' TODO: make sure that baseSimPrediction is dealt with somehow
#' TODO: fix the theta*p logic
#' This will run a simulation for some changing value for each lambda and omega inputted; it can
#' be either a vector or a single value.
#' plot_file: where the data is stored as a txt file.
#'  TODO: put a comment that has the parameter values at the top of the file.
#' changing_name: the name of the variable that is being changed, either "alpha" or "thetaP"
#' changing_step_size: by how much the changing value increments
#' change_start: where the changing value starts its sequence at
#' changing_length: the length of the changing value sequence
#' nA, alpha, theta, p, lambda, and omega are parameters for the model
#' save_base: when true, the base data will be saved alongside the equilibriumSimulation data
#'  TODO: implement something that will make the base.txt file easy to read and access; it 
#'  currently just stores a mess of data with no parameter values being stored.
#' base_file: where the base data is saved.
#' bt_size: how many times the fitness vectors from equilibriumSimulation are bootstrapped
#' reps: how many replicates are run.
#' RETURN: the fitness dataframe, which can also be found at plot_file
runChanging <- function(plot_file, changing_name, changing_step_size, 
                        change_start, changing_length=10, nA, alpha, theta, p, 
                        lambda, omega, inc_past=0.15, bt_size=1e4, 
                        reps=5e5) {

  # The changing value
  changing <- seq(from=change_start, by=changing_step_size, length.out=changing_length)
  
  if(max(changing) > 1)
    stop("Illegal probability vector")
  
  # Where the base data is stored for access by equilibriumSimulation
  base <- c()
  
  for(curr in changing) {
    
    # Sets the current changing value
    if(changing_name == "alpha")
      alpha <- curr
    else if(changing_name == "theta")
      theta <- curr
    else 
      stop("Invalid changing name")
    
    
    for(l in lambda){
      for(o in omega){
        for(prop in p){
          current_base <- baseSimulation(nA=nA, alpha=alpha, lambda=l,
                                         omega=o, theta=theta, p=prop, reps=reps)
          
          current_base_wS <- current_base[,seq(from=1, to=ncol(current_base), by=2)]
          current_base_wG <- current_base[,seq(from=2, to=ncol(current_base), by=2)]
          
          base <- cbind(base, current_base)
          equilibriumSimulation(plot_file=plot_file, nA=nA, omega=o, alpha=alpha,
                                theta=theta, p=prop, lambda=l, fitness_s=current_base_wS,
                                fitness_g=current_base_wG, changing=changing_name, bt_size=bt_size)
        }
      }
    }
  }
  
  fitness <- read.table(file=plot_file, sep=",", header=TRUE)
  
  # This finds where W_g = W_s and finishes filling the main dataframe by filling r_star,
  # r_star_mean, lower_quantile, and upper_quantile.
  for(curr in changing){
    for(l in lambda){
      for(o in omega){
        for(per in p){
          current <- which(floatMatch(fitness$changing_value, curr) & floatMatch(fitness$lambda, l) &
                             floatMatch(fitness$omega, o) & floatMatch(fitness$p, per))
          
          final_values <- as.data.frame(zeroHunter(fitness[current,]))
          
          fitness$lower_quantile[current] <- final_values$lower_quantile
          fitness$upper_quantile[current] <- final_values$upper_quantile
          fitness$r_star[current] <- final_values$r_star
          fitness$r_star_mean[current] <- final_values$r_star_mean
        }
      }
    }
  }
  write.table(fitness, file=plot_file, sep=",", append=FALSE)
  
  return(fitness)
} #end of runChanging


safeRunChanging <- function(safety=1.5e6, reps=5e5, ...) {
  result <- NA
  safety <- reps*2
  while(reps != safety){
    if(is.data.frame(result))
      break
    
    result <- tryCatch({
      runChanging(..., reps=reps)
    }, error = function(err){
      message(err)
      return(NA)
    }, warning = function(war){
      message(war)
    })
    
    reps <- reps + 5e5
    print(paste("reps is ", reps))
  }
  
  if(!is.data.frame(result)){
    print("nope")
  }
  
  return(result)
}


gatherOverlap <- function(nA, alpha, theta, p, lambda, omega,
                          bt_size=1e4, reps=5e5) {
  
  s_over <- c()
  g_over <- c()
  
  base <- baseSimulation(nA=nA, alpha=alpha, theta=theta,
                         p=p, lambda=lambda, omega=omega, reps=reps)
  
  base_s <- base[,seq(from=1, to=ncol(base), by=2)]
  base_g <- base[,seq(from=2, to=ncol(base), by=2)]
  
  
  for(i in 1:ncol(base_s)){
    boot <- bootstrap(base_g[,i], base_s[,i], gather.overlap=TRUE)
    s_over[i] <- boot[1]
    g_over[i] <- boot[2]
  }
  
  return(cbind(s_over, g_over))
}


#' Standalone function to be used with the output of runChanging
#' TODO: make it generate a ggplot with confidence intervals
#' Returns a dataframe easy to use for plotting, with each row representing
#' a separate found value for the current changing_value, lambda, and omega.
plotEquilibriumSimulation <- function(fitness) {
  
  plotting <- c()
  changing_value <- as.factor(fitness$changing_value)
  lambda <- as.factor(fitness$lambda)
  omega <- as.factor(fitness$omega)
  
  
  changing_value_levels <- levels(changing_value)
  lambda_levels <- levels(lambda)
  omega_levels <- levels(omega)
  
  for(curr in changing_value_levels){
    for(l in lambda_levels){
      for(o in omega_levels){
        
        plotting <- rbind(plotting, fitness[which(fitness$changing_value == curr &
                                                    fitness$lambda == l &
                                                    fitness$omega == o)[1],])
      }
    }
  }
  
  # have this generate a plot and not just return
  
  return(plotting)
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
