
# TODO: reformat this comment block.
#       It will until the the fourth changing point (alpha = 0.07) at which point it will fail to collect any rows from the fitness df 
#       send to zeroHunter.



###The main function, runChanging, calls equilibriumSimulation 
###length(changing_variable)*length(omega)*length(lambda) times. Each 
###time it is called, it finds the highest point where every entry of a bootstrapped 
###fitness vector of the generalist is higher than each entry of a bootstrapped 
###fitness vector of the specialist. The simulation then reverses and finds the lowest 
###point where each entry of a bootstrapped fitness vector of the specialist is 
###higher than that of the generalist. These four fitness vectors are recorded. 
###The positions of the found vectors (their x values) are recorded as well. 
###Once this loop finishes, the equilibrium point is found by determining 
###the slopes between both pairs of fitness vectors and solving for when the 
###slope of the specialist equals the 
###slope of the generalist. The mean and a 95% confidence interval is taken and 
###recorded for each value in alpha_change, this is what is plotted, along with
###the particular value's mathematical prediction.


# Global error messages:
no_upper_found <- "Reached the end of total_treatments without finding an upper limit"
no_lower_found <- "Reached the bottom without finding a lower limit"


#' Helper function used throughout
#' For verifying floating numbers
floatMatch <- function(x, y){
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
bootstrap <- function(fitness_g, fitness_s, bt_size=10000, reverse=FALSE, gather.overlap=FALSE) {
  
  g_means <- rep(0, bt_size)
  s_means <- rep(0, bt_size)
  
  # indicates the amount of overlapping that takes place; used when gather.overlap = T
  g_overlap <- 0
  s_overlap <- 0
  
  
  for (i in 1:bt_size) {
    g_means[i] <- mean(sample(fitness_g, length(fitness_g), replace=TRUE))
    s_means[i] <- mean(sample(fitness_s, length(fitness_s), replace=TRUE))
    
    if(!gather.overlap){
      if(!reverse){
        if(g_means[i] < s_means[i]){
          return(FALSE)
        }
      } 
      else if(reverse){
        if(s_means[i] < g_means[i]){
          return(FALSE)
        }
      }
    }
    else if(gather.overlap){
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
  else
    return(cbind(g_means, s_means))
}


#' means[,1] is the vector that is supposed to be larger than means[,2]
#' NOT CURRENTLY RECORDING ANYTHING
switcheroo <- function(means, perc_cap){
  perc <- 0
  if(any(means[,1] < means[,2]))
    perc <- (length(which(means[,1] < means[,2])) / nrow(means))*100
  
  if(perc < perc_cap | floatMatch(perc, perc_cap))
    return(means)
  
  
  
  
  to_fix <- which(means[,1] < means[,2] | floatMatch(means[,1], means[,2]))

  legal.perc <- FALSE
  
  for(index in to_fix){
    if(legal.perc)
      return(means)
    
    viable <- which(!1:nrow(means) %in% to_fix)

    overlap <- means[index,]
    
    count <- 1 # records how many times a redraw happens (how many items from viable were used)
    while(!legal.perc & length(viable) > 0){
      
      rand <- sample(viable, 1)
      means_rand <- means[rand,]
      
      if((means_rand[1] > overlap[2] | floatMatch(means_rand[1], overlap[2])) &
         (overlap[1] > means_rand[2] | floatMatch(overlap[1], means_rand[2]))){
        
        temp <- means[rand,1]
        means[rand,1] <- overlap[1]
        means[index,1] <- temp
        
        
        to_fix <- to_fix[which(!to_fix %in% index)]
        
        # If count ever gets used, then it will get wrongfully incremented if this chunk triggers
      }
      
    
      # Update viable to exclude the just-used value
      viable <- viable[which(!viable %in% rand)]
      

      curr_failure_perc <- (length(to_fix) / nrow(means)) * 100
      legal.perc <- curr_failure_perc < perc_cap | floatMatch(curr_failure_perc, perc_cap)
      
      count <- count + 1
    }
  }
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
slopeEqual <- function(min_x, min_wg, min_ws, max_x, max_wg, max_ws){
  
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
zeroHunter <- function(fitness){
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


#' Helper function used in runChanging
#' This runs the simulation and finds all the necessary data to find R* and the 
#' bootstrapped confidence intervals. base_outcomes are two matrices which hold nA columns,
#' one for each point in the base plot.
#' changing: a string of the current changing variable name
#' bt_size: length of the bootstrap vector
equilibriumSimulation <- function(plot_file, nA, omega, alpha, theta, p, lambda, 
                                  fitness_s, fitness_g, changing, bt_size) {
  
  # The names of the dataframe generated by equilibriumSimulaiton
  legend <- c('x_lower', 'bt_ws_lower','bt_wg_lower','x_upper', 'bt_ws_upper', 
              'bt_wg_upper', 'lower_quantile', 'upper_quantile', 'r_star',
              'r_star_mean', 'changing_value', 'lambda', 'omega', 'p', 'prediction')
  
  
  #' These hold the fitnesses during the length of total_treatments; the fitness vectors
  #' that are not needed by the end are discarded.
  g_outcomes_dataframe <- c()
  s_outcomes_dataframe <- c()
  
  current_value <- theta
  if(changing == "alpha") current_value <- alpha
  
  nB_vec <- 1:nA
  
  total_treatments <- length(nB_vec)
  
  
  # When this is true, the parent loop will know to break. It will be true once
  # all four fitness vectors are found.
  halt <- FALSE
  
  for(i in nB_vec){
    if(halt)
      break
    
    
    # The nB value to be used in this loop
    nB <- nB_vec[i]
    
    
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
    

    #' This will check if the matrix is viable by seeing if each found bootstraped
    #' value for fitness_g is either larger or lower than the found bootstrapped 
    #' value for fitness_s, depending on if reverse is false or true.
    #' If the given fitness vectors are not viable, then matrix_or_bool will be
    #' false, if they are viable then it will be a matrix with two columns
    #' representing the bootstrapped values of the inputted fitness_g and 
    #' fitness_s with g,s being the ordering of the columns.
    matrix_or_bool <- bootstrap(fitness_g[,i], fitness_s[,i], bt_size=bt_size, 
                                reverse=FALSE)

    
    if(isFALSE(matrix_or_bool)){
      
      # If fitness_g and fitness_s are not viable, then they will be saved
      # so that they can be accessed when reverse=TRUE
      g_outcomes_dataframe <- cbind(g_outcomes_dataframe, fitness_g[,i])
      s_outcomes_dataframe <- cbind(s_outcomes_dataframe, fitness_s[,i])
      
      # This check just makes sure that the user is notified when the simulation
      # doesn't find a value when it should.
      if(i < total_treatments)
        next
      else 
        stop(no_upper_found, ", ", current_value)
      
     } 
    else if(is.matrix(matrix_or_bool)){
      print(paste("Uppermost found at", i, "out of", total_treatments))
      
      # Using bound to record the bootstrapped viable upper-limit fitness vectors,
      # as well as the nB they were found at. Dividing by nA so R* will be found.
      bound[,4:6] <- cbind(rep(nB, bt_size), 
                           matrix_or_bool[,2], # bootstrapped specialist
                           matrix_or_bool[,1]) # bootstrapped generalist
      
      for(y in (nB - 1):1) {
        
        matrix_or_bool <- bootstrap(g_outcomes_dataframe[,y],
                                    s_outcomes_dataframe[,y], bt_size=bt_size,
                                    reverse=TRUE)
        
        # Will continue reversing if bootstrap() says the vectors aren't viable by
        # returning a boolean
        if(isFALSE(matrix_or_bool) & y == 1){
          stop(no_lower_found, ", ", current_value)
        }
        else if(isFALSE(matrix_or_bool) & y > 1){
          next
        }

        # Place the lower entries into the bound vector
        bound[,1:3] <- cbind(rep(y, bt_size), #y is the found lower value of nB
                             matrix_or_bool[,2], matrix_or_bool[,1])
        
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
baseSimulation <- function(nA, alpha, lambda, omega, theta, p, reps=500000){
  # Simulated fitnesses
  base <- c()

  
  for(i in 0:nA){
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
baseSimPrediction <- function(nA, alpha, theta, p, lambda, omega){
  
  nB <- 1:nA
  
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
plotBaseSimulation <- function(base, prediction=NA, with.prediction=FALSE, colors=c('coral4', 'darkorchid4')){
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
runChanging <- function(plot_file, changing_name, changing_step_size, change_start, changing_length=10,
                         nA, alpha, theta, p, lambda, omega, save_base=FALSE, 
                        base_file="", bt_size=10000, reps=500000, gather.overlap=FALSE) {
  
  if(save_base & !file.exists(base_file))
    stop("Use a proper file path for the base")
  
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
        for(per in p){
          current_base <- baseSimulation(nA=nA, alpha=alpha, lambda=l,
                                         omega=o, theta=theta, p=per, reps=reps)
          
          current_base_wS <- current_base[,seq(from=1, to=ncol(current_base), by=2)]
          current_base_wG <- current_base[,seq(from=2, to=ncol(current_base), by=2)]
          
          base <- cbind(base, current_base)
          equilibriumSimulation(plot_file=plot_file, nA=nA, omega=o, alpha=alpha,
                                theta=theta, p=per, lambda=l, fitness_s=current_base_wS,
                                fitness_g=current_base_wG, changing=changing_name, bt_size=bt_size)
        }
      }
    }
    
    # Currently, the base is not well organized and would be silly to save.
    if(save_base)
      write.table(base, file=base_file)
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


safeRunChanging <- function(safety=1.5e6, reps=5e5, ...){
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
                          bt_size=10000, reps=500000){

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
plotEquilibriumSimulation <- function(fitness){

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

#' 
#' #' Standalone function
#' #' Makes a grid of of the two parameters which are left NA. Throws an error if two parameters
#' #' aren't NA.
#' #' final_file: where the grid is stored using write.table
#' #' reps: the size of the fitness vectors which baseSimulation creates
#' #' seq_step_size: the width and height of each square in the grid
#' #' margin: the LOWER significant percentage value. mean(WS)/nA less than this is significant
#' #'         for WG and mean(WS)/nA greater than 1-this is significant for WS.
#' #' colors: a character vector of length 4, with [1] being not-significant either way
#' #'         [2] being WG significant, [3] being WS significant, and [4] being very WS
#' #'         significant (>1)
#' quantileGrid <- function(final_file, nA, alpha=NA, theta=NA, p=NA, lambda=NA, omega=NA,
#'                         reps=500000, seq_step_size, colors){
#' 
#'   if(!file.exists(final_file))
#'     stop("The file path is not valid.")
#' 
#'   if(length(colors) != 4)
#'     stop("Incorrent length of your colors vector")
#' 
#'   if(sum(is.na(c(nA, alpha, theta, p, lambda, omega))) != 2)
#'     stop("Too few or too many axes chosen")
#' 
#' 
#' 
#'   params <- c(nA, alpha, theta, p, lambda, omega)
#'   var_names <- c('nA','alpha','theta','p','lambda','omega')
#'   names(params) <- var_names
#' 
#' 
#'   xs <- seq(0, 1, seq_step_size)
#'   ys <- seq(0, 1, seq_step_size)
#' 
#' 
#'   grid <- matrix(colors[1], nrow=length(ys), ncol=length(xs),
#'                  dimnames=list(as.character(ys), as.character(xs)))
#' 
#'   # Used for plotting labels
#'   attributes(grid)$x_name <- var_names[which(is.na(params))[1]]
#'   attributes(grid)$y_name <- var_names[which(is.na(params))[2]]
#'   attributes(grid)$held_steady <- var_names[which(!is.na(params))]
#'   for(x_val in xs){
#' 
#'     params[[attributes(grid)$x_name]] <- x_val
#' 
#'     for(y_val in ys){
#' 
#'       params[[attributes(grid)$y_name]] <- y_val
#' 
#'       base <- baseSimulation(params[['nA']], params[['alpha']], params[['lambda']],
#'                              params[['omega']], params[['theta']], params[['p']], reps)
#' 
#'       mean_wS <- mean(base[,seq(from=1, to=ncol(base), by=2)])
#'       mean_wG <- mean(base[,seq(from=2, to=ncol(base), by=2)])
#' 
#'       perc <- mean_wS / mean_wG
#' 
#' 
#'       if(is.na(perc))
#'         next
#'       else if(perc > 1)
#'         grid <- fillLocation(grid, x_val, y_val, colors, FALSE)
#'       else if(perc < 1)
#'         grid <- fillLocation(grid, x_val, y_val, colors, TRUE)
#'       else if(floatMatch(perc, 1))
#'         grid <- fillLocation(grid, x_val, y_val, colors, NULL)
#'     }
#'   }
#' 
#'   # #for reading comments
#'   # content <- readLines(final_file)
#'   # comments <- content[grep("^#", content)]
#'   # #
#' 
#'   write.table(grid, file=final_file, append=FALSE, quote=FALSE, sep=",")
#' 
#'   return(grid)
#' } #end of quantileGrid
#' 
#' 
#' #' Helper function for quantileGrid
#' #' Fills in the appropriate location on the grid with the appropriate color.
#' #' TODO: fix the inclusive/exclusive logic
#' fillLocation <- function(grid, x_val, y_val, colors, is.lower){
#'   x_loc <- as.numeric(rownames(grid))
#'   y_loc <- as.numeric(colnames(grid))
#' 
#' 
#'   x <- which(sapply(x_loc, floatMatch, x_val))
#'   y <- which(sapply(y_loc, floatMatch, y_val))
#' 
#' 
#'   grid[y,x] <- colors[2] # < 1
#'   if(is.null(is.lower))
#'     grid[y,x] <- colors[4] # > 1
#'   else if(is.lower)
#'     grid[y,x] <- colors[3] # = 1
#' 
#'   return(grid)
#' }
#' 
#' 
#' #' Standalone function to be used with the output of quantileGrid.
#' #' Plots the given matrix as a visual grid.
#' plotQuantileGrid <- function(grid, colors, nA, alpha, theta, p, lambda, omega){
#' 
#' 
#'   rows <- as.numeric(rownames(grid))
#'   cols <- as.numeric(colnames(grid))
#' 
#'   dimension <- cols[2] - cols[1]
#' 
#'   plot(cols~rows, xlab=attributes(grid)$x_name,
#'        ylab=attributes(grid)$y_name, pch="",
#'        ylim=c(0, 1 + dimension), xlim=c(0, 1 + dimension))
#' 
#'   par(mar=c(5, 4, 4, 10), xpd=TRUE)
#'   legend("topright", inset=c(-0.5, 0),
#'          legend=c("WG prominent",
#'                   "WS prominent", "Equal mean"),
#'          fill=colors[2:4])
#' 
#'   y <- 1
#'   for(row in rows){
#'     x <- 1
#'     for(col in cols){
#'       x1 <- col
#'       x2 <- col + dimension
#'       y1 <- row
#'       y2 <- row + dimension
#'       rect(x1,y1,x2,y2, col=grid[y,x])
#'       x <- x+1
#'     }
#'     y <- y+1
#'   }
#' 
#' 
#'   params <- c(nA, alpha, theta, p, lambda, omega)
#'   var_names <- c('nA','alpha','theta','p','lambda','omega')
#'   names(params) <- var_names
#' 
#' 
#' 
#'   title_string <- ""
#' 
#'   for(steady in attributes(grid)$held_steady)
#'     title_string <- paste(title_string, steady, "=", params[[steady]], ";")
#' 
#'   title_string <- substr(title_string, 1, nchar(title_string) - 1)
#' 
#'   title(main=title_string, cex.main=0.9)
#' }
