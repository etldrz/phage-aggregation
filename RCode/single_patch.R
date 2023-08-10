burst.size.A <- 50 # burst size of host A
burst.sizes.B <- 0:(burst.size.A + as.integer(burst.size.A * 0.05))
reps <- 1.5e6 # how large each fitness vector is
bt.size <- 1e4 # the size of the bootstrapped fitness vector
allowed.overlap <- 0.01 # proportion of allowed overlap between bootstrapped vectors



#' #' Helper function used throughout
#' #' Used to verify floating numbers
#' floatMatch <- function(x, y) { abs(x - y) < 1e-6 }


viableParams <- function(alpha, theta, lambda) {
  worst.case <- (5^theta) / (theta + lambda + alpha) * exp(-10*lambda)
  
  best.case <- (45^theta) / (theta + lambda + alpha) * exp(-0.1*lambda)
  cbind(worst.case, best.case)
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
  
  for(b in burst.sizes.B){
    ws.fitness <- sample(c(0, 1, 2), reps, replace=TRUE, prob=c(lambda, alpha, theta*p))
    wg.fitness <- sample(c(0, 1, 2, 3), reps, replace=TRUE, prob=c(lambda, alpha, theta*p, theta*(1-p)))
    
    ws.fitness <- findLyseFitness(ws.fitness, burst.size.B=b, alpha, theta, p, 
                                  lambda, omega, is.specialist=TRUE)
    wg.fitness <- findLyseFitness(wg.fitness, burst.size.B=b, alpha, theta, p,
                                  lambda, omega, is.specialist=FALSE)
    
    base <- cbind(base, ws.fitness, wg.fitness)
  }
  
  return(base)   
}


#' Helper function used by baseSimulation
#' This will be called to deal with phage interactions with hosts, in order to
#' calculate the fitnesses of children and grandchildren.
findLyseFitness <- function(outcomes, burst.size.B, alpha, theta, p, lambda, omega, 
                            is.specialist) {
  # infected.bacteria represents a phage interacting with a host
  infected.bacteria <- which(outcomes %in% c(2, 3))
  if(length(infected.bacteria) == 0) return(outcomes)
  
  # logical vector indicating whether the infected host bursts inside the aggregation
  burst.chances <- rbinom(infected.bacteria, 1, omega)

  # The loop covers both specialist and generalist by dealing with
  # outcomes having 2 and 3 via an if else statement. The fitnesses of
  # children and grandchildren are then found and put into outcomes, which is 
  # then returned.
  for(k in 1:length(infected.bacteria)){
    loc <- infected.bacteria[k]
    
    if(outcomes[loc] == 2){
      burst.size <- burst.size.A
      if(burst.chances[k] == TRUE){
        inside.lyse.fitness <- NULL
        
        if(is.specialist){
          inside.lyse.fitness <- sample(c(0, 1, burst.size.A), burst.size, replace=TRUE,
                                        prob=c(lambda, alpha, theta*p))
        }else{
          inside.lyse.fitness <- sample(c(0, 1, burst.size.A, burst.size.B), burst.size, replace=TRUE,
                                        prob=c(lambda, alpha, theta*p, theta*(1-p)))
        }
        
        outcomes[loc] <- sum(inside.lyse.fitness)
      }else{ # if the bacterium bursts outside the aggregation
        outcomes[loc] <- burst.size
      }
    }else if(outcomes[loc] == 3){
      burst.size <- burst.size.B

      if(burst.chances[k] == TRUE){
        inside.lyse.fitness <- sample(c(0, 1, burst.size.A, burst.size.B), burst.size, replace=TRUE, 
                                      prob=c(lambda, alpha, theta*p, theta*(1-p)))
        outcomes[loc] <- sum(inside.lyse.fitness)
      }else{
        outcomes[loc] <- burst.size
      }
    }
  }
  return(outcomes)
}


#' Generates a pair of bootstrapped fitness vectors
basicBootstrap <- function(s, g){
  
  size <- length(g) #length of g == length of s
  
  df.g <- as.data.frame(table(g))
  unique.g <- as.numeric(levels(df.g[,1]))
  freq.g <- df.g[,2] / size
  
  df.s <- as.data.frame(table(s))
  unique.s <- as.numeric(levels(df.s[,1]))
  freq.s <- df.s[,2] / size
  
  g.means <- replicate(bt.size, sum(rmultinom(1, size, freq.g) * unique.g) / size)
  s.means <- replicate(bt.size, sum(rmultinom(1, size, freq.s) * unique.s) / size)
  
  # g.means <- replicate(bt.size, mean(sample(unique.g, size, replace=T, prob=freq.g)))
  # s.means <- replicate(bt.size, mean(sample(unique.s, size, replace=T, prob=freq.s)))
  
  return(cbind(s.means, g.means))
}


#' Returns a list with length(changing) sublists, each one of those sublists 
#' contains 13 entries
#' The entry names are:
#'  lower.bound: the burst.sizes.B point where the lower boundary was found
#'  lower.boot.s: the bootstrapped fitness vector for S at the lower bound
#'  lower.boot.g: same as above except for G
#'  upper.bound: the burst.sizes.B point where the upper boundary was found.
#'  upper.boot.s: the bootstrapped fitness vector for S at the upper bound
#'  upper.boot.g: same as above except for G
#'  r.star: the values of the burst.sizes.B point where WG = WS. 
#'    Solved for using the linear slope equation found in equalityPoint
#'  r.star.mean: the mean of r.star
#'  changing.value: the current value of the changing variable
#'  lower.quantile: lower bound of a 95% confidence interval for the r.star vector
#'  upper.quantile: upper bound of a 95% confidence interval for the r.star vector
#'  swap.count.lower: how many times the swapping algorithm was performed on the 
#'    lower point.
#'  swap.count.upper: how many times the swapping algorithm was performed on the 
#'    upper point.
#'  The name of the current changing variable will be stored in the main list's
#'    attributes.
preprocessed <- function(files, changing, changing.name) {
  
  exist <- sapply(files, file.exists)
  if(sum(exist) != length(files))
    stop(cat("the following files do not exist",
             files[!1:length(files) %in% exist],
             sep="\n"))
  
  if(length(files) != length(changing))
    stop("mismatch between files length and changing length")
  
  fitness <- vector("list", length=length(changing))
  for(f in 1:length(files)){
    fitness[[f]] <- rStar(files[f], changing[f])
  }
  
  attributes(fitness)$changing.name <- changing.name
  
  return(fitness)
}


#' Finds R* and all accompanying data for a single file and returns it as a list
#' to preprocessed().
rStar <- function(file, current.changing) {
  boot <- read.table(file, header=TRUE, sep=",")
  
  # Files are organized by repeating the sequence W.S W.G
  s <- boot[,seq(1, ncol(boot), 2)]
  g <- boot[,seq(2, ncol(boot), 2)]
  
  # Containers for possibly useful upper/lower bounds in the file from which
  # R* can be found
  viable.lower <- c()
  viable.upper <- c()
  
  # For each column in the pair of matrices, find the proportion of overlap and 
  # add it to the appropriate container if it is less than or equal to allowed.overlap
  for(i in 1:ncol(g)){
    s.greater <- which(s[,i] > g[,i])
    g.greater <- which(g[,i] > s[,i])
    
    if(1 - (length(s.greater) / nrow(boot)) <= allowed.overlap){
      viable.lower <- c(viable.lower, i)
    }
    else if(1 - (length(g.greater) / nrow(boot)) <= allowed.overlap){
      viable.upper <- c(viable.upper, i)
    }
  }
  
  # Upper and lower points used to calculate R* via the linear slope equation
  lower.pair <- NULL
  upper.pair <- NULL
  lower.burst.B <- NULL
  upper.burst.B <- NULL
  r.stars <- c()
  
  # Finding the lower bound
  while(is.null(lower.pair)){
    # This if block contains placeholder logic for what to do if the simulation 
    # cannot find a viable lower point.
    if(length(viable.lower) == 0){
      message(paste("No lower found for ", file, "\n",
                    "Setting R* equal to 0", sep=""))
      r.stars <- 0
      lower.burst.B <- -Inf
      break
    }
    # Taking the max so it will be closer to viable.upper's chosen point
    lower.burst.B <- max(viable.lower) 
    
    # swap() will return the matrix unaltered if there is nothing to swap and null 
    # if the algorithm fails
    lower.pair <- swap(cbind(s[,lower.burst.B], g[,lower.burst.B]))
    
    viable.lower <- viable.lower[!viable.lower %in% lower.burst.B]
  }
  
  # Finding the upper bound
  while(is.null(upper.pair)){
    if(length(viable.upper) == 0){
      message(paste("No upper found for ", file, "\n",
                    "Setting R* equal to 1", sep=""))
      r.stars <- 1
      upper.burst.B <- Inf
      break
    }
    upper.burst.B <- min(viable.upper)
    upper.pair <- swap(cbind(g[,upper.burst.B], s[,upper.burst.B]))
    viable.upper <- viable.upper[!viable.upper %in% upper.burst.B]
  }
  
  # equalityPoint() finds the point where the slope of W.G equals W.S
  if(is.null(r.stars)) r.stars <- mapply(equalityPoint, 
                                         lower.burst.B, lower.pair[,2], 
                                         lower.pair[,1], upper.burst.B, 
                                         upper.pair[,1], upper.pair[,2])
  
  quant <- quantile(r.stars, probs=c(0.025, 0.975))
  lower.quantile <- quant[1]
  upper.quantile <- quant[2]
  
  data <- list(lower.bound = lower.burst.B, lower.boot.s = lower.pair[,1], 
               lower.boot.g = lower.pair[,2], upper.bound = upper.burst.B,
               upper.boot.s = upper.pair[,2], upper.boot.g = upper.pair[,1], 
               r.star = r.stars, r.star.mean = mean(r.stars), 
               changing.value = current.changing, lower.quantile = lower.quantile, 
               upper.quantile = upper.quantile, 
               swap.count.lower = attributes(lower.pair)$swap.count, 
               swap.count.upper = attributes(upper.pair)$swap.count)
  return(data)
}


#' This attempts to swap overlapping locations with non-overlapping locations 
#' from the matrix so that there is no overlap at all. The first column is 
#' expected to be larger than the second; locations where the second is larger
#' than the first is considered to be overlap.
swap <- function(data){
  
  # Vector of locations inside data where the second column is > than the first
  to.fix <- which(data[,1] <= data[,2])
  to.fix.size <- length(to.fix)
  
  if(to.fix.size == 0) return(data)
  
  
  for(index in to.fix){
    # The current data location trying to be fixed.
    overlap <- data[index,]
    
    # Locations that could potentially swap with the current overlap row
    viable <- which(!1:nrow(data) %in% to.fix)
    
    # Records how many times a redraw happens
    count <- 1 
    
    while(length(viable) > 0){
      
      rand <- sample(viable, 1)
      
      current.viable <- data[rand,]
      
      if(current.viable[1] > overlap[2] & overlap[1] > current.viable[2]){
        
        # Swapping the first column of the two rows
        temp <- current.viable[1]
        data[rand,1] <- overlap[1]
        data[to.fix[index],1] <- temp
        
        to.fix.size <- to.fix.size - 1
        
        if(to.fix.size > 0){
          break
        }
        attributes(data)$swap.count <- count
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


#' Helper function used by rStar
#' Finds the point where wG = wS and return that point by solving the 
#' linear slope equations for both.
equalityPoint <- function(min.x, min.wg, min.ws, max.x, max.wg, max.ws) {
  
  g.slope <- (max.wg - min.wg) / (max.x - min.x)
  s.slope <- (max.ws - min.ws) / (max.x - min.x)
  
  g.intercept <- min.wg - (g.slope * min.x)
  s.intercept <- min.ws - (s.slope * min.x)
  
  
  # Now solving g.slope*x + g.intercept = s.slope*x + s.intercept
  x <- (s.intercept - g.intercept) / (g.slope - s.slope)
  
  return((x - 1) / (burst.size.A - 1)) # Dividing nB by nA will return R*
}


#' Standalone function
#' Generates a matrix with two columns and burst.size.A rows which contains the 
#' predictions for W_s and W_g of the base simulation.
baseSimPrediction <- function(alpha, theta, p, lambda, omega) {
  
  ws.prediction <- alpha/(alpha + lambda + theta*p) + 
    theta*p/(alpha + lambda + theta*p)*((1 - omega)*(burst.size.A - 1) + 
                                          omega*((burst.size.A - 1) + 1)*(alpha/(alpha + lambda + theta*p) + 
                                                                            (burst.size.A - 1)*theta*p/(alpha + lambda + theta*p)))
  
  wg.prediction <- alpha/(alpha + lambda + theta) + 
    theta*p/(alpha + lambda + theta)*((1 - omega)*(burst.size.A - 1) + 
                                        omega*((burst.size.A - 1) + 1)*(alpha/(alpha + lambda + theta) + 
                                                                          theta/(alpha + lambda + theta) * (p * (burst.size.A - 1) + (1 - p) * (burst.sizes.B- 1)))) +
    theta*(1 - p)/(alpha + lambda + theta)*((1 - omega)*(burst.sizes.B - 1) + 
                                              omega*((burst.sizes.B - 1) + 1)*(alpha/(alpha + lambda + theta) + 
                                                                                theta/(alpha + lambda + theta) * (p * (burst.size.A - 1) + (1 - p) * (burst.sizes.B- 1))))
  data <- cbind(ws.prediction, wg.prediction)
  
  # plot(x=burst.sizes.B, y=data[,2], type='l', col='firebrick', lwd=1.5, ylab="fitness")
  # lines(x=burst.sizes.B, y=data[,1], col='darkblue', lwd=1.5)
  # legend("topleft", legend=c("fitness.s", "fitness.g"), lty=1, 
  #        col=c('darkblue', 'firebrick'), lwd=1.5)
  return(data)
}


complexSimPrediction <- function(alpha, theta, p, lambda){ 1 / (1 + (alpha + lambda) / (theta * p)) }


#' Standalone function that plots the results generated by baseSimulation
#' and baseSimPrediction
plotBaseSimulation <- function(base, colors=c('firebrick', 'darkorchid4')) {
  current.ws <- colMeans(base[,seq(from=1, to=ncol(base), by=2)])
  current.wg <- colMeans(base[,seq(from=2, to=ncol(base), by=2)])
  
  y.min <- min(current.wg)
  y.max <- max(current.wg)
  if(max(current.ws) > max(current.wg))
    y.max <- max(current.ws)
  
  
  # total <- data.frame(rbind(current.wg, current.ws))
  # total$x <- burst.sizes.B
  # names(total) <- c("fitness", "type", "x")
  # suppressMessages(library('ggplot2'))
  # plot <- ggplot(total, aes(x=x, y=fitness, color=type, group=type)) +
  #   geom_line() +
  #   scale_y_continuous(n.breaks=8) +
  #   theme_classic()
  # return(plot)
  
  plot(y=current.ws, x=1:length(current.ws), type='l', col=colors[1], 
       xlim=c(1, length(current.ws)), ylim=c(y.min, y.max), lwd=1.5, ylab="fitness",
       xlab="burst.sizes.B", cex.lab = .75)
  lines(y=current.wg, x=1:length(current.wg), col=colors[2], lwd=1.5)
  #legend('topleft', legend=c("WS","WG"), fill=colors)
}


#' When given two fitness matrices, this function will return a matrix containing
#' values comparing overlap for each ith column.
checkOverlap <- function(s, g) {
  out <- mapply(function(x, y) {
    out <- matrix(ncol=2, nrow=1)
    colnames(out) <- c("S.greater", "G.greater")
    
    out[1,1] <- length(which(x > y))
    out[1,2] <- bt.size - out[1,1]
    if(out[1,1] + out[1,2] != bt.size)
      message("at least one pair of values are equal")
    return(out)
  }, 
  s, g, SIMPLIFY=FALSE)
  
  return(do.call(rbind, out))
}

#' Takes a fitness list generated by preprocessed() and returns a ggplot object
plotFitness <- function(fitness, generate.plot=TRUE) {
  plotting <- data.frame(matrix(ncol=4, nrow=length(fitness)))
  names(plotting) <- c("r.star.mean", "lower.quantile",
                       "upper.quantile", attributes(fitness)$changing.name)
  
  for(i in 1:length(fitness)){
    plotting[i,] <- c(fitness[[i]]$r.star.mean, fitness[[i]]$lower.quantile,
                      fitness[[i]]$upper.quantile, fitness[[i]]$changing.value)
  }
  
  if(!generate.plot) return(plotting)
  
  suppressMessages(library(ggplot2))
  plot <- ggplot(plotting, aes(x=plotting[,4], y=plotting[,1])) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=plotting[,2], ymax=plotting[,3]), width=0.01, alpha=.75) +
    xlab(attributes(fitness)$changing.name) +
    ylab("R*") +
    theme_classic()
  unloadNamespace('ggplot2')
  
  
  return(plot)
}


#' Technical overview
#' Parameters are first tested to see if they are biologically plausible via 
#' the function viableParams. Then, a matrix of non-altered fitness data is
#' is created by baseSimulation and then saved to a txt file. This fitness data
#' is bootstrapped and saved to its own respective txt file. Once a desired
#' number and set of parameters have been chosen and bootstrapped, the function
#' preprocessed iterates through each file and calls rStar on it. rStar finds 
#' the equilibrium point by using the linear slope equation for each matrix 
#' where W.S equals W.G and returns that and other relevant data, such as 95%
#' quantiles and the upper and lower bounds used for the slope equation. 
#' preprocessed will return a list of the found values for each file.
#' 
#' rStar finds viable points for the linear slope equation by noting B burst 
#' sizes (columns) where all of the generalist fitness data are greater
#' than the specialist fitness (upper bound) and where all of the specialist
#' fitness data are greater than the generalist fitness (lower bound). If there
#' is a small amount of incorrect overlap, less than or equal to 1 percent, then 
#' an algorithm which randomly swaps rows of the current pair of bootstrapped
#' fitness vectors is called to try and reduce the amount of incorrrect overlap 
#' to 0%. 
#' 
#' This swapping algorithm is called for two reasons. Firstly, to try and
#' decrease the amount of times that a suitable upper or lower bound is not
#' found for the linear slope equation. Secondly, to potentially offset failures
#' of finding the actual R* value due to W.G having a nonlinear slope--the 
#' closer that the upper and lower bounds are together, the greater the
#' likelihood of the line segment being linear. Should an upper or lower bound
#' not be found, the R* value is set to be 1 or 0 respectively and a note is 
#' made on the plot.
