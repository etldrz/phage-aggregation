###This simulation plots the found equilibrium for each changing value of alpha against the changing value
###of alpha. It runs through two different values of lambda, records all the necessary data, and plots 
###both lines on the same graph.

###The main loop calls equilibriumSimulation() length(alpha_change) times. Each time it is
###called, it finds the lowest point where each entry of a bootstrapped fitness vector of the generalist
###is higher than each entry of a bootstrapped fitness of the specialist. 
###The simulation then reverses and finds the uppermost point where each entry of a bootstrapped
###fitness vector of the specialist is higher than that of the generalist. These 
###four fitness vectors are recorded. The positions of the found vectors (their 
###x values) are recorded as well. This process is done for two different lambdas,
###so 2*length(alpha_change) times. Once this loop finishes, the equilibrium point
###is found by determining the slopes of each and solving for when the slope of the
###specialist equals the slope of the generalist. The mean and a 95% confidence interval
###is taken and recorded for each value in alpha_change, these are what are plotted, along
###with the mathematical prediction that is made for each changing value.

###The dataframe columns are:
###'x_lower', 'bt_ws_lower','bt_wg_lower','x_upper', 'bt_ws_upper', 'bt_wg_upper', 
###''lower_quantile', 'upper_quantile', 'bootsrapped_r_star', 'r_star_mean', 
###'changing_value', 'lambda', 'prediction')

library(ggplot2)
plot_file <- "C:\\Users\\Evan\\Desktop\\R Projects\\phageAggregation\\plottingData\\alphaChangePlot.txt"
#Where the dataframe is stored:

bootstrap <- function(vec) {
  bt_size <- 10000
  bt_means <- rep(0, bt_size)
  for (i in 1:bt_size) {
    bt_means[i] <- mean(sample(vec, length(vec), replace=TRUE))
  }
  return(bt_means)
}

#this is used to find R* by solving the equation m*wG + b = wS
##The start of the vector is the first location where wS is bigger than 
  ##wG and the end is the first location where wG is bigger than wS
zeroHunter <- function(fitness){

  slopeEqual <- function(min_x, min_wg, min_ws, max_x, max_wg, max_ws){
    
    g_slope <- (max_wg - min_wg) / (max_x - min_x)
    s_slope <- (max_ws - min_ws) / (max_x - min_x)

    #To catch dividing by zero
    if(!is.integer(s_slope)){
      s_slope <- 0
    } else if(!is.integer(g_slope)){
      g_slope <- 0
    }
    
    
    intercept_g <- min_wg - (g_slope * min_x)
    intercept_s <- min_ws - (s_slope * min_x)

    
    #Now solving g_slope*x + intercept_g = s_slope*x + intercept_s
    x <- (intercept_s - intercept_g) / (g_slope - s_slope)
    if(x < 0){
      stop("slopeEqual just gave negative rStar")
    }
    return(x)
  }
  
  initial_r_star <- rep(0, nrow(fitness))
  
  for(i in 1:nrow(fitness)){
    initial_r_star[i] <- slopeEqual(min_x = fitness$x_lower[i], min_wg = fitness$bt_wg_lower[i], 
                                 min_ws = fitness$bt_ws_lower[i], max_x = fitness$x_upper[i], 
                                 max_wg = fitness$bt_wg_upper[i], max_ws = fitness$bt_ws_upper[i])
  }
  
  bt_r_star <- bootstrap(initial_r_star)
  quant <- quantile(bt_r_star, probs=c(0.025, 0.975))
  
  
  lower_quantile <- rep(quant[1], nrow(fitness))
  upper_quantile <- rep(quant[2], nrow(fitness))
  bootstrapped_r_star <- bt_r_star
  r_star_mean <- rep(mean(bt_r_star), nrow(fitness))
  combined <- cbind(lower_quantile, upper_quantile, bootstrapped_r_star, r_star_mean)
  
  return(combined)
}



#this runs the simulation and finds R*, as well as the bootstrapped confidence intervals
equilibriumSimulation <- function(nA, lyse, alpha, theta, p,
                                  lambda) {


  replaceFitness <- function(outcomes) {
    # outcomes: 1 = death, 2 = leave, 3 = host A, 4 = host B


      found_death <- which(outcomes == 1)
      found_alpha <- which(outcomes == 2)

      #found_theta represents a phage interacting with a host
      found_theta <- which(outcomes != 1 & outcomes != 2)

      #The loop covers both specialist and generalist by dealing with
        #outcomes having 3 and 4 via an if else statement. The fitnesses of
        #children and grandchildren are then found and put into outcomes.
      for(k in found_theta){
        if(outcomes[k] == 3){
          burst_size <- rpois(1, nA)
          if(runif(1) < lyse){
              outcomes_lyse <- sample(1:3, burst_size, replace=TRUE,
                                      prob=c(lambda, alpha, theta*p))

              tab_lyse <- tabulate(outcomes_lyse, 3)
              outcomes[k] <- sum(tab_lyse*c(0, 1, nA + 1))
          } else {
            outcomes[k] <- burst_size
          }
        } else if(outcomes[k] == 4){
          burst_size <- rpois(1, nB)
          if(runif(1) < lyse){
              outcomes_lyse <- sample(1:3, burst_size, replace=TRUE,
                                      prob=c(lambda, alpha, theta*(1-p)))

              tab_lyse <- tabulate(outcomes_lyse, 3)
              outcomes[k] <- sum(tab_lyse*c(0, 1, nB + 1))
            } else {
            outcomes[k] <- burst_size
          }
        }
      }

      outcomes[found_death] <- 0
      outcomes[found_alpha] <- 1

      bootstrap <- bootstrap(outcomes)
      print(paste("The mean for the found bootstrapped outcomes is ", mean(bootstrap)))
      return(bootstrap)
    }


  #This is used to find the bootstrapped values of wS and wG. When reverse changes signs,
  #then the script knows to stop and save the used S_fitness and G_fitness
  zeroHunterPrecision <- function(s_fitness, g_fitness, reverse=FALSE) {

    #Reverse starts out as false as the first iteration calls it. It remains false
    #until wG has greater fitness than wS. Then reverse is returned true, at which
    #point the script calling this knows to record the values of s_fitness and g_fitness
    #. reverse stays true until it finds the point where all boostrapped values
    #of wS are greater than wG. Then these values of s_fitness and g_fitness are recorded

    if(all(g_fitness > s_fitness)){
      reverse <- TRUE
      return(reverse)
    }
    if(reverse & all(s_fitness > g_fitness)){
      reverse <- FALSE
      return(reverse)
    }

    return(reverse)
  }


  #These hold the fitnesses during the length of treatments; the fitness vectors
    #that are not needed by the end are discarded.
  g_outcomes_dataframe <- c()
  s_outcomes_dataframe <- c()


  nBs <- 1:nA

  treatments <- length(nBs)

  reps <- 500000

  stop <- FALSE

  for(i in 1:treatments){
    #stop is true once all four appropriate fitness vectors are located.
    if(stop){
      break
    }

    nB <- nBs[i]


    s_outcomes <- sample(1:3, reps, replace=TRUE,
                        prob=c(lambda, alpha, theta * p))
    g_outcomes <- sample(1:4, reps, replace=TRUE,
                         prob=c(lambda, alpha, theta * p, theta * (1-p)))
    bootstrap_fitness_g <- replaceFitness(g_outcomes)
    bootstrap_fitness_s <- replaceFitness(s_outcomes)



    #This dataframe is added to so that past outcomes can be accessed by
      #zeroHunterPrecision() when it starts to reverse.
    g_outcomes_dataframe <- cbind(g_outcomes_dataframe, bootstrap_fitness_g)
    s_outcomes_dataframe <- cbind(s_outcomes_dataframe, bootstrap_fitness_s)

    #If b is true, then the bootstrapped vector of wG is larger than the wS vector.
    b <- zeroHunterPrecision(s_outcomes_dataframe[,i], g_outcomes_dataframe[,i])

    if(i == treatments & !b){
      stop("Reached the end of treatments without finding an upper limit.")
    }
    if(b) {
      print(paste("Uppermost found at", i, "in", treatments))

      #Mathematical prediction
      prediction <- nA / (1 + (alpha + lambda) / (theta * p))

      size <- length(s_outcomes_dataframe[,i])
      #This holds all of the relevant data that will be saved and used later.
      bound <- cbind(NA, NA, NA, rep(nB, size),
                     s_outcomes_dataframe[,i],  g_outcomes_dataframe[,i], NA,
                     NA, NA, NA,
                     rep(alpha, size), rep(lambda, size),
                     rep(prediction, size))

      legend <- c('x_lower', 'bt_ws_lower','bt_wg_lower','x_upper',
                  'bt_ws_upper', 'bt_wg_upper', 'lower_quantile',
                  'upper_quantile', 'bootstrapped_r_star', 'r_star_mean',
                  'changing_value', 'lambda', 'prediction')

      for(y in nB:1) {
        #A boolean that states whether bootstrap_fitness_g < bootstrap_fitness_s
        b <- zeroHunterPrecision(s_outcomes_dataframe[,y],
                                 g_outcomes_dataframe[,y], reverse = TRUE)

        if(!b){
          #Place the lower entries into the bound vector
          bound[,1:3] <- cbind(rep(y, length(s_outcomes_dataframe[,y])),
                               s_outcomes_dataframe[,y], g_outcomes_dataframe[,y])

          bound <- as.data.frame(bound)
          names(bound) <- legend
          print(head(bound))
          #These ifs deal with moving the values into a text file
          if(file.info(plot_file)$size != 0){

            fitness_bootstrap <- read.table(file=plot_file, header=TRUE, sep=",")
            names(fitness_bootstrap) <- legend
            fitness_bootstrap <- rbind(fitness_bootstrap, bound)
            write.table(fitness_bootstrap, file=plot_file, append=FALSE, sep=",", quote=FALSE)
            stop <- TRUE
            break
          } else{

            fitness_bootstrap <- bound
            names(fitness_bootstrap) <- legend
            write.table(fitness_bootstrap, file=plot_file, append=FALSE, sep=",", quote=FALSE)
            stop <- TRUE
            break
          }
        }
      }
    }
  }
}


nA <- 30

#alpha is the variable that is changing in this simulation
alpha_change <- 1:20 / 1000
iterations <- length(alpha_change)


lyse <- 0.001
theta <- 0.2
p <- 0.3
lambda1 <- 0.001
lambda2 <- 0.01


### Be sure to take mean of bootstrapped data for new R*
for (i in 1:iterations) {
  #6, lambda2
  print(c("Starting " = i))

  #This finds the data needed to find R* and the bootstrapped confidence intervals
    #for changing values of alpha.
  equilibriumSimulation(nA=nA, lyse=lyse, alpha=alpha_change[i],
                                     theta=theta, p=p, lambda=lambda1)

 equilibriumSimulation(nA=nA, lyse=lyse, alpha=alpha_change[i],
                                    theta=theta, p=p, lambda=lambda2)

}


fitness <- read.table(file=plot_file, sep=",", header=TRUE)

#This finds the equilibrium point and finishes filling the dataframe with r_star,
  #lower_quantile, and upper_quantile.
for(i in alpha_change){
  #Seperating the dataframe by lambda1 and lambda2:
  curr_lambda1 <- which(fitness$changing_value == i & fitness$lambda == lambda1)
  curr_lambda2 <- which(fitness$changing_value == i & fitness$lambda == lambda2)
  
  changed_lambda1 <- as.data.frame(zeroHunter(fitness[curr_lambda1,]))
  changed_lambda2 <- as.data.frame(zeroHunter(fitness[curr_lambda2,]))
  
  
  fitness$lower_quantile[curr_lambda1] <- changed_lambda1$lower_quantile
  fitness$upper_quantile[curr_lambda1] <- changed_lambda1$upper_quantile
  fitness$bootstrapped_r_star[curr_lambda1] <- changed_lambda1$bootstrapped_r_star
  fitness$r_star_mean[curr_lambda1] <- changed_lambda1$r_star_mean
  
  fitness$lower_quantile[curr_lambda2] <- changed_lambda2$lower_quantile
  fitness$upper_quantile[curr_lambda2] <- changed_lambda2$upper_quantile
  fitness$bootstrapped_r_star[curr_lambda2] <- changed_lambda2$bootstrapped_r_star
  fitness$r_star_mean[curr_lambda2] <- changed_lambda2$r_star_mean
}
write.table(fitness, file=plot_file, sep=",", quote=FALSE, append=FALSE)



#This loop is just for making the plot; it extracts each different changing_value 
  #for both lambdas and uses this for the plot. This is just a work-around for
  #the 10,000 long bootstrapped fitness vector stored in the main dataframe.
plotting_dataframe_lambda1 <- c()
plotting_dataframe_lambda2 <- c()
for(i in alpha_change){
  find_lambda1 <- which(fitness$lambda == lambda1 & fitness$changing_value == i)
  find_lambda2 <- which(fitness$lambda == lambda2 & fitness$changing_value == i)
  
  plotting_dataframe_lambda1 <- rbind(plotting_dataframe_lambda1, fitness[find_lambda1[1],])
  plotting_dataframe_lambda2 <- rbind(plotting_dataframe_lambda2, fitness[find_lambda2[1],])
  next
}


ggplot(data=plotting_dataframe_lambda1, aes(x=changing_value)) +
  geom_line(aes(y=r_star_mean)) +
  geom_errorbar(aes(ymin=lower_quantile, ymax=upper_quantile))
