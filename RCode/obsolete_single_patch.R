
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
