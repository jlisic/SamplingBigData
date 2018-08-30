

split_sample <- function(
  prob,
  delta = exp(-16)
) {

  r_prob <- range(prob)
  n <- length(prob)

  # send our data to the C program
  r.result <- .C("R_split_sample",
                 as.double( prob ),   # probability vector
                 as.integer( n ),                   # length of prob and nrow of x 
                 as.double(delta)
  )
    
  return( which( round(r.result[[1]]) == 1  )) 
}







