#' @include autognet.R
NULL

# Internal function that potentially expands
# a data frame given a multinomial covariate
# and further classifies the appropriate covariate
# function for the specific covariate model

# Return a two item list:
# 1) expanded data frame; 2) zero_list
# 3) group lengths vector; 4) group functions vector

#' @importFrom stats model.matrix
#' @importFrom utils combn
#' @importFrom MASS mvrnorm
.covariate_process <- function(test_df){

  # Initialize values
  final_index <- 0; zero_list_internal <- c()
  group_lengths <- c(); group_functions <- c()

  #----------------------------
  # group_functions dictionary
  # 1 - binary
  # 2 - multinomial
  #----------------------------

  # Loop over supplied covariates and establish output vectors
  lapply(1:dim(test_df)[2], function(cov_index){

    # Pull column names
    column_name <- colnames(test_df)[cov_index]

    # If already binary, don't do anything
    if(length(unique(test_df[,cov_index])) == 2 & is.numeric(test_df[,cov_index])){
      final_index <<- final_index + 1
      group_lengths <<- c(group_lengths, 1)
      group_functions <<- c(group_functions, 1)

      test_df[,cov_index, drop = FALSE]

    } else {

      # Multiple levels -- figure out how to binarize
      levels <- levels(test_df[,cov_index])
      total_out_columns <- length(levels) -1

      # Blacklist rho terms if we are expanding here
      if(total_out_columns > 1){
        new_zero_list_items <- apply(combn(final_index + 1:total_out_columns, 2), 2, function(pair){
          paste(pair, collapse = "_")
        })
        zero_list_internal <<- c(zero_list_internal, new_zero_list_items)
      }
      final_index <<- final_index + total_out_columns

      group_lengths <<- c(group_lengths, total_out_columns)

      # If user supplied a vector with two levels, we have to do these
      # steps to recode it 0/1 but we want to treat it as a logistic
      # in the downstream; thus, this group function should be 1
      classify_out <- ifelse(is.factor(test_df[,cov_index]) & length(unique(test_df[,cov_index])) == 2, 1, 2)
      group_functions <<- c(group_functions, classify_out)

      # Establish covariate matrix
      vecToBeRenamed <-  test_df[,cov_index]
      binary_matrix <- model.matrix(~vecToBeRenamed)[,-1, drop = FALSE]
      colnames(binary_matrix) <- gsub("vecToBeRenamed", paste0(column_name, "_"), colnames(binary_matrix))

      binary_matrix
    }
  }) -> list_cov_expand

  # Prepare an S3 object to be returned
  listOut <- list(
    do.call(cbind,list_cov_expand),
    zero_list_internal,
    group_lengths,
    group_functions
  )
  names(listOut) <- c("df", "zero_list", "group_lengths", "group_functions")

  return(listOut)
}


# create coverage intervals
cov_int <- function(x) {
  return(quantile(x,probs=c(.025,.975)))
}

# to check convergence
psr <- function(samples) {
  ##
  R <- nrow(samples)
  ##
  B <- R * var(apply(samples, 2, mean))
  W <- mean(apply(samples, 2, var))
  ##
  value <- sqrt((B + (R-1)*W) / (R*W))
  return(value)
}

# Helper function for getting a sum weighted
# by neighbor values -- specific to these internal data structures
get_sum <- function(idx, l_grid, tau, rho, n){

  # Establish Ls for this iteration
  l_vec <- l_grid[idx,]

  # Izzie / Caleb formulation of sum--
  # Make a temporary matrix; fill by row so have to transpose later
  # see https://stackoverflow.com/questions/30787317/a-vector-to-an-upper-triangle-matrix-by-row-in-r
  temp_m <- matrix(0, n, n); temp_m[lower.tri(temp_m, diag=FALSE)] <- rho

  # intercept             main
  exp((tau%*%l_vec)[1,1] + sum(l_vec %*% t(l_vec) * t(temp_m)))
}

# vectorize the above over the index of the individual in the network
v_get_sum <- Vectorize(get_sum, "idx")

