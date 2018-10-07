#' @include helperFunctions.R
NULL

################################################
##### SET 1: AUXILARY VARIABLE CALCULATION #####
################################################

## 1.1 COVARIATE MODELS ##
aux.var.cl <- function(tau, rho, nu, N, R, J, rho_mat,
                       adjacency, cov.i, weights, group_lengths, group_functions){

  cov.mat <- cov.i

  # Number of iterations
  for (r in 1:R){

    # Number of people
    for (i in 1:N){

      # Index covariate
      j <- 1

      # Number of groups (of covariants)
      for (group_index in 1:length(group_lengths)){

        # Values associated with each group
        group_length <- group_lengths[group_index]
        group_function <- group_functions[group_index]

        # group_length is the number of binarized covariates
        # group_function is a numeric dictionary key that utilizes a value from the covariate_process function
        # if group_length is > 1, meaning that we have multiple binarized covariates associated with a specific
        # entity, then it is definitely a multinomial
        if(group_length > 1){

          # Multinomial case
          prob_vec <- sapply(0:(group_length -1), function(m){
            j_prime <- j + m
            exp(tau[j_prime] + sum(rho_mat[,j_prime]*cov.mat[i,]) + nu[j_prime]*sum(cov.mat[adjacency[[i]],j_prime]/weights[i]))
          })
          #
          pro_vec <- c(prob_vec,1)/sum(c(prob_vec,1))
          #message(pro_vec)
          # for the rmultinom call, have to append a 1 and remove the last value; update several values
          cov.mat[i, j + (0:(group_length -1))] <- (rmultinom(1,1,pro_vec)[,1])[-1*group_length]

        } else if(group_function == 1){

          # Logistic / binary case
          prob_Lj <- plogis(tau[j] + sum(rho_mat[,j]*cov.mat[i,]) + nu[j]*sum(cov.mat[adjacency[[i]],j]/weights[i]))
          cov.mat[i,j] <- rbinom(1,1,prob_Lj)
        } # add in normal here once form is decide

        j <- j + group_length

      }
    }
  }
  return(cov.mat)
}

# 1.2 OUTCOME MODEL #

aux.var.outcome.cl <- function(beta,trt,cov,N,R,adjacency,start,weights){
  # beta is vector of parameters for outcome model (beta.p in 04a)
  # trt is Nx1 matrix of treatment values (trt.i in 04a)
  # cov is a Nxp matrix of covariate values (cov.i in 04a)
  # N is number of individuals in network THAT HAVE AT LEAST ONE NEIGHBOR
  # R is number of iterations to run through
  # adjacency is a list of length N with neighbor IDs for each individual
  # start is starting point for chain (outcome.i in 04a)
  # weights is a vector of length N with the number of neighbors for each individual

  vec <- start
  for (r in 1:R){
    for (i in 1:N){
      prob_outcome <- plogis(beta[1] + # intercept term
                               beta[2]*trt[i] + # individual treatment term
                               beta[3:(2+ncol(cov))]%*%cov[i,] + # individual covariate term(s)
                               beta[(3+ncol(cov))]*sum(vec[adjacency[[i]]]/weights[i]) + # neighbor outcome
                               beta[(4+ncol(cov))]*sum(trt[adjacency[[i]]]/weights[i]) + # neighbor treatment
                               sum(beta[(4+ncol(cov)+(1:ncol(cov)))]*colSums(cov[adjacency[[i]],(1:ncol(cov)),drop = FALSE]/weights[i])) # neighbor covariate
      )


      vec[i] <- rbinom(1,1,prob_outcome)
    }
  }
  return(vec)
}

#######################################################
##### SET 3: GIBBS SAMPLER FOR CAUSAL ESTIMANDS #######
#######################################################


## NEW FUNCTION FOR COVARIATES ##

#This will generate an array of covariates to use in the next functions
network.gibbs.cov <-  function(tau, rho, nu, ncov, R, N, adjacency, weights,
                               group_lengths, group_functions){

  # tau is a vector of length ncov
  # rho is a vector of lenth choose(ncov,2)
  # nu is a vector of length ncov
  # beta is a vector of length P
  # ncov is the number of covariates
  # R is the number of iterations
  # N is the total network size (includes indep and non indep)
  # p is the treatment allocation
  # adjacency is a list of neighbors for each pesron (only non indep nodes)
  # weights are the weights for each person (includes indep and non indep)

  J <- ncov #easier to use this notation

  #storage
  cov.save <- list()
  cov_mat <- matrix(rbinom(N*J, 1, runif(J,.1,.9)),N,J)

  #Make symmetric matrix of the rho values (needed for covariate iterations)
  rho_mat <- matrix(0, nrow = J, ncol = J)
  rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho; rho_mat <- rho_mat + t(rho_mat)

  for (r in 1:R){

    #Gibbs (conditional probabilities for each person)
    for (i in 1:N){

      if(weights[i] > 0){
        #Generate L
        j <- 1
        # Number of groups (of covariants)
        for (group_index in 1:length(group_lengths)){

          # Values associated with each group
          group_length <- group_lengths[group_index]
          group_function <- group_functions[group_index]

          # group_length is the number of binarized covariates
          # group_function is a numeric dictionary key that utilizes a value from the covariate_process function
          # if group_length is > 1, meaning that we have multiple binarized covariates associated with a specific
          # entity, then it is definitely a multinomial
          if(group_length > 1){

            # Multinomial case
            prob_vec <- sapply(0:(group_length -1), function(m){
              j_prime <- j + m
              exp(tau[j_prime] + sum(rho_mat[,j_prime]*cov_mat[i,]) + nu[j_prime]*sum(cov_mat[adjacency[[i]],j_prime]/weights[i]))
            })

            # for the rmultinom call, have to append a 1 and remove the last value; update several values
            cov_mat[i, j + (0:(group_length -1))] <- (rmultinom(1,1,c(prob_vec,1))[,1])[-1*group_length]

          } else if(group_function == 1){

            # Logistic / binary case
            prob_Lj <- plogis(tau[j] + sum(rho_mat[,j]*cov_mat[i,]) + nu[j]*sum(cov_mat[adjacency[[i]],j]/weights[i]))
            cov_mat[i,j] <- rbinom(1,1,prob_Lj)
          } # add in normal here once form is decide

          j <- j + group_length

        }
      } else {
        #Generate L
        j <- 1
        # Number of groups (of covariants)
        for (group_index in 1:length(group_lengths)){

          # Values associated with each group
          group_length <- group_lengths[group_index]
          group_function <- group_functions[group_index]

          # group_length is the number of binarized covariates
          # group_function is a numeric dictionary key that utilizes a value from the covariate_process function
          # if group_length is > 1, meaning that we have multiple binarized covariates associated with a specific
          # entity, then it is definitely a multinomial
          if(group_length > 1){

            # Multinomial case
            prob_vec <- sapply(0:(group_length -1), function(m){
              j_prime <- j + m
              exp(tau[j_prime] + sum(rho_mat[,j_prime]*cov_mat[i,]))
            })

            # for the rmultinom call, have to append a 1 and remove the last value; update several values
            cov_mat[i, j + (0:(group_length -1))] <- (rmultinom(1,1,c(prob_vec,1))[,1])[-1*group_length]

          } else if(group_function == 1){

            # Logistic / binary case
            prob_Lj <- plogis(tau[j] + sum(rho_mat[,j]*cov_mat[i,]))
            cov_mat[i,j] <- rbinom(1,1,prob_Lj)
          } # add in normal here once form is decide

          j <- j + group_length

        }
      }
    }
    cov.save[[r]] <- cov_mat
  }
  return(cov.save)
}

#This is just the loop over individuals for some precified gamma
network.gibbs.out1 <- function(cov.list,beta,p,
                               R,N,adjacency,weights){
  # beta is a vector of length P

  J <- ncov #easier to use this notation


  #storage
  prob_outcome_1 <- matrix(NA,R,N) #this is the E(Y_i(a_bar)), or psi

  #starting values
  outcome_1 <- rbinom(N,1,runif(1,.1,.9))
  a_bar <- rbinom(N,1,p)

  for (r in 1:R){
    #Gibbs (conditional probabilities for each person)
    cov_mat <- cov.list[[r]]

    for (i in 1:N){
    if (p != 0){
      a_bar[i] <- rbinom(1,1,p)

      #Generate Y given A,L
      prob_outcome_1[r,i] <- plogis(beta[1] + # intercept term
                                  beta[2]*a_bar[i] + # individual treatment term
                                  beta[3:(2+J)]%*%cov_mat[i,] + # individual covariate term(s)
                                  beta[(3+J)]*sum(outcome_1[adjacency[[i]]]/weights[i]) + # neighbor outcome
                                  beta[(4+J)]*sum(a_bar[adjacency[[i]]]/weights[i]) + # neighbor treatment
                                  sum(beta[(4+J+(1:J))]*colSums(cov_mat[adjacency[[i]],(1:J), drop = FALSE]/weights[i]))) # neighbor covariate

      outcome_1[i] <- rbinom(1,1,prob_outcome_1[r,i])
    } else {
      prob_outcome_1[r,i] <- plogis(beta[1] + # intercept term
                                      beta[3:(2+J)]%*%cov_mat[i,] + # individual covariate term(s)
                                      beta[(3+J)]*sum(outcome_1[adjacency[[i]]]/weights[i]) + # neighbor outcome
                                      sum(beta[(4+J+(1:J))]*colSums(cov_mat[adjacency[[i]],(1:J), drop = FALSE]/weights[i]))) # neighbor covariate

      outcome_1[i] <- rbinom(1,1,prob_outcome_1[r,i])
    }
    }
  }
  output <- mean(prob_outcome_1[R,]) #could just output the end result for each person (then take mean outside)
  return(output)
}


network.gibbs.out2 <- function(cov.list,beta,p,
                               R,N,adjacency,weights,
                               person,treatment_value){
  # person is the index for the person whose covariate value is fixed
  # treatment_value is the value that treatment should take (1 or 0)

  J <- ncov #easier to use this notation


  #storage
  prob_outcome <- matrix(NA,R,N)

  #starting values
  outcome <- rbinom(N,1,runif(1,.1,.9))
  a_bar <- rbinom(N,1,p)

  for (r in 1:R){
    #Gibbs (conditional probabilities for each person)
    cov_mat <- cov.list[[r]]
    for (i in 1:N){
        a_bar[i] <- rbinom(1,1,p)
        a_bar[person] <- treatment_value

        #Generate Y given A,L
        prob_outcome[r,i] <- plogis(beta[1] + # intercept term
                                        beta[2]*a_bar[i] + # individual treatment term
                                        beta[3:(2+J)]%*%cov_mat[i,] + # individual covariate term(s)
                                        beta[(3+J)]*sum(outcome[adjacency[[i]]]/weights[i]) + # neighbor outcome
                                        beta[(4+J)]*sum(a_bar[adjacency[[i]]]/weights[i]) + # neighbor treatment
                                        sum(beta[(4+J+(1:J))]*colSums(cov_mat[adjacency[[i]],(1:J), drop = FALSE]/weights[i]))) # neighbor covariate

        outcome[i] <- rbinom(1,1,prob_outcome[r,i])
    }
  }
  output <- prob_outcome[R,person]
  return(output)
}
