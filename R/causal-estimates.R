#' @include gibbs-factors_DMH.R
NULL

#' Evaluate the out model for the Bayesian auto-g network computation.
#'
#' \code{agcEffect} takes a covariate model fit from \code{agcParam},
#' parameters for an MCMC, and the specified treatment effect and
#' determines the causal estimates (direct and spillover) generated
#' from the network structure.
#'
#'
#' @param mod An \code{agcParamClass} object from \code{agcParam}
#' containing variables to be considered for the covariate model.
#'
#' @param burnin The index to start evaluation as one would normally
#' have for a burnin for a Bayesian computation. Default = 1
#'
#' @param thin The rate at which to evaluate the outcome model.
#' The closer to 1, the more values to compute (supplying 1 will
#' not thin at all).
#'
#' @param treatment_allocation The proportion of the individuals
#' in the network to receive the treatment. Should be a number between
#' 0 and 1. Default = 0.5.
#'
#' @param R The number of iterations for the Gibbs inner loop.
#' Default = 10.
#'
#' @param seed An integer value for \code{set.seed}.
#' This will be applied uniformly to all chains. Default = 1.
#'
#' @return An S3 object of type \code{agcEffectClass} that contains essential
#' values for the outcome model.
#'
#'
#' @examples
#'
#' rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
#' adjmat <- rdsIn[[1]]
#' data <- rdsIn[[2]]
#' treatment <- "treatment"
#' outcome <- "outcome"
#' B = 10
#' R = 5
#' seed = 1
#'
#' mod <- agcParam(data, treatment, outcome, adjmat,
#'                          B = B, R = R, seed = c(1))
#' effects <- agcEffect(mod, burnin = 1, thin = 0.5, treatment_allocation = 0.5, R = 10, seed = 1)
#'
#'
#' @export
setGeneric(name = "agcEffect",
           def = function(mod, burnin = 1, thin = 0.5, treatment_allocation = 0.5, R = 10, seed = 1)
             standardGeneric("agcEffect"))

#' @rdname agcEffect
setMethod("agcEffect", signature("list", "ANY", "ANY", "ANY", "ANY", "ANY"),
          definition = function(mod, burnin = 1, thin = 0.2, treatment_allocation = 0.5, R = 10, seed = 1){

            stopifnot(treatment_allocation > 0 & treatment_allocation < 1)
            stopifnot(thin> 0 & thin <= 1)
            stopifnot("agcParamClass" %in% class(mod))
            stopifnot(R >= 1)

            # Determine which indices to actually compute
            total <- dim(mod[[1]][[1]])[1]
            indices <- seq(from = 1 + burnin, to = total, by = round(1/thin))

            if(length(indices) < 1){
              stop("Change burnin and thin values -- no valid indices will currently be computed")
            }

            # Now evaluate the model evaluating each chain in parallel
            out_names <- names(mod)
            mclapply(out_names, function(chain_name){

              # Define parameters from the covariate model
              alpha <- mod[[chain_name]][[1]]
              beta <- mod[[chain_name]][[2]]
              group_lengths <- mod[[chain_name]][[6]]
              group_functions <- mod[[chain_name]][[7]]
              adjmat <- mod[[chain_name]][[8]]
              weights <- apply(adjmat,1,sum)

              # Define some other parameters using IF's logic
              ncov <- (ncol(beta)-4)/2 #this should just flow from 04a R script
              nrho <- choose(ncov,2) #ibid
              L <- ncov*2 + nrho #ibid

              estimates <- matrix(NA,total,3)

              # Set up the R adj. list
              N <- nrow(adjmat)
              adjacency <- list(NA)
              for (i in 1:N){
                adjacency[[i]] <- which(adjmat[i,]==1)
              }

              # Loop over only value nominated
              for (b in indices){

                tau <- alpha[b+1,1:ncov]
                rho <- alpha[b+1,(ncov+1):(ncov+nrho)]
                nu  <- alpha[b+1,(ncov+nrho+1):L]

                ## CAUSAL EFFECT ##
                estimates[b,] <- network.gibbs3.cl(tau, rho, nu, beta[b+1,],
                                                   ncov, R, N, treatment_allocation, adjacency, weights,
                                                   group_lengths, group_functions)

              }
              estimates
            }) -> outoutlist


            # Prepare the return values
            names(outoutlist) <- out_names

            class(outoutlist) <- append(class(outoutlist),"agcEffectClass")
            return(outoutlist)
          })
