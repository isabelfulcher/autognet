#' @include gibbs-factors_DMH.R
NULL

#' Evaluate the network causal effects from estimated Gibbs factors.
#'
#' \code{agcEffect} takes a covariate and outcome model fit
#' from \code{agcParam}, parameters for an MCMC, and the
#' specified treatment effect and determines the causal estimates (direct and spillover)
#' generated from the network structure.
#'
#' ADD DISCUSSION OF DIFFERENT TREATMENT REGIME OPTIONS by linking to URL
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
#' @param subset The indices of the individuals, as they appear in the
#' adjacency matrix, to be included in the network causal effects estimates.
#' Default = 0 (include everyone).
#'
#' @param a_fixed Izzie to Update
#'#'
#' @param dynamic_coef_vec Izzie to Update
#'
#' @param dynamic_among_treated Izzie to Update
#'
#' @param dynamic_single_edge Izzie to Update
#'
#' @param R The number of iterations for the Gibbs inner loop.
#' Default = 10.
#'
#' @param burnin_R The index to start evaluation as one would normally
#' have for a burnin for a Bayesian computation. Default = 10.
#'
#' @param burnin_cov The index to start saving covariates for use in the
#' network causal effect chains. Default = 10.
#'
#' @param average An indicator of whether to evaluate the causal effects as an average
#' of the R iterations. Default = "TRUE".
#'
#' @param index_override The MCMC outer loop iteration numbers to include
#' in the evaluation of the causal effects. This will override the burnin and
#' thin parameters. Default = 0 (all iterations included).
#'
#' @param return_effects An indicator of whether to return direct and spillover effects
#' or just the counterfactual expectations. Use the latter if you wish to construct overall
#' effects or different types of spillover effects. Default = 1.
#'
#' @param p_vec Izzie to do
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
#' effects <- agcEffect(mod)
#'
#'
#' @export
setGeneric(name = "agcEffect",
           def = function(mod, burnin = 1, thin = 0.5, treatment_allocation = 0.5, subset = 0,
                          a_fixed = c(0), dynamic_coef_vec = c(0), dynamic_among_treated = 0, dynamic_single_edge = 0,
                          R = 10, burnin_R = 10, burnin_cov = 10, average = TRUE, index_override = 0,
                          return_effects = 1, p_vec = c(0.5,0.5))
             standardGeneric("agcEffect"))

#' @rdname agcEffect
setMethod("agcEffect", signature("list", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY"),
          definition = function(mod, burnin = 1, thin = 0.2, treatment_allocation = 0.5, subset = 0,
                                a_fixed = c(0), dynamic_coef_vec = c(0), dynamic_among_treated = 0, dynamic_single_edge = 0,
                                R = 10, burnin_R = 10, burnin_cov = 10, average = TRUE, index_override = 0,
                                return_effects = 1, p_vec = c(0.5,0.5)){

            stopifnot(treatment_allocation >= 0 & treatment_allocation <= 1)
            stopifnot(return_effects == 0 | return_effects == 1)
            stopifnot(thin> 0 & thin <= 1)
            stopifnot("agcParamClass" %in% class(mod))
            stopifnot(R >= 1)

            # Determine which indices to actually compute
            total <- dim(mod[[1]])[1]
            indices <- seq(from = 1 + burnin, to = total, by = round(1/thin))
            noprog <- as.numeric(length(indices)==1)

            # make up value so that the variable is defined
            # avoids an error but does nothing
            treated_indicator = c(0,0)

            # Manually establish the indices that will be evaluted in the loop
            if(index_override != 0){
              indices <- index_override
            }

            if(length(indices) < 1){
              stop("Change burnin and thin values -- no valid indices will currently be computed")
            }

            # Define parameters from the covariate model
            alpha <- mod[[1]]
            beta <- mod[[2]]
            group_lengths <- mod[[6]]
            group_functions <- mod[[7]]
            adjmat <- mod[[8]]
            weights <- apply(adjmat,1,sum)

            # Define some other parameters using IF's logic
            ncov <- (ncol(beta)-4)/2
            nrho <- choose(ncov,2)
            L <- ncov*2 + nrho
            N <- nrow(adjmat)

            # Modify fake data to test nu off diagonal terms
            if(FALSE){
              alpha <- cbind(alpha, rnorm(11,0, 0.1), rnorm(11,0, 0.1))
              colnames(alpha)[4:7] <- c("new_11", "new_22","new_12", "new_21")
              alpha <- alpha[,c(1,2,3,4,6,5,7)]
            }

            if(ncol(alpha) > L){
                L <- ncov + nrho + ncov^2
            }

            # Establish a C++-friendly adjacency list
            adjacency <- list(NA)
            for (i in 1:N){
              adjacency[[i]] <- unname(which(adjmat[i,]==1)) - 1 # the -1 is for zero indexing in Cpp
            }

            # Establish C++-friendly dynamic single edge case for edge case of one covariate
            dynamic_single_edge <- dynamic_single_edge - 1

            # Figure out the subset logic
            # WE ALREADY DO THE ZERO-based indexing in C++ DON'T DO IT AGAIN HERE
            if(subset == 0){
              subset <- 1:N
            } else {
              subset <- subset
            }

            # Now evaluate each chain that results from the burnin and thin
            if (noprog==0){pb <- txtProgressBar(min = 1, max = length(indices), style = 3)}
            if (return_effects == 1){
              mat <- matrix(NA, nrow = length(indices), ncol = 3)
            } else {
              mat <- matrix(NA, nrow = length(indices), ncol = 4)
            }
            for(idxx in 1:length(indices)){
              b <- indices[idxx]
              tau <- alpha[b,1:ncov]
              rho <- alpha[b,(ncov+1):(ncov+nrho)]
              nu  <- alpha[b,(ncov+nrho+1):L]

              # Reformat to matrix for off-diagonal nu terms
              if(length(nu) == ncov){
                nu_mat <- diag(nu)
                additional_nu <- 0
              } else {
                nu_mat <- matrix(nu, nrow = ncov, byrow = TRUE)
                additional_nu <- 1
              }
              nu <- nu_mat

              rho_mat <- matrix(0, nrow = ncov, ncol = ncov)
              rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho; rho_mat <- rho_mat + t(rho_mat)


              # Make starting values
              cov_mat <- matrix(rbinom(ncov*N, 1, 0.5), ncol = ncov, nrow = N )


              ## MAKE COVARIATE ARRAY ##
              cov.list <- networkGibbsOutCovCpp(tau, rho, nu,
                                                ncov, R + burnin_R, N, burnin_cov, rho_mat,
                                                adjacency, weights, cov_mat,
                                                group_lengths, group_functions, additional_nu)

              ## For dynamic among treated setup ##
              if(dynamic_among_treated==1){
                treated_indicator = 1:N %in% subset
              }
              ## NON-INDIVIDUAL GIBBS (overall effects) ##
              psi_gamma <-mean(networkGibbsOuts1Cpp(cov_list = cov.list, beta = beta[b,], p = treatment_allocation,
                                                    a_fixed = a_fixed,  dynamic_coef_vec = dynamic_coef_vec,
                                                    dynamic_among_treated = dynamic_among_treated, dynamic_single_edge = dynamic_single_edge,
                                                    ncov = ncov, R= R + burnin_R, N = N,  # have to do R + burnin for weird C++ error
                                                    adjacency = adjacency, weights = weights, treated_indicator = treated_indicator,
                                                    burnin = burnin_R, average = as.numeric(average), p_vec = p_vec)[subset])

              psi_zero <-mean(networkGibbsOuts1Cpp(cov_list = cov.list, beta = beta[b,], p = 0,
                                                   a_fixed = a_fixed,  dynamic_coef_vec = dynamic_coef_vec,
                                                   dynamic_among_treated = dynamic_among_treated, dynamic_single_edge = dynamic_single_edge,
                                                   ncov = ncov, R = R + burnin_R, N = N, # have to do R + burnin for weird C++ error
                                                   adjacency = adjacency, weights = weights, treated_indicator = treated_indicator,
                                                   burnin = burnin_R, average = as.numeric(average), p_vec = p_vec)[subset])

              psi_1_gamma <- mean(networkGibbsOuts2Cpp(cov_list = cov.list, beta = beta[b,], p = treatment_allocation,
                                                       a_fixed = a_fixed,  dynamic_coef_vec = dynamic_coef_vec,
                                                       dynamic_among_treated = dynamic_among_treated, dynamic_single_edge = dynamic_single_edge,
                                                       ncov = ncov, R = R + burnin_R, N = N,   # have to do R + burnin for weird C++ error
                                                       adjacency = adjacency, weights = weights, treated_indicator = treated_indicator,
                                                       subset = subset,
                                                       treatment_value = 1.0,
                                                       burnin = burnin_R, average = as.numeric(average), p_vec = p_vec))

              psi_0_gamma <- mean(networkGibbsOuts2Cpp(cov_list = cov.list, beta = beta[b,], p = treatment_allocation,
                                                       a_fixed = a_fixed,  dynamic_coef_vec = dynamic_coef_vec,
                                                       dynamic_among_treated = dynamic_among_treated, dynamic_single_edge = dynamic_single_edge,
                                                       ncov = ncov, R = R + burnin_R, N = N,   # have to do R + burnin for weird C++ error
                                                       adjacency = adjacency, weights = weights, treated_indicator = treated_indicator,
                                                       subset = subset,
                                                       treatment_value = 0,
                                                       burnin = burnin_R, average = as.numeric(average), p_vec = p_vec))

              if (return_effects == 1){
                return <- c(psi_gamma, psi_1_gamma - psi_0_gamma, psi_0_gamma - psi_zero)
              } else {
                return <- c(psi_gamma, psi_zero, psi_1_gamma, psi_0_gamma)
              }

              mat[idxx,] <- return
              if(noprog==0){setTxtProgressBar(pb, idxx)}
            }

            output <- mat

            if (return_effects == 1){
              colnames(output) <- c("average", "direct", "spillover")
            } else {
              colnames(output) <- c("average", "psi_zero", "psi_1_gamma","psi_0_gamma")
            }

            return(output)
          })
