#' @include bayesModelFunctions.R
NULL

#' Fit the covariate model for the Bayesian auto-g network computation.
#'
#' \code{agcParam} takes a data frame containing a treatment,
#' outcome, and various covariates to be modeled via an MCMC.
#' The network of observations (rows in the data frame) should
#' be contained in the adjacency matrix (adjmat).
#'
#' B and R specify the number of iterations to be run for the
#' various Bayesian sampling procedures. Defaults are conservative.
#'
#' The seed variable requires a vector of integers to
#' run potentially multiple chains.
#'
#' @param data A data.frame containing variables to be considered
#' for the covariate model.
#'
#' @param treatment A string specifying the column in data
#' that is the treatment effect in the underlying model.
#' Must be a binary (0/1) value.
#'
#' @param outcome A string specifying the column in data
#' that is the outcome of interest in the underlying model.
#' Must be a binary (0/1) value.
#'
#' @param adjmat An adjacency matrix (0/1s) specifying
#' the network structure. The number of rows and columns should
#' match the number of rows in the data object.
#'
#' @param B The number of iterations for the MCMC outer loop.
#' Default = 25,000
#'
#' @param R The number of iterations for the Gibbs inner loop.
#' Default = 10.
#'
#' @param seed An integer vector of the values for \code{set.seed}.
#' By default, only c(1), which runs one chain with a seed of 1. If
#' the vector is of length greater than 1, then multiple chains with
#' both values as a seed will be run in parallel.
#'
#' @return An S3 object of type agcParamClass that contains essential
#' values for the covariate model.
#'
#' @importFrom stats plogis quantile rbinom rmultinom runif var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import methods
#' @import Rcpp
#' @importFrom parallel mclapply
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
#'
#' @export
setGeneric(name = "agcParam",
           def = function(data, treatment, outcome, adjmat,
                          B = 25000, R = 10, seed = c(1))

             standardGeneric("agcParam"))

#' @rdname agcParam
setMethod("agcParam", signature("data.frame", "character", "character", "matrix",
                                "numeric", "numeric", "numeric"),
          definition = function(data, treatment, outcome, adjmat,
                                B, R, seed){

            "%ni%" <- Negate("%in%")

            # Sanity checks for user-specified parameters
            possibilities <- colnames(data)
            stopifnot(outcome %in% possibilities)
            stopifnot(treatment %in% possibilities)
            stopifnot(length(outcome) == 1)
            stopifnot(length(treatment) == 1)

            # Setup covariate dataframe
            covariate_data_frame <- data[ , !(names(data) %in% c(treatment, outcome))]
            stopifnot(dim(covariate_data_frame)[2] > 1)

            # Verify binary outcome and treatment
            outcome_vec <- data[,outcome]
            treatment_vec <- data[,treatment]
            stopifnot(all(outcome_vec %in% c(0,1)))
            stopifnot(all(treatment_vec %in% c(0,1)))

            # Verify dimensions of adjacency matrix
            nadj1 <- dim(adjmat)[1]
            nadj2 <- dim(adjmat)[2]
            ndata <- dim(data)[1]
            if(!all.equal(nadj1, nadj2, ndata)){
              stop("Error: ensure the same number of rows/columns matches the number of rows in the supplied data frame")
            }

            # Preprocess covariates and process output of helper function
            process_covariate <-.covariate_process(covariate_data_frame) # categorical --> binary

            covariate <- process_covariate[[1]]
            zero_pairs <- process_covariate[[2]]
            group_lengths <- process_covariate[[3]]
            group_functions <- process_covariate[[4]]

            data <- cbind(outcome_vec, treatment_vec, covariate) # final dataframe

            # create a new dataframe of indep persons (i.e. no neighbors)
            weights <- apply(adjmat,1,sum) # number of neighbors for everyone
            data.indep <- data[weights==0,] # independent data
            n.indep <- nrow(data.indep)

            # create a new dataframe for non-indep persons (i.e. at least 1 neighbor)
            data <- data[weights>0,] # non-indpendent data
            adjmat <- adjmat[weights>0,weights>0] # non-independent adjacency
            weights <- apply(adjmat,1,sum) # number of neighbors
            N <- nrow(adjmat) # number of non-independent persons
            adjacency <- list(NA) #adjacency list (lists ID for each persons neighbors)

            for (i in 1:N){
              adjacency[[i]] <- which(adjmat[i,]==1)
            }

            #Parameters
            ncov <- ncol(data)-2 #number of covariates / confounders
            nrho <- choose(ncov,2) #number of rho terms in covariate model
            L <- ncov*2 + nrho #number of params in covariate model
            P <- 2 + 2 + ncov*2 #number of params in outcome model

            #individual terms
            ##interconnected units
            outcome.i <- as.matrix(data[1])
            trt.i <- as.matrix(data[2])
            cov.i <- as.matrix(data[3:ncol(data)])

            ##independent units
            outcome.ind <- as.matrix(data.indep[1])
            trt.ind <- as.matrix(data.indep[2])
            cov.ind <- as.matrix(data.indep[3:ncol(data.indep)])

            #sum of neighbors terms
            outcome.n <- (adjmat%*%outcome.i)/weights
            trt.n <- (adjmat%*%trt.i)/weights
            cov.n <- apply(cov.i,2,function(x) {(adjmat%*%x)/weights})

            #sufficient statistics for covariate model

            ##establish which rhos we want to keep or not (i.e. categorical b)
            grid <- t(combn(1:ncov, 2))
            use_rho <- apply(grid, 1, function(pair){
              char <- paste(pair, collapse = "_")
              as.numeric(char %ni% zero_pairs)
            })
            use_rho_all <- c(rep(1,ncov),use_rho,rep(1,ncov))

            ##interconnected units
            sum.l <- c(unname(colSums(cov.i)),
                       unname(colSums(cov.i[,grid[,1], drop = FALSE] * cov.i[,grid[,2], drop = FALSE]))*use_rho,
                       unname(colSums(cov.i*cov.n)))

            ##independent units
            sum.l.indep <- c(unname(colSums(cov.ind)),
                             unname(colSums(cov.ind[,grid[,1], drop = FALSE] * cov.ind[,grid[,2], drop = FALSE]))*use_rho)


            #sufficient statistics for outcome model
            ##interconnected units
            sum.y <- c(sum(outcome.i), #sum individuals outcome
                       sum(outcome.i*trt.i), #sum individuals outcome * treatment
                       (t(outcome.i)%*%cov.i)[1,], #sum individuals outcome * all covs
                       (t(outcome.i)%*%outcome.n)[1,1], #sum neighbors outcome
                       (t(outcome.i)%*%trt.n)[1,1], #sum neighbors treatment
                       (t(outcome.i)%*%cov.n)[1,]#sum covariates treatment
            )

            ##independent units
            sum.y.indep <- c(sum(outcome.ind),
                             sum(outcome.ind*trt.ind),
                             apply(cov.ind,2,function(x) {sum(outcome.ind*x)}))

            # See if there are independent units?
            indepedents_present <- dim(data.indep)[1] != 0

            # Includes the independent individual terms from the outcome model
            if(indepedents_present) design.mat.indep <- as.matrix(cbind(1,data.indep[-1]))

            # Establish parameters for the MCMC
            alpha <- matrix(NA,B+1,L)
            beta <- matrix(NA,B+1,P)
            accept_alpha <- NA
            accept_beta <- NA

            #for independent terms denominator
            l_grid <- unname(data.matrix(expand.grid(replicate(ncov, c(0,1), simplify=FALSE))))

            # Have to change the index of adjacency
            for (i in 1:N){
              adjacency[[i]] <- adjacency[[i]] - 1
            }

            # Last parameters needed for MCMC
            J <- dim(covariate)[2] # number of covariates

            # Parallelize over the chains
            mclapply(seed, function(seed_value){

              set.seed(seed_value)

              # Main loop for the MCMC
              pb <- txtProgressBar(min = 1, max = B, style = 3)

              for (b in 1:B){
                #Step 0. Starting values
                if (b==1) {
                  alpha[b,] <- MASS::mvrnorm(1,rep(0,L),.1*diag(L))
                  beta[b,] <- MASS::mvrnorm(1,rep(0,P),.1*diag(P))
                }

                ##ALPHA##
                #Step 1. Proposal
                scale <- .005
                alpha.p <- MASS::mvrnorm(1,alpha[b,],scale*diag(L))*use_rho_all

                #assign tau, rho, nu
                tau.p <- alpha.p[1:ncov]
                rho.p <- alpha.p[(ncov+1):(ncov+nrho)]
                nu.p <- alpha.p[(ncov+nrho+1):L]

                # Make symmetrical matrix of the rho values
                rho_mat <- matrix(0, nrow = J, ncol = J)
                rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho.p; rho_mat <- rho_mat + t(rho_mat)
                aux.p.mat <- auxVarCpp(tau.p,rho.p,nu.p,N,R,J, rho_mat,
                                       adjacency,cov.i,weights,group_lengths,group_functions)
                aux.p.cov.n <- apply(aux.p.mat,2,function(x) {(adjmat%*%x)/weights})
                sum.aux.p <- c(unname(colSums(aux.p.mat)),
                               unname(colSums(aux.p.mat[,grid[,1, drop = FALSE], drop = FALSE] * aux.p.mat[,grid[,2, drop = FALSE], drop = FALSE]))*use_rho,
                               unname(colSums(aux.p.mat*aux.p.cov.n)))

                #Step 3. Calculate accept-reject ratio (everything is log-ed!)
                #interconnected units
                h.l1.p <- alpha.p%*%sum.l
                h.l1.c <- alpha[b,]%*%sum.l
                h.aux.p <- alpha.p%*%sum.aux.p
                h.aux.c <- alpha[b,]%*%sum.aux.p

                #independent units
                if(indepedents_present){
                  f.p.num <- alpha.p[1:(ncov+nrho)]%*%sum.l.indep
                  f.p.denom <- log(sum(v_get_sum(1:dim(l_grid)[1], l_grid, tau.p, rho.p, ncov)))
                  f.c.num <- alpha[b,1:(ncov+nrho)]%*%sum.l.indep
                  f.c.denom <- log(sum(v_get_sum(1:dim(l_grid)[1], l_grid, alpha[b,1:ncov], alpha[b,(ncov+1):(ncov+nrho)], ncov)))

                  #priors
                  #prior.p <- mvtnorm::dmvt(alpha.p,delta=rep(0,L),sigma=4*diag(L),df=3,log=TRUE)
                  #prior.c <- mvtnorm::dmvt(alpha[b,,c],delta=rep(0,L),sigma=4*diag(L),df=3,log=TRUE)

                  ratio <- h.l1.p + h.aux.c - h.l1.c - h.aux.p + f.p.num - n.indep*f.p.denom - f.c.num + n.indep*f.c.denom
                } else {
                  ratio <- h.l1.p + h.aux.c - h.l1.c - h.aux.p
                }

                accept <- rbinom(1,1,min(1,exp(ratio)))
                if (accept == 1){alpha[b+1,] <- alpha.p} else { alpha[b+1,] <- alpha[b,] }
                accept_alpha[b] <- accept

                ##BETA##
                #Step 1. Proposal
                scale <- .005
                beta.p <- MASS::mvrnorm(1,beta[b,],scale*diag(P))

                #Step 2. Auxilary variable based on proposal
                aux.p <- aux.var.outcome.cl(beta.p,trt.i,cov.i,N,10,adjacency,outcome.i,weights)
                sum.aux.p <- c(sum(aux.p),
                               sum(aux.p*trt.i),
                               apply(cov.i,2,function(x) {sum(aux.p*x)}),
                               sum(aux.p*(adjmat%*%aux.p)/weights),
                               t(aux.p)%*%trt.n,
                               apply(cov.n,2,function(x) {t(aux.p)%*%x}))

                #Step 3. Calculate accept-reject ratio (everything is log-ed!)
                #interconnected units
                h.p <- beta.p%*%sum.y
                h.c <- beta[b,]%*%sum.y
                h.aux.p <- beta.p%*%sum.aux.p
                h.aux.c <- beta[b,]%*%sum.aux.p

                if(indepedents_present){

                  #independent units

                  f.p.num <- beta.p[1:(2+ncov)]%*%sum.y.indep
                  f.p.denom <- sum(log(1+ exp(design.mat.indep%*%beta.p[1:(2+ncov)])))
                  f.c.num <- beta[b,1:(2+ncov)]%*%sum.y.indep
                  f.c.denom <- sum(log(1+exp(design.mat.indep%*%beta[b,1:(2+ncov)])))

                  #priors
                  #prior.p <- mvtnorm::dmvt(beta.p,delta=rep(0,P),sigma=4*diag(P),df=3,log=TRUE)
                  #prior.c <- mvtnorm::dmvt(beta[b,,c],delta=rep(0,P),sigma=4*diag(P),df=3,log=TRUE)

                  ratio <- h.p + h.aux.c - h.c - h.aux.p + f.p.num - f.p.denom - f.c.num + f.c.denom
                } else {
                  ratio <- h.p + h.aux.c - h.c
                }
                accept <- rbinom(1,1,min(1,exp(ratio)))
                if (accept == 1){beta[b+1,] <- beta.p} else { beta[b+1,] <- beta[b,] }
                accept_beta[b] <- accept
                setTxtProgressBar(pb, b)
              }

              #save labels
              rho.terms <- apply(grid, 1, function(pair){char <- paste(pair, collapse = "_")})
              cov.lab <- c(colnames(covariate),rho.terms,paste("neighbors",colnames(covariate)))
              outcome.lab <- c(colnames(data[,1:2]),colnames(covariate),
                               paste("neighbors",colnames(data[,1:2])),
                               paste("neighbors",colnames(covariate)))

              colnames(alpha) <- cov.lab
              colnames(beta) <- outcome.lab

              outlist <- list(alpha,beta,accept_alpha,accept_beta,use_rho_all,group_lengths,group_functions)
              names(outlist) <- c("alpha", "beta", "accept_alpha", "accept_beta", "use_rho_all", "group_lengths", "group_functions")

              return(outlist)

            }) -> outoutlist

            # Name according to the chain
            names(outoutlist) <- paste0("chain", as.character(1:length(outoutlist)))

            attr(outoutlist, "class") <- "agcParamClass"
            return(outoutlist)
          })
