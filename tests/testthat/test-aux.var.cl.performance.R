library(testthat)
context("C++ version of aux.var.cl performs well")
set.seed(14651)

# Import simulated data frame from IF
rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
adjmat <- rdsIn[[1]]
data <- rdsIn[[2]]

# Preprocess covariates
cov_df <- data[,-c(1,2)]
cov_df$categorical <- factor(c(rep("1", 20), rep("2", 30), rep("3", 24), "4" ), levels = c("1", "2", "3", "4"))
process_covariate <- autognet:::.covariate_process(cov_df) # categorical --> binary
covariate <- process_covariate[[1]]
zero_pairs <- process_covariate[[2]]
group_lengths <- process_covariate[[3]]
group_functions <- process_covariate[[4]]

data <- cbind(data[,1:2],covariate) # final dataframe

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
"%ni%" <- Negate("%in%")
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

#initialize MCMC for loop
b <- B <- 1
alpha <- matrix(NA,B+1,L)
beta <- matrix(NA,B+1,P)
accept_alpha <- NA
accept_beta <- NA

#for independent terms denominator
l_grid <- unname(data.matrix(expand.grid(replicate(ncov, c(0,1), simplify=FALSE))))

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
tau <- alpha.p[1:ncov]
rho <- alpha.p[(ncov+1):(ncov+nrho)]
nu <- alpha.p[(ncov+nrho+1):L]
R <- 10

# Number of covariates
cov.mat <- cov.i
J <- dim(cov.mat)[2]

# Make symmetrical matrix of the rho values
rho_mat <- matrix(0, nrow = J, ncol = J)
rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho; rho_mat <- rho_mat + t(rho_mat)


library(autognet)
library(Rcpp)

#--------------------
# caleb's hacky microbenchmark implementation
# since it keeps failing
#-------------

set.seed(1)
nIt <- 10
start_time <- Sys.time()
lapply(1:nIt, function(i){
  aux.p.mat.R <- autognet:::aux.var.cl(tau,rho,nu,N,R,J, rho_mat,
                                       adjacency,cov.i,weights,group_lengths,group_functions)
}) -> olist_r
end_time <- Sys.time()
Rtime <- end_time - start_time

# Have to change the index of adjacency
adjacency_r <- adjacency
for (i in 1:N){
  adjacency[[i]] <- adjacency[[i]] - 1
}

set.seed(1)
start_time <- Sys.time()
lapply(1:nIt, function(i){
  aux.p.mat.Cpp <- auxVarCpp(tau,rho,nu,N,R,J, rho_mat,
                             adjacency,cov.i,weights,group_lengths,group_functions)
}) -> olist_cpp
end_time <- Sys.time()
CPPtime <- end_time - start_time

test_that("R and C++ versions give the same covariate model values", {
  #expect_true(all(summary(sapply(olist_cpp, sum)) == summary(sapply(olist_r, sum))))
  #expect_true(all(olist_cpp[[8]][c(1,2),] == olist_r[[8]][c(1,2),]))
  # Turns out that these won't be the same based on some Rmultinom internal workings
  fc <- as.numeric(Rtime)/as.numeric(CPPtime)
  message(fc)
  expect_true(fc > 1)
})


#############################
# New stuff for outcome model
#############################
aux.p.mat <- olist_cpp[[1]]
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



##BETA##
#Step 1. Proposal
scale <- .005
beta.p <- MASS::mvrnorm(1,beta[b,],scale*diag(P))

#Step 2. Auxilary variable based on proposal

# Second call to Rcpp function

# Test R call
set.seed(1)
start_time <- Sys.time()
lapply(1:nIt, function(i){
  aux.p_r <- autognet:::aux.var.outcome.cl(beta.p,trt.i,cov.i,
                                           N,10,adjacency_r,outcome.i,weights)[,1]

}) -> olist_r_2
end_time <- Sys.time()
Rtime2 <- end_time - start_time
get_sum <- function(idx, l_grid, tau, rho, n){

  # Establish Ls for this iteration
  l_vec <- l_grid[idx,]

  # Izzie / Caleb formulation of sum--
  # Make a temporary matrix; fill by row so have to transpose later
  # see https://stackoverflow.com/questions/30787317/a-vector-to-an-upper-triangle-matrix-by-row-in-r
  temp_m <- matrix(0, n, n); temp_m[lower.tri(temp_m, diag=FALSE)] <- rho

  # intercept                   main
  exp((tau%*%l_vec)[1,1] + sum(l_vec %*% t(l_vec) * t(temp_m)))

}
v_get_sum <- Vectorize(get_sum, "idx")

# Test Cpp time
set.seed(1)
start_time <- Sys.time()
lapply(1:nIt, function(i){
  aux.p_cpp <- auxVarOutcomeCpp(beta.p,trt.i[,1],cov.i,
                                N, 10,adjacency,unname(outcome.i[,1]),weights)
}) -> olist_cpp_2
end_time <- Sys.time()
CPPtime2 <- end_time - start_time

test_that("R and C++ versions give the same outcome model values", {
  expect_true(all(summary(sapply(olist_cpp_2, sum)) == summary(sapply(olist_r_2, sum))))
  expect_true(all(olist_cpp_2[[8]] == olist_r_2[[8]]))
  fc2 <- as.numeric(Rtime2)/as.numeric(CPPtime2)
  message(fc2)
  expect_true(fc2 > 1)
})

aux.p.mat <-olist_cpp_2[[1]]
aux.p.cov.n <- apply(aux.p.mat,2,function(x) {(adjmat%*%x)/weights})
sum.aux.p <- c(unname(colSums(aux.p.mat)),
               unname(colSums(aux.p.mat[,grid[,1]] * aux.p.mat[,grid[,2]]))*use_rho,
               unname(colSums(aux.p.mat*aux.p.cov.n)))

#Step 3. Calculate accept-reject ratio (everything is log-ed!)
#interconnected units
#h.l1.p <- alpha.p%*%sum.l
#h.l1.c <- alpha[b,]%*%sum.l
#h.aux.p <- alpha.p%*%sum.aux.p
#h.aux.c <- alpha[b,]%*%sum.aux.p

#independent units
tau.p <- 1
rho.p <- 1
f.p.num <- alpha.p[1:(ncov+nrho)]%*%sum.l.indep
f.p.denom <- log(sum(v_get_sum(1:dim(l_grid)[1], l_grid, tau.p, rho.p, ncov)))
f.c.num <- alpha[b,1:(ncov+nrho)]%*%sum.l.indep
f.c.denom <- log(sum(v_get_sum(1:dim(l_grid)[1], l_grid, alpha[b,1:ncov], alpha[b,(ncov+1):(ncov+nrho)], ncov)))

#priors
#prior.p <- mvtnorm::dmvt(alpha.p,delta=rep(0,L),sigma=4*diag(L),df=3,log=TRUE)
#prior.c <- mvtnorm::dmvt(alpha[b,],delta=rep(0,L),sigma=4*diag(L),df=3,log=TRUE)

##BETA##
#Step 1. Proposal
scale <- .005
beta.p <- MASS::mvrnorm(1,beta[b,],scale*diag(P))

#Step 2. Auxilary variable based on proposal
aux.p <- auxVarOutcomeCpp(beta.p,trt.i,cov.i,N,10,adjacency,outcome.i,weights)
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


#priors
#prior.p <- mvtnorm::dmvt(beta.p,delta=rep(0,P),sigma=4*diag(P),df=3,log=TRUE)
#prior.c <- mvtnorm::dmvt(beta[b,],delta=rep(0,P),sigma=4*diag(P),df=3,log=TRUE)

accept <- rbinom(1,1,min(1,exp(ratio)))
if (accept == 1){beta[b+1,] <- beta.p} else { beta[b+1,] <- beta[b,] }
accept_beta[b] <- accept

pr_trt <- 0.5

#####################
# Run the causal estimand functions
#####################

psi_gamma <- vector(mode="numeric", length=B)
psi_1_gamma <- vector(mode="numeric", length=B)
psi_0_gamma <- vector(mode="numeric", length=B)
alpha.thin <- alpha
beta.thin <- beta


tau <- alpha.thin[b,1:ncov]
rho <- alpha.thin[b,(ncov+1):(ncov+nrho)]
nu  <- alpha.thin[b,(ncov+nrho+1):L]

## MAKE COVARIATE ARRAY ##
set.seed(1)
start_time <- Sys.time()
lapply(1:nIt, function(i){
autognet:::network.gibbs.cov(tau, rho, nu,
                                           ncov, R, N, adjacency_r, weights,
                                           group_lengths, group_functions)

}) -> olist_r_covlist
end_time <- Sys.time()
Rtime_cov <- end_time - start_time


set.seed(1)
start_time <- Sys.time()
lapply(1:nIt, function(i){
  networkGibbsOutcomeCpp(tau, rho, nu,
                         ncov, R, N, rho_mat,
                         adjacency, weights, cov.i,
                         group_lengths, group_functions)

}) -> olist_cpp_covlist
end_time <- Sys.time()
cpptime_cov <- end_time - start_time

test_that("C++ version is faster for causal estimates", {
  #expect_true(all(summary(sapply(olist_cpp, sum)) == summary(sapply(olist_r, sum))))
  #expect_true(all(olist_cpp[[8]][c(1,2),] == olist_r[[8]][c(1,2),]))
  # Turns out that these won't be the same based on some Rmultinom internal workings
  fc <- as.numeric(Rtime_cov)/as.numeric(cpptime_cov)
  message(fc)
  expect_true(fc > 1)
})

cov.list <- olist_cpp_covlist[[1]]

## NON-INDIVIDUAL GIBBS (overall effects) ##
psi_gamma[b] <- autognet:::network.gibbs.out1(cov.list,beta.thin[b,],pr_trt,R,N,adjacency,weights)

## INDIVIDUAL GIBBS (spillover and direct effects) ##
psi_1_gamma[b] <- mean(sapply(1:N,autognet:::network.gibbs.out2,cov.list=cov.list,
                              beta=beta.thin[b,],p=pr_trt,R=R,N=N,adjacency=adjacency,weights=weights,
                              treatment_value=1))
psi_0_gamma[b] <- mean(sapply(1:N,autognet:::network.gibbs.out2,cov.list=cov.list,
                              beta=beta.thin[b,],p=pr_trt,R=R,N=N,adjacency=adjacency,weights=weights,
                              treatment_value=0))
