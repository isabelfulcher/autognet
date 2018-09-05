
context("C++ version of aux.var.cl performs well")
set.seed(14651)

# Import simulated data frame from IF
rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
adjmat <- rdsIn[[1]]
data <- rdsIn[[2]]

# Preprocess covariates
process_covariate <- autognet:::.covariate_process(data[,-c(1,2)]) # categorical --> binary
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

start_time <- Sys.time()
for(i in 1:100){
  aux.p.mat.Cpp <- auxVarCpp(tau,rho,nu,N,R,J, rho_mat,
                             adjacency,cov.i,weights,group_lengths,group_functions)
}
end_time <- Sys.time()
CPPtime <- end_time - start_time

start_time <- Sys.time()
for(i in 1:100){
aux.p.mat.R <- autognet:::aux.var.cl(tau,rho,nu,N,R,J, rho_mat,
                                     adjacency,cov.i,weights,group_lengths,group_functions)
}
end_time <- Sys.time()
Rtime <- end_time - start_time

as.numeric(Rtime)/as.numeric(CPPtime)

