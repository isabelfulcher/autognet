library(testthat)
context("We can flexibly handle additional nu terms in the model")

context("Verify autog computation works for covariate model")

# Use sample data baked into package
rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
adjmat <- rdsIn[[1]]
data <- rdsIn[[2]]
treatment <- "treatment"
outcome <- "outcome"
B = 10
R = 5

set.seed(123)
yes_nu <- agcParam(data, treatment, outcome, adjmat,
                   B = B, R = R, seed = c(1), additional_nu = TRUE)


set.seed(123)
no_nu <- agcParam(data, treatment, outcome, adjmat,
                       B = B, R = R, seed = c(1))


tau.p <- c(1,2)
rho.p <- c(-.5)
nu.p.mat <- diag(c(.5,.5))
N <- 75

group_lengths <- c(1,1)
group_functions <- c(1,1)
cov.i <- data.matrix(data[,c(3,4)])

# Test two scenarios about whether additional_nu makes sense
J = 2
rho_mat <- matrix(0, nrow = J, ncol = J)
rho_mat[lower.tri(rho_mat, diag=FALSE)] <- rho.p
rho_mat <- rho_mat + t(rho_mat)
adjacency <- data.matrix(no_nu$adjmat)
weights <- apply(adjacency,1,sum) # number of neighbors for everyone

set.seed(134)
w_add_nu <- auxVarCpp(tau.p, rho.p, nu.p.mat, N, R, J=2, rho_mat,
          adjacency,cov.i,weights,group_lengths,group_functions,
          as.numeric(TRUE))
set.seed(134)
no_add_nu  <- auxVarCpp(tau.p, rho.p, nu.p.mat, N, R, J=2, rho_mat,
          adjacency,cov.i,weights,group_lengths,group_functions,
          as.numeric(FALSE))



test_that("Utilizing additional nu terms doesn't change anything", {
  expect_true(all(w_add_nu == no_add_nu))
})

test_that("Matrix setup with off-diag works as expected", {
  expect_true(all(w_add_nu == no_add_nu))

  # Now look at the
  aux.p.mat <- w_add_nu
  aux.p.cov.n <- apply(aux.p.mat,2,function(x) {(adjmat%*%x)/weights})
  mat_mult_version <- diag(matrix(as.vector(t(t(aux.p.mat) %*% aux.p.cov.n)), ncol = 2))

  aux.p.mat <- no_add_nu
  aux.p.cov.n <- apply(aux.p.mat,2,function(x) {(adjmat%*%x)/weights})
  vec_version <- unname(colSums(aux.p.mat*aux.p.cov.n))
  expect_true(all(round(vec_version,1) == round(mat_mult_version,1)))

})




