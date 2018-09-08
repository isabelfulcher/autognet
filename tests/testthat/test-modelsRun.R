
context("Verify autog computation works for covariate model")
set.seed(14651)

# Use sample data baked into package
rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
adjmat <- rdsIn[[1]]
data <- rdsIn[[2]]
treatment <- "treatment"
outcome <- "outcome"
B = 10
R = 5

test_that("Supplying multiple values in the seed vector yields multiple chains", {

  # Run the covariate model with only one chain
  mod_1chain <- agcParam(data, treatment, outcome, adjmat,
                         B = B, R = R, seed = c(1))

  # Run the covariate model with two chains (in parallel)
  mod_2chain <- agcParam(data, treatment, outcome, adjmat,
                         B = B, R = R, seed = c(1, 2))

  # Run the covariate model with three chains (still in parallel)
  mod_3chain <- agcParam(data, treatment, outcome, adjmat,
                         B = B, R = R, seed = c(1, 2, 12))

  # Verify outcomes of the function calls work for the numbers of chains anticipated
  expect_true(length(mod_1chain) == 1)
  expect_true(length(mod_2chain) == 2)
  expect_true(length(mod_3chain) == 3)
  expect_true(all(names(mod_3chain) == c("chain1", "chain2", "chain3")))
})

context("Verify autog computation works for outcome model")
mod <- agcParam(data, treatment, outcome, adjmat,
                B = 100, R = R, seed = c(1, 2))
effects <-  agcEffect(mod, burnin = 1, thin = 0.2, treatment_allocation = 0.5, R = 10, seed = 1)

test_that("Outcome model did something", {

  expect_true(sum(!is.na(effects[["chain1"]])) == 60)
  expect_true(sum(!is.na(effects[["chain2"]])) == 60)

})