
context("Covariate parsing for data frames")
set.seed(14651)

test_df <- data.frame(
  a = factor(c("black", "hispanic", "white", "white", "black", "white", "asian", "white"),
             levels = c("white", "asian", "black", "hispanic")),
  b = c(0,1,1,1,1,0,0,1),
  c =  factor(c("male", "female", "female", "male", "male", "male", "female", "other"),
              levels = c("female", "male", "other")),
  d = c(rep("hi", 7), "bye"),
  e = factor(c("black", "hispanic", "white", "white", "black", "white", "asian", "white"),
             levels = c("white", "asian", "black", "hispanic"))
)

expand_out <- autognet:::.covariate_process(test_df)

test_that("Covariate expand_out works", {
  expect_true(length(expand_out) == 4)
  expect_true(dim(expand_out[[1]])[2] == 10)
  expect_true(length(expand_out[[2]]) == 7)
  expect_true(length(expand_out[[3]]) == 5)
  expect_true(length(expand_out[[4]]) == 5)
})




