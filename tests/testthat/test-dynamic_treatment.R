
context("Verify autog computation works for various treatment regimes when estimating the causal effects")
set.seed(1)
library(autognet)
# Use sample data baked into package
rdsIn <- readRDS(paste0(system.file('extdata',package='autognet'),"/agc-example.rds"))
adjmat <- rdsIn[[1]]
data <- rdsIn[[2]]
treatment <- "treatment"
outcome <- "outcome"
B = 10
R = 5

# Run the covariate model with only one chain
mod_1chain <- agcParam(data, treatment, outcome, adjmat,
                       B = B, R = R, seed = c(1))
set.seed(1)
dyn_effect_base1 <- autognet::agcEffect(mod_1chain)
set.seed(1)
dyn_effect_se <- agcEffect(mod_1chain, dynamic_single_edge = 2)
set.seed(1)
dyn_effect_base1_again <- agcEffect(mod_1chain)
set.seed(1)
dyn_effect_dat <- agcEffect(mod_1chain, subset = which(data$treatment == 1),
                            dynamic_single_edge = 2)
set.seed(1)
dyn_effect_dat_oo <- agcEffect(mod_1chain, subset = which(data$treatment == 1),
                            dynamic_single_edge = 2, overall_effects_only = 1)

