# RECKMON
RECKMON (REluCtant Kernel-based Modeling On Non-linearity) is a reluctant modeling framework tailored for high-dimensional genomic prediction. RECKMON reluctantly models main, second-order, and higher-order effects, leveraging an efficient kernel-based feature selection method to identify non-linear predictors with computational efficacy. 


# Installation
devtools::install_github("RuiCao34/RECKMON")

# Example code

library(RECKMON)

set.seed(123)

n <- 100

p <- 10

x_data <- matrix(rnorm(n * p), nrow = n, ncol = p)

colnames(x_data) <- paste0("V", 1:p)

y_data <- 2 * x_data[, 1] - 1.5 * x_data[, 2] + x_data[, 3] * x_data[, 4] + sin(pi * x_data[, 5]) + rnorm(n, 0, 0.5)                          

model_fit <- main_RECKMON(x = x_data, y = y_data, poly_features = 5, step_gaussian = 1, cv = 3, mc.cores = 1)

cat("Selected interaction features:\n")

print(model_fit$poly_feature)

cat("Selected Gaussian kernel features:\n")

print(model_fit$gaussian_feature)

predict_RECKMON(newx = matrix( rnorm(p*2), nrow = 2), model_fit)
