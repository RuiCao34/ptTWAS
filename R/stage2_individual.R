#' @title Stage 2: Regression with GWAS Outcome
#'
#' @description
#' This function performs the second stage regression (scalar-on-function) in pseudotime-dependent
#' Transcriptome-Wide Association Study (pt-TWAS) using individual-level GWAS data. Still, it uses continuous orthonormal splines to model the
#' gene effect curve \eqn{\alpha}(t).
#'
#' @details
#' This function regress the GWAS outcome onto predicted gene expression curve (scalar-on-function regression).
#' Inferences are based on hypothesis tests (global, non-flat) and simultaneous confidence band.
#'
#'
#' @param Z A numeric vector of GWAS outcome.
#' @param G_til A numeric matrix of genotypes, with individuals in rows and SNPs
#'   in columns. Row names are the subject IDs.
#' @param model_stage1 The gene expression model fitted by \code{ptTWAS_stage1}.
#' @param outcome Either `binary` or `continuous`.
#' @param simul_band_cutoffs A numeric vector of alpha levels to test for the simultaneous confidence band.
#' @param boot_iter Number of bootstrap iterations.
#'
#' @return A list containing the results of the Stage 2 regression:
#'    \item{`alpha_hat`}{The estimate of vector \eqn{\alpha} under orthonormal basis splines.}
#'    \item{`alpha_var`}{The variance-covariance estimate of vector \eqn{\alpha} under orthonormal basis splines.}
#'    \item{`Global_pval`}{P-value for global test.}
#'    \item{`nonflat_pval`}{P-value for nonflat test.}
#'    \item{`gamma_til`}{The calibrated confidence factor.}
#' @export
ptTWAS_stage2_individual <- function(Z, G_til, model_stage1, outcome = "binary", simul_band_cutoffs = c(5e-2, 4e-2, 3e-2, 2e-2, 1e-2, 5e-5, 1e-3), boot_iter = 1e3){
  order = 4
  ortho_splines = matrix(0, nrow = order, ncol = order)
  ortho_splines[,1] = c(1,0,0,0)
  for(i in 2:order){
    e_i = ortho_splines[,i-1]
    h_ip1 = rep(0,order)
    h_ip1[i] = 1
    e_ip1 = h_ip1 - Reduce("+", lapply(1:(i-1), function(col){
      polynom::integral(polynom::polynomial(ortho_splines[,col]) * polynom::polynomial(h_ip1), limits = c(0,1)) /
        integral(polynom::polynomial(ortho_splines[,col]) * polynom::polynomial(ortho_splines[,col]), limits = c(0,1)) *
        ortho_splines[,col]
    }))
    ortho_splines[,i] = e_ip1 / sqrt(polynom::integral(polynom::polynomial(e_ip1) * polynom::polynomial(e_ip1), limits = c(0,1)))
  }

  p_snps = model_stage1$p_snps
  M = model_stage1$n_splines
  ortho_splines = ortho_splines[1:M, 1:M]
  integral_btilbT = diag(1,M)
  B_hat = t(matrix(model_stage1$model$beta, ncol = M))
  B_hat = B_hat[,1:p_snps]
  W = t(B_hat) %*% integral_btilbT

  ## hypothesis test
  if(outcome == "binary"){
    logit_reg = glm(Z ~ G_til%*%W, family = binomial)
    alpha_hat = coef(logit_reg)[-1]
    alpha_var = vcov(logit_reg)[-1,-1,drop=F]
  }else{
    linear_reg = lm(Z ~ G_til%*%W)
    alpha_hat = coef(linear_reg)[-1]
    alpha_var = vcov(linear_reg)[-1,-1,drop=F]
  }
  Global_pval = 1- pchisq(t(alpha_hat) %*% solve(alpha_var) %*% alpha_hat, df = length(alpha_hat))
  if(M>1){
    nonflat_pval = 1- pchisq(t(alpha_hat[-1]) %*% solve(alpha_var[-1,-1]) %*% alpha_hat[-1], df = length(alpha_hat)-1)
  }else{
    nonflat_pval = NA
  }


  ## confidence band
  time_intervals = seq(0,1, length.out = 100)
  bt = do.call("cbind", lapply(0:(M-1), function(x) time_intervals^x)) %*% ortho_splines
  var_bt = diag(bt %*% alpha_var %*% t(bt))
  point_bt = bt %*% alpha_hat
  mydata = as.data.frame(cbind(Z, G_til %*% W))
  colnames(mydata)[1] = "Z"

  # The bootstrap statistic function. It re-fits the model on bootstrapped data
  # and checks if the new estimated curve lies within a candidate confidence band.
  bootstrap.fun <- function(data, alpha0) {
    if(outcome == "binary"){
      fit <- glm(Z ~. , data = data, family = "binomial")
    }else{
      fit <- lm(Z ~. , data = data)
    }
    # Handle cases where model fitting fails, e.g., due to singularity
    if (any(is.na(coef(fit)))) {
      return(rep(FALSE, length(time_intervals)))
    }
    alpha_hat_new = coef(fit)[-1]
    point_bt_new = bt %*% alpha_hat_new
    alpha = abs(qnorm(alpha0/2))
    # 'point_bt' and 'var_bt' are from the original model fit in the parent environment
    whether_in_band = (point_bt_new > point_bt - alpha*sqrt(var_bt)) & (point_bt_new < point_bt + alpha*sqrt(var_bt))
    return(whether_in_band)
  }

  # The random data generation function for parametric bootstrap.
  # 'mle' is the model object fitted to the original data.
  ran_gen <- function(data, mle){
    out <- data
    if(outcome == "binary"){
      # For glm (binary), fitted values are probabilities.
      probs <- fitted(mle)
      out$Z <- rbinom(n = length(probs), size = 1, prob = probs)
    }else{
      # For lm (continuous), fitted values are means. We also need the residual standard error.
      means <- fitted(mle)
      sigma_hat <- summary(mle)$sigma
      out$Z <- rnorm(n = length(means), mean = means, sd = sigma_hat)
    }
    out
  }

  # Select the correct model fit to pass to the bootstrap's 'mle' argument.
  if (outcome == "binary") {
    fit_model <- logit_reg
  } else {
    fit_model <- linear_reg
  }

  # Initialize gamma_til to NA. It will be updated if the loop finds a suitable value.
  gamma_til <- NA

  # Loop through candidate significance levels to find the one that gives >= 95% coverage.
  for(cutoff in simul_band_cutoffs){
    message(paste0("bootstrapping under modified band width factor ", round(abs(qnorm(cutoff/2)),2), "..." ))
    boot.sim <- boot::boot(data=mydata, sim = "parametric", statistic=bootstrap.fun,
                           R=boot_iter, ran.gen = ran_gen, mle = fit_model, alpha0 = cutoff)

    # Calculate the proportion of replicates where the *entire* curve was within the band.
    coverage_rate = mean(rowMeans(boot.sim$t, na.rm = TRUE) == 1, na.rm = TRUE)
    if(is.na(coverage_rate)) coverage_rate <- 0 # Handle cases where all bootstrap runs fail.

    if(coverage_rate > 0.95){
      gamma_til = abs(qnorm(cutoff/2))
      break # Exit loop once desired coverage is achieved.
    }
  }

  return(list('alpha_hat' = alpha_hat,
              'alpha_var' = alpha_var,
              'Global_pval' = Global_pval,
              'nonflat_pval' = nonflat_pval,
              'gamma_til' = gamma_til))
}
