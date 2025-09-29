#' @title Stage 2: Regression with GWAS Outcome
#'
#' @description
#' This function performs the second stage regression (scalar-on-function) in pseudotime-dependent
#' Transcriptome-Wide Association Study (pt-TWAS) using summary-level GWAS data. Still, it uses continuous orthonormal splines to model the
#' gene effect curve \eqn{\alpha}(t).
#'
#' @details
#' This function regress the GWAS outcome onto predicted gene expression curve (scalar-on-function regression) using summary statistics.
#' Inferences are based on hypothesis tests (global, non-flat) only.
#'
#' @param GWAS_beta_hat GWAS regression coefficients.
#' @param GWAS_beta_se GWAS regression coefficient sds.
#' @param GWAS_sample_size GWAS sample size.
#' @param outcome Either `binary` or `continuous`.
#' @param model_stage1 The gene expression model fitted by \code{ptTWAS_stage1}.
#' @param ref_panel A reference panel.
#' @param control_case_ratio Control-case ratio. Only available when \code{outcome='binary'}.
#' @param af_vec Allele frequencies of SNPs. If \code{NULL}, estimated from the reference panel.
#' @param var_Z Variance of outcome \code{Z}. Only available when \code{outcome='continuous'}.
#'
#' @return A list containing the results of the Stage 2 regression:
#'    \item{`alpha_hat`}{The estimate of vector \eqn{\alpha} under orthonormal basis splines.}
#'    \item{`alpha_var`}{The variance-covariance estimate of vector \eqn{\alpha} under orthonormal basis splines.}
#'    \item{`Global_pval`}{P-value for global test.}
#'    \item{`nonflat_pval`}{P-value for nonflat test.}
#' @export
ptTWAS_stage2_summary <- function(GWAS_beta_hat, GWAS_beta_se, GWAS_sample_size, outcome = "binary", model_stage1, ref_panel, control_case_ratio = NULL, af_vec = NULL, var_Z = NULL){
  M = model_stage1$n_splines
  p_snps = model_stage1$p_snps
  B_hat = t(matrix(model_stage1$model$beta, ncol = M))
  B_hat = B_hat[,1:p_snps]

  if(is.null(af_vec)) af_vec = colMeans(ref_panel, na.rm = T) /2
  N_til = GWAS_sample_size


  ref_panel = scale(ref_panel, center = T, scale = F)

  # snp variance
  s2 = 2 * af_vec * (1-af_vec)

  # transform logistic to linear
  if(outcome == "continuous"){
    xi_hat = GWAS_beta_hat
    xi_se = GWAS_beta_se
  }else if(outcome == "binary"){
    xi_hat = control_case_ratio/(1+control_case_ratio)^2 * GWAS_beta_hat
    xi_se = control_case_ratio/(1+control_case_ratio)^2 * GWAS_beta_se
  }else{
    warning('GWAS model must be either linear or logistic')
  }


  scale_GTG_ref = t(ref_panel) %*% ref_panel / nrow(ref_panel)
  scale_BGTZ = B_hat %*% (xi_hat * s2) ## B_hat %*% t(X) %*% Y/nrow(X)
  # scale_BGTZ = B_hat %*% scale_GTG_ref %*% xi_hat ## equal to B_hat %*% t(X) %*% Y/nrow(X)
  if(is.null(var_Z)){
    scale_ZTZ = median(N_til * s2 * xi_se^2 + s2 * xi_hat^2)
  }else{
    scale_ZTZ = var_Z
  }

  scale_H = solve(B_hat %*% scale_GTG_ref %*% t(B_hat))
  sigma_til_square = ((scale_ZTZ - t(scale_BGTZ) %*% scale_H %*% scale_BGTZ)  / ((N_til-nrow(B_hat))/N_til))[1,1]
  alpha_hat = scale_H %*% scale_BGTZ
  alpha_var = sigma_til_square * scale_H / GWAS_sample_size

  if(outcome == "binary"){
    alpha_hat = alpha_hat /(control_case_ratio/(1+control_case_ratio)^2)
    alpha_var = alpha_var /(control_case_ratio/(1+control_case_ratio)^2)^2
  }

  Global_pval = 1 - pchisq(t(alpha_hat) %*% solve(alpha_var) %*% alpha_hat, df = length(alpha_hat))
  if(M>1){
    nonflat_pval = 1- pchisq(t(alpha_hat[-1]) %*% solve(alpha_var[-1,-1]) %*% alpha_hat[-1], df = length(alpha_hat)-1)
  }else{
    nonflat_pval = NA
  }

  return(list(alpha_hat = alpha_hat,
              alpha_var = alpha_var,
              Global_pval = Global_pval,
              nonflat_pval = nonflat_pval))
}
