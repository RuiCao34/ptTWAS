#' @title Stage 1: Build Gene Imputation Model in pt-TWAS
#'
#' @description
#' This function performs the first stage regression (function-on-scalar) in pseudotime-dependent
#' Transcriptome-Wide Association Study (pt-TWAS). It fits the single cell gene expression, a function of pseudotime, using continuous
#' orthonormal splines, which are selected through a forward stepwise algorithm. A group lasso penalty term is imposed
#' onto the genotype coefficients to ensure SNP-wise sparsity.
#'
#' @details
#' The function first generates a set of orthonormal polynomial basis functions
#' on the interval [0, 1] using the Gram-Schmidt process. The order of the
#' polynomials is fixed at 4 (i.e., cubic polynomials).
#'
#' It then employs a forward selection strategy. It starts by fitting a predictive
#' model (via the helper function `function_on_scalar`) with only the first
#' basis spline. It progressively adds one spline at a time, refitting the model
#' and evaluating its performance using cross-validated R-squared ($R^2_{CV}$).
#' The process stops when adding a new spline does not improve the mean $R^2_{CV}$
#' by at least a predefined threshold (0.01). The function returns the model
#' corresponding to the optimal number of splines selected.
#'
#' @param Y A numeric vector of cell-level gene expression values, which can be log transformation of aggregated cell read counts.
#' @param G A numeric matrix of genotypes, with individuals in rows and SNPs
#'   in columns. Row names are the subject IDs.
#' @param U A numeric matrix of covariates to be adjusted for in the model,
#'   with individuals in rows and covariates in columns. Row names are the subject IDs.
#' @param id A cell-level vector of individual identifiers. The IDs should match those from input \code{G} and \code{U}.
#' @param t_vec A numeric vector representing the cell-level pseudotime values
#'   (scaled to be within [0, 1]).
#' @param weights A vector of cell-level weights, typically
#'   as the cell numbers per aggregation.
#'
#' @return A list containing the results of the model selection:
#'   \item{`model`}{The final model object
#'   that corresponds to the optimal number of splines.}
#'   \item{`n_splines`}{An integer representing the optimal number of basis
#'   splines selected.}
#'
#' @note The packages `gglasso`, `polynom` and `lme4` are dependencies for this function.
#' @export
ptTWAS_stage1 <- function(Y, G, U, id, t_vec, weights) {
  # require(c("polynom", "gglasso", "lme4"), quietly = T)
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
  bt_order4 = do.call("cbind", lapply(0:(order-1), function(x) t_vec^x)) %*% ortho_splines

  last_r2 = -Inf
  current_r2 = 0
  for(j in 1:ncol(bt_order4)){
    bt_temp = bt_order4[,1:j, drop=F]
    main_order4_noknot = function_on_scalar_regression(Y, G, U, id, bt = bt_temp, weights, t_vec)
    current_r2 = mean(main_order4_noknot$r2_cv)
    # print(main_order4_noknot$r2_cv)
    if(current_r2 >= last_r2 + 0.01 | j == 1){
      n_splines = j
      output_model = main_order4_noknot
      last_r2 = mean(main_order4_noknot$r2_cv)
    }else{
      break
    }
  }
  output = list(model = output_model, p_snps = ncol(G), p_covar = ncol(U), n_splines = n_splines)

  return(output)
}
