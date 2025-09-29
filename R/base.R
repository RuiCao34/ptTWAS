function_on_scalar_regression <- function(Y, G, U, id, bt, weights, t_vec){
  set.seed(1)
  n = length(Y)
  n_s = length(unique(id))

  chunks = split(id, id)
  chunks_size = lapply(chunks, length)
  chunks_seq = lapply(chunks, function(x){
    which(id %in% x)
  })

  X = do.call('rbind', lapply(1:n_s, function(indi){
    # print(indi)
    indi_name = names(chunks_seq)[indi]
    kronecker(bt[chunks_seq[[indi]],], cbind(G, U, 1)[indi_name,,drop = F])
  } ))

  weights_sqrt = scale(sqrt(weights), center = F)
  X_weighted = sweep(X, 1, weights_sqrt, "*")
  Y_weighted = Y * weights_sqrt

  group = c(rep(1:(ncol(G)+ncol(U)+1), ncol(bt)))

  pf = rep(0, ncol(G)+ncol(U)+1)
  pf[1:ncol(G)] = 1

  model_indep_cv = gglasso::cv.gglasso(x = X_weighted, y = Y_weighted, group = group, pf = pf)
  model_indep = gglasso::gglasso(x = X_weighted, y = Y_weighted, lambda = model_indep_cv$lambda.min, group = group, pf = pf)

  resid_glmnet_indep_weighted = Y - predict(model_indep, X)
  suppressWarnings(
    {model_lmer = lme4::lmer(resid_glmnet_indep_weighted ~ -1 + (1|id))}
  )
  sigma2_r_hat = lme4::VarCorr(model_lmer)$`id`
  sigma2_error_hat = attr(lme4::VarCorr(model_lmer), "sc")
  Sigma_a_sqrtinverse_block = lapply(1:length(chunks_size), function(i){
    # print(i)
    chunks_size_i = chunks_size[[i]]
    eigen_decomp = eigen( diag(1/sqrt(weights[chunks_seq[[i]]]))%*% matrix(sigma2_r_hat, nrow = chunks_size_i, ncol = chunks_size_i)%*% diag(1/sqrt(weights[chunks_seq[[i]]])) + sigma2_error_hat * diag(1/weights[chunks_seq[[i]] ], chunks_size_i))
    return(eigen_decomp$vectors %*% diag(eigen_decomp$values^(-1/2))  %*% t(eigen_decomp$vectors) )
  })
  # Sigma_a_sqrtinverse_block2 = eigen_variance_hat_block2$vectors %*% diag(eigen_variance_hat_block2$values^(-1/2))  %*% t(eigen_variance_hat_block2$vectors)
  newX_list = lapply(1:length(unique(id)), function(indi){
    # print(indi)
    Sigma_a_sqrtinverse_block[[indi]] %*% X[chunks_seq[[indi]],]
  })
  newX = do.call("rbind", newX_list)
  newY_list = lapply(1:length(unique(id)), function(indi){
    Sigma_a_sqrtinverse_block[[indi]] %*% Y[chunks_seq[[indi]]]
  })
  newY = do.call("rbind", newY_list)
  model_glmnet2_final_cv = cv.gglasso.modified.v19(x = newX, y = newY, group = group, t_vec = t_vec, pf = pf, nsnps = ncol(G), Q = ncol(bt))
  model_glmnet2_final = gglasso::gglasso(x = newX, y = newY, lambda = model_glmnet2_final_cv$lambda.min, group = group, pf = pf)
  model_glmnet2_final$r2_cv_intervals = model_glmnet2_final_cv$r2_cv_intervals
  model_glmnet2_final$r2_cv = model_glmnet2_final_cv$r2_cv
  model_glmnet2_final$h2_cv = model_glmnet2_final_cv$h2_cv
  model_glmnet2_final$sigma2_r_hat = sigma2_r_hat
  model_glmnet2_final$sigma2_error_hat = sigma2_error_hat

  return(model_glmnet2_final)
}

cv.gglasso.modified.v19 = function(x, y, group, lambda = NULL, pred.loss = c("misclass", "loss", "L1", "L2"), nfolds = 5, foldid, t_vec, delta, nsnps, Q, ...){
  if (missing(pred.loss))
    pred.loss <- "default"
  else pred.loss <- match.arg(pred.loss)
  N <- nrow(x)
  y <- drop(y)
  if (missing(delta))
    delta <- 1
  if (delta < 0)
    stop("delta must be non-negtive")
  gglasso.object <- gglasso::gglasso(x, y, group, lambda = lambda, delta = delta) ########### modified
  lambda <- gglasso.object$lambda
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = N))
  else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))

  r2_matrix = array(dim = c(length(lambda), nfolds, 10)) ########### modified
  h2_vec = r2_vec = matrix(0, ncol = nfolds, nrow = length(lambda)) ########### modified
  intervals = (0:9)/10 ########### modified
  snp_seq = rep(c(rep(T, nsnps), rep(F, ncol(x)/Q-nsnps)), Q)########### modified

  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    outlist[[i]] <- gglasso::gglasso(x = x[!which, , drop = FALSE],
                            y = y_sub, group = group, lambda = lambda, delta = delta)  ########### modified

    r2_vec[,i] <- (cor(predict(outlist[[i]], newx = x[which, , drop = FALSE]), y[which] ))^2  ########### modified
    h2_vec[,i] <- (cor(x[which, snp_seq, drop = FALSE] %*% outlist[[i]]$beta[snp_seq,], y[which] ))^2########### modified
    t_groups <- findInterval(t_vec[which],intervals)   ########### modified
    for(j in 1:10){   ########### modified
      if(sum(t_groups==j) ==0){ ########### modified
        r2_matrix[,i,j] = 0
      }else{
        r2_matrix[,i,j] = apply(predict(outlist[[i]], newx = x[which, , drop = FALSE][which(t_groups==j),]),
                                2,
                                function(x){
                                  r2 = (cor(x, y[which][which(t_groups==j)] ))^2
                                  r2 = ifelse(is.na(r2), 0, r2)
                                  r2
                                }  ) ########### modified
      }
    }
    r2_vec[is.na(r2_vec)] = 0  ########### modified
    h2_vec[is.na(h2_vec)] = 0  ########### modified
  }
  fun <- paste("cv", class(gglasso.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid,
                               pred.loss, delta))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm +
                cvsd, cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.gglasso"
  obj$r2_cv_intervals = colMeans(r2_matrix[which(lambda == lamin$lambda.min),,]) ########### modified
  obj$r2_cv = r2_vec[which(lambda == lamin$lambda.min),]
  obj$h2_cv = h2_vec[which(lambda == lamin$lambda.min),]
  obj
}

getmin <- function(lambda, cvm, cvsd){
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}
