# pt-TWAS
![](/Fig/FlowChart.png)
RECKMON (REluCtant Kernel-based Modeling On Non-linearity) is a reluctant modeling framework tailored for high-dimensional genomic prediction. RECKMON reluctantly models main, second-order, and higher-order effects, leveraging an efficient kernel-based feature selection method to identify non-linear predictors with computational efficacy. 



# Installation
devtools::install_github("RuiCao34/ptTWAS")

# Example code

library(ptTWAS)

library(polynom)

library(gglasso)

library(ggplot2)

set.seed(123)

## Initialize parameters
n_s = 1e2

n_c = 1e1

n_step2 = 1e4

p = 2e1

p_true = 5e0

p_null = p - p_true

p_covar = 5

maf = 0.4

beta = rep(0.2,2)

beta_covar = 0

sigma2_error = 1

sigma2_r = 0

alpha1 = c(0.2,0.2)

alpha1_limit = c(0,1)

knot = NULL

t0 = 1

ge_baseline = c(0.5,1)

M = 2 # true gene expression curve is linear

## Initialize orthonormal splines

order = 4

ortho_splines = ortho_splines2 = matrix(0, nrow = order, ncol = order)

ortho_splines[,1] = c(1/sqrt(t0),0,0,0)

for(i in 2:order){

  e_i = ortho_splines[,i-1]
  
  h_ip1 = rep(0,order)
  
  h_ip1[i] = 1
  
  e_ip1 = h_ip1 - Reduce("+", lapply(1:(i-1), function(col){
  
    integral(polynomial(ortho_splines[,col]) * polynomial(h_ip1), limits = c(0,t0)) /
    
      integral(polynomial(ortho_splines[,col]) * polynomial(ortho_splines[,col]), limits = c(0,t0)) *
      
      ortho_splines[,col]
      
  }))
  
  ortho_splines[,i] = e_ip1 / sqrt(integral(polynomial(e_ip1) * polynomial(e_ip1), limits = c(0,t0)))
  
}

ortho_splines = ortho_splines[1:M,1:M]

## Stage 1 regression
### Define coefficients
beta_list = vector('list', length = 2)

beta_list[[1]] = c(rep(beta[1], p_true), rep(0, p_null), rep(0,p_covar)) # the beta of y coordinate of the node

beta_list[[2]] = c(rep(0, p_true), rep(beta[2], p_true), rep(0, p_null-p_true), rep(beta_covar,p_covar))

B_matrix = do.call("rbind", lapply(1:2, function(m){

  return(beta_list[[m]])
  
} ))

### Initialize genotype and covariate matrix
G = matrix(rbinom(n_s*p, size = 2, prob = maf), nrow = n_s)

U = matrix(rnorm(n_s*p_covar), nrow = n_s)

GU = cbind(G, U)

### Define individual gene expression function
ge_curve = X = t_vec = Y = vector("list", length = n_s)

for(subject in 1:n_s){

  ge_curve[[subject]] = as.vector(ge_baseline + ortho_splines %*% t(GU[subject,,drop = F]%*%t(B_matrix)))
  
  t_vec[[subject]] = runif(n_c, 0, 1)
  
  t_order_matrix = do.call("rbind", lapply(t_vec[[subject]], function(t) c(1,t,t^2,t^3)))[,1:M]
  
  Y[[subject]] = t_order_matrix %*% ge_curve[[subject]] + rnorm(n_c, sd = sqrt(sigma2_r)) + rnorm(n_c, sd = sqrt(sigma2_error))
  
  X[[subject]] = kronecker(X = do.call("cbind", lapply(0:(M-1), function(x) t_vec[[subject]]^x)) %*% ortho_splines, Y = cbind(1,GU)[subject,,drop = F])
  
}

ge_curve = do.call("rbind", ge_curve)

Y = unlist(Y)

t_vec = unlist(t_vec)

id = as.factor(rep(1:n_s, each = n_c))

X = do.call("rbind", X)

rownames(G) = rownames(U) = 1:n_s

weights = rep(1, nrow(X))

model_stage1 = ptTWAS_stage1(Y, G, U, id, t_vec, weights)

## Stage 2 using individual data
G_til = matrix(rbinom(n_step2*p, size = 2, prob = maf), nrow = n_step2)

U_til = matrix(rnorm(n_step2*p_covar), nrow = n_step2)

GU_til = cbind(G_til, U_til)

ge_curve = vector("list", length = n_step2)

mu = Z = rep(0, n_step2)

for(subject in 1:n_step2){

  ge_curve[[subject]] = as.vector(ge_baseline + ortho_splines %*% t(GU_til[subject,,drop = F]%*%t(B_matrix)))
  
  mu[subject] = integral( polynomial(ge_curve[[subject]]) * polynomial(ortho_splines %*% alpha1), limits = alpha1_limit )
  
}

mu = scale(mu, center = T, scale = F)

Z = sapply(mu, function(mu_i){

  rbinom(n = 1, size = 1, prob = exp(mu_i)/(1+exp(mu_i)))
  
})

### Hypothesis tests
model_stage2 = ptTWAS_stage2_individual(Z, G_til, model_stage1, outcome = "binary")

model_stage2$Global_pval # global p-value

model_stage2$nonflat_pval # non-flat p-value

### Simultaneous Confidence Band Inference

time_intervals = (0:1000)/1000

tilde_bqt = do.call("cbind", lapply(0:(model_stage1$n_splines-1), function(x) time_intervals^x)) %*% ortho_splines

alpha_hat = model_stage2$alpha_hat

alpha_var = model_stage2$alpha_var

var_bt = diag(tilde_bqt %*% alpha_var %*% t(tilde_bqt))

point_bt = tilde_bqt %*% alpha_hat

var_function = 0 # initial the denominator to scale

mean_genotype = colMeans(G_til)

B_hat = t(matrix(model_stage1$model$beta, ncol = M))

B_hat = B_hat[,1:model_stage1$p_snps]

for(i in 1:nrow(G_til)){

  diff_function = polynomial(ortho_splines %*% B_hat %*% (G_til[i,] - mean_genotype))
  
  var_function = var_function + diff_function*diff_function
  
}

var_function = var_function/(nrow(G_til)-1)

f = function(t){

  return(sqrt(predict(var_function, t)))
  
}

sd_integral = stats::integrate(f, lower = 0, upper = 1)$ value

df_ggplot <- data.frame(

  t = time_intervals,
  
  mean = point_bt*sd_integral,
  
  upper2 = (point_bt-qnorm(3e-2/2)*sqrt(var_bt))*sd_integral,
  
  lower2 = (point_bt+qnorm(3e-2/2)*sqrt(var_bt))*sd_integral
  
)

ggplot() +

  geom_line(data = df_ggplot, aes(x = t, y = upper2), color = "navy", linetype = "dashed", linewidth = 1) +
  
  geom_line(data = df_ggplot, aes(x = t, y = lower2), color = "navy", linetype = "dashed", linewidth = 1) +
  
  geom_line(data = df_ggplot, aes(x = t, y = mean), color = "black", linewidth = 2) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  
  labs(x = expression(paste('Pseudotime ',"(",italic(t),")")), y = expression(paste(hat(italic(alpha)),"(",italic(t),")"))) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  
  theme_classic()

## Stage 2 using summary data
GWAS_beta_hat = GWAS_beta_se = rep(0, p)

for(i in 1:p){

  GWAS_logistic_model = glm(Z~G_til[,i], family = "binomial")
  
  GWAS_beta_hat[i] = summary(GWAS_logistic_model)$coefficients[2,1]
  
  GWAS_beta_se[i] = summary(GWAS_logistic_model)$coefficients[2,2]
  
}

ptTWAS_stage2_summary(GWAS_beta_hat, GWAS_beta_se, GWAS_sample_size = length(Z), model_stage1 = model_stage1, ref_panel = G_til, varZ = var(Z))
