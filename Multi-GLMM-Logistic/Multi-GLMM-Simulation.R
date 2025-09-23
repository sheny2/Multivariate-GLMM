library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(glmmTMB)
library(TMB)
library(ggplot2)

source("Multi-GLMM-Logistic/MGLMM-Logistic-DGP.R")
TMB::compile("Multi-GLMM-Logistic/multivariate_logistic.cpp")

# parameters
px <- 4
beta <- matrix(runif(px * 2, -0.7, 0.7), nrow = px, ncol = 2)
H <- 60
m_hosp <- rep(1000, 60)
p_bin <- 1
sigma_a <- 0.7
sigma_b <- 0.3
rho <- 0.2

# parallel execution
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# prep results to be stored
N <- 100 # number of simulations
result_beta = matrix(nrow = (px*2), ncol = N) # manually store for two outcomes
rownames(result_beta) <- paste0("Beta", 1:(px*2))
result_sigma = matrix(nrow = 3, ncol = N)
rownames(result_sigma) <- c("sigma_a", "sigma_b", "rho")

results <- foreach(k = 1:N, .packages = c("dplyr", "MASS")) %dopar% {
  
  dat <- generate_log_data_mv(seed = sample(1:1e7, 1), H = H, m_hosp = m_hosp, px = px, p_bin = p_bin, 
                              beta = beta, sigma_a = sigma_a, sigma_b = sigma_b, rho = rho)
  
  Y <- as.matrix(dat[, paste0("Y", 1:2)])
  X <- as.matrix(dat[, paste0("X", 1:4)])
  h_id0 <- dat$site - 1 # c++ uses 0-based indexing
  
  H <- length(unique(dat$site))
  Ntrials <- matrix(1, nrow(Y), ncol(Y))
  p <- ncol(X)
  r <- ncol(Y)
  
  base::dyn.load(TMB::dynlib("Multi-GLMM-Logistic/multivariate_logistic"))
  
  parameters <- list(
    beta  = matrix(0, p, r),   
    U     = matrix(0, H, r),  
    rho   = rep(0, r*(r-1)/2), 
    sigma = rep(1, r))
  
  data <- list(X = X, Y = Y, Ntrials = Ntrials,
    hospital_id = as.integer(h_id0))
  
  obj <- TMB::MakeADFun(data = data,
    parameters = parameters,
    random = "U",                      
    DLL = "multivariate_logistic")
  
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
  rep <- TMB::sdreport(obj)
  fit <- summary(rep, "fixed")
  betas <- fit[rownames(fit) == "beta", "Estimate"]
  beta_mat <- matrix(betas, nrow = px, ncol = 2, byrow = FALSE)
  rho_est <- fit[rownames(fit) == "rho", "Estimate"] 
  
  # estimated parameters
  list(beta_res = beta_mat,
       sigma_res = c(fit[10,1],
                     fit[11,1],
                     rho_est)
  )
}

for (k in 1:N) {
  result_beta[, k] <- results[[k]]$beta_res
  result_sigma[, k] <- results[[k]]$sigma_res
}

stopCluster(cl)

true_beta <- c(beta)
true_sigma <- c(sigma_a, sigma_b, rho)

beta_df <- reshape2::melt(as.data.frame(result_beta))
colnames(beta_df) <- c("Simulation", "Estimate")
beta_df$Parameter <- rep(rownames(result_beta), ncol(result_beta))
beta_df$Parameter <- factor(beta_df$Parameter, levels = paste0("Beta", 1:(px*2)))

sigma_df <- reshape2::melt(as.data.frame(result_sigma))
colnames(sigma_df) <- c("Simulation", "Estimate")
sigma_df$Parameter <- rep(rownames(result_sigma), ncol(result_sigma))

beta_df$True_Value <- rep(true_beta, times = ncol(result_beta))
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma))

beta_est_bias <- beta_df %>% mutate(Bias = Estimate - True_Value) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "mediumvioletred", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()
ggsave("logistic_beta_est_bias.png", beta_est_bias, width = 10, height = 4, dpi = 300)

sigma_est_bias <- sigma_df %>% mutate(Bias = Estimate - True_Value) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "mediumvioletred", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()
ggsave("logistic_sigma_est_bias.png", sigma_est_bias, width = 4, height = 4, dpi = 300)

