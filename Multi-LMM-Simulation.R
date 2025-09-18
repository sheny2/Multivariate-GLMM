library(doParallel)
library(foreach)
library(data.table)
library(dplyr)

source("DGP.R")
source("Multi-LMM.R")

# Parameters
H = 10                                # number of sites
m_hosp = sample(10:50, H, replace=TRUE) # patients per site (vector length H)
px = 6                               # number of covariates
p_bin = 3                            # number of binary covariates
py = 3                               # number of outcomes
beta = matrix(runif(px*py, 1, 10), nrow = px, ncol = py) # fixed effects
sigma_u = 0.2                        # site RE sd
sigma_e = 0.5                        # residual sd (per outcome)
rho = 0.3                            # exchangeable residual correlation across outcomes

# Define number of cores for parallel execution
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

N = 100


result_beta = matrix(nrow = (px*py), ncol = N)
rownames(result_beta) <- paste0("Beta", 1:(px*py))
result_sigma = matrix(nrow = 3, ncol = N)
rownames(result_sigma) <- c("sigma_u", "sigma_e", "rho")


# Run simulations in parallel
results <- foreach(k = 1:N, .packages = c("data.table", "dplyr", "MASS")) %dopar% {

  source("Multi-LMM.R")

  dat <- generate_data_mv(seed = sample(1:1e7, 1), H, m_hosp, px, p_bin, py, beta, sigma_u, sigma_e, rho)

  Y <- as.matrix(dat[, paste0("Y", 1:3)])
  X <- as.matrix(dat[, paste0("X", 1:6)])
  id.site <- dat$site

  fit <- peal_profile_exch_RI(Y = Y, X = X, id.site = id.site)

  fit <- tryCatch(
    {
      peal_profile_exch_RI(Y = Y, X = X, id.site = id.site)
    },
    error = function(e) {
      message("Error in fitting model: ", e$message)
      return(NULL)   # or return(list(b = NA, theta = NA, rho = NA, s2 = NA))
    },
    warning = function(w) {
      message("Warning in fitting model: ", w$message)
      invokeRestart("muffleWarning")  # keep going, ignore warning
    }
  )

  # Return the estimated parameters
  list(beta_res = fit$beta %>% matrix(nrow = px, ncol = py, byrow = F),
       sigma_res = c(sqrt(fit$sigma_u2),
                     sqrt(fit$sigma2),
                     fit$rho)
  )
}

# Store results into matrices
for (k in 1:N) {
  result_beta[, k] <- results[[k]]$beta_res
  result_sigma[, k] <- results[[k]]$sigma_res
}


stopCluster(cl)


# Sample true parameter values
true_beta <- c(beta)
true_sigma <- c(sigma_u, sigma_e, rho)


beta_df <- reshape2::melt(as.data.frame(result_beta))
colnames(beta_df) <- c("Simulation", "Estimate")
beta_df$Parameter <- rep(rownames(result_beta), ncol(result_beta))
beta_df$Parameter <- factor(beta_df$Parameter, levels = paste0("Beta", 1:(px*py)))

sigma_df <- reshape2::melt(as.data.frame(result_sigma))
colnames(sigma_df) <- c("Simulation", "Estimate")
sigma_df$Parameter <- rep(rownames(result_sigma), ncol(result_sigma))

beta_df$True_Value <- rep(true_beta, times = ncol(result_beta))
sigma_df$True_Value <- rep(true_sigma, times = ncol(result_sigma))


beta_est_bias <- beta_df %>% mutate(Bias = Estimate - True_Value) %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("beta_est_bias.png", beta_est_bias, width = 10, height = 4, dpi = 300)


sigma_est_bias <- sigma_df %>% mutate(Bias = Estimate - True_Value) %>%
  # filter(Parameter == "sigma_e") %>%
  ggplot(aes(x = Parameter, y = Bias)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(fill = "lightblue", alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sigma_est_bias.png", sigma_est_bias, width = 4, height = 4, dpi = 300)


