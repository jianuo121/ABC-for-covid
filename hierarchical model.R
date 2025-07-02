################ run model ###############
run_model <- function(mu, gamma, tau, N1, N2) {
  
  # random effects for each studies
  delta <- rnorm(length(N1), 0, tau)
  
  logit_p1 <- mu + delta # implies logit(prob1)~rnorm(mu,tau)
  
  p1 <- 1 / (1 + exp(-logit_p1)) # inverse logit
  n1_star <- rbinom(length(N1), N1, p1) # simulation
  
  logit_p2 <- mu + delta + gamma #gamma is log(OR)
  p2 <- 1 / (1 + exp(-logit_p2))
  n2_star <- rbinom(length(N2), N2, p2)
  
  return(list(n1_star = n1_star, n2_star = n2_star))
}


############ weighted distance #############
weighted_distance <- function(p, q, N,tau) {
  
  se <- sqrt(p * (1 - p) / N)
  w<-1/(se^2+tau)
  normalized_w<-w/sum(w)
  dist <- sum((p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)))*normalized_w)
  
  return(dist)
}


################# some fixed values setting ##############
# given N1 N2
N1 <- bin_data_anonymised$N1
N2 <- bin_data_anonymised$N2

# number of accepted particles
N <- 3000 

# observed p1 and p2
p1_hat <- bin_data_anonymised$n1 / N1
p2_hat <- bin_data_anonymised$n2 / N2




################## ABC Alogrithm ###################

run_abc <- function(N) {
  # presample to choose epsilon
  distances_presample <- numeric(1000)
  for (k in 1:1000) {
    
    # Sampling from prior
    mu <- rnorm(1, -2, 2)
    gamma <- rnorm(1, 1, 1)   
    tau <- runif(1, 0.1, 1)  
    
    
    # simulation data
    D_star <- run_model(mu, gamma, tau, N1, N2)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    dist_arm1 <- weighted_distance(p1_hat, p1_star,N1,tau)
    dist_arm2 <- weighted_distance(p2_hat, p2_star,N2,tau)
    #dist_arm3<-weighted_distance(pmin(pmax(p1_hat+p2_hat,eps),1-eps), pmin(pmax(p1_star+p2_star,eps),1-eps), N1+N2,tau)
    
    distances_presample[k] <- dist_arm1 + dist_arm2 #+dist_arm3
  }
  # pick epsilon: 1th quantile of distance_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  
  accepted <- list(mu = c(), gamma = c(), tau = c(), distance = c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
   while(i <= N) {
    
    mu <- rnorm(1, -2, 1)
    gamma <- rnorm(1, 1, 1)   
    tau <- runif(1, 0.1, 1)  
    
    D_star <- run_model(mu, gamma, tau, N1, N2)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    
    dist_arm1 <- weighted_distance(p1_hat, p1_star,N1,tau)
    dist_arm2 <- weighted_distance(p2_hat, p2_star,N2,tau)
    #dist_arm3<-weighted_distance(pmin(pmax(p1_hat+p2_hat,eps),1-eps), pmin(pmax(p1_star+p2_star,eps),1-eps), N1+N2,tau)
    distance<-dist_arm1+dist_arm2#+dist_arm3
    
    
    if (distance < epsilon) {
      accepted$mu <- c(accepted$mu, mu)
      accepted$gamma <- c(accepted$gamma, gamma)
      accepted$tau <- c(accepted$tau, tau)
      accepted$distance <- c(accepted$distance, distance)
      
      i <- i + 1
    }
   
     j <- j + 1
     acc_rate <- i / j  # acceptance rate
     cat("current acceptance rate = ", round(acc_rate * 100, 2), "%\r")
   
     }
  
  return(accepted)
}

################# checking ################   
####### method 1 
library(ggplot2)

outputs<-run_abc(N)
mu_post <- median(outputs$mu)
gamma_post <- median(outputs$gamma)
tau_post <- median(outputs$tau)

sim_result <- run_model(mu = mu_post, gamma = gamma_post, tau = tau_post, 
                        N1 = N1, N2 = N2)

n1_sim <- sim_result$n1_star
n2_sim <- sim_result$n2_star

study <- 1:length(N1)
df_prob <- data.frame(
  study = rep(study, 4),
  prob = c(p1_hat, n1_sim / N1, p2_hat, n2_sim / N2),
  group = rep(c("p1_obs", "p1_sim", "p2_obs", "p2_sim"), each = length(study))
)

ggplot(df_prob, aes(x = study, y = prob, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Observed vs Simulated Infection probabilities", y = "Infection probability") +
  theme_minimal()

####### method 2
library(ggplot2)
library(dplyr)
library(tidyr)

# 1000 random samples of simulation results from all accepted parameter combinations

n_draws <- 1000
sim_list <- replicate(n_draws, {
  idx <- sample(1:length(outputs$mu), 1)
  params <- list(
    mu = outputs$mu[idx],
    gamma = outputs$gamma[idx],
    tau = outputs$tau[idx]
  )
  run_model(params$mu, params$gamma, params$tau, N1, N2)
}, simplify = FALSE)


sim_matrix_p1 <- sapply(sim_list, function(x) x$n1_star / N1)  
sim_matrix_p2 <- sapply(sim_list, function(x) x$n2_star / N2)

# calculate 2.5% and 97.5% quantiles
calc_quantiles <- function(matrix) {
  data.frame(
    lower = apply(matrix, 1, quantile, 0.025),
    upper = apply(matrix, 1, quantile, 0.975)
  )
}
q_1 <- calc_quantiles(sim_matrix_p1)
q_2 <- calc_quantiles(sim_matrix_p2)

# construct a dataframe for ggplot
study <- 1:length(N1)
df_plot <- data.frame(
  study = study,
  p1_obs = p1_hat,
  p2_obs = p2_hat,
  p1_lower = q_1$lower,
  p1_upper = q_1$upper,
  p2_lower = q_2$lower,
  p2_upper = q_2$upper
)

# ggplot with credible band
ggplot(df_plot) +
  # vaccine 1
  geom_ribbon(aes(x = study, ymin = p1_lower, ymax = p1_upper), fill = "red", alpha = 0.2) +
  geom_point(aes(x = study, y = p1_obs), color = "red", size = 2) +
  
  # vaccine 2
  geom_ribbon(aes(x = study, ymin = p2_lower, ymax = p2_upper), fill = "blue", alpha = 0.2) +
  geom_point(aes(x = study, y = p2_obs), color = "blue", size = 2) +
  
  labs(title = "Observed vs Posterior Simulated Infection Probabilities",
       y = "Infection probability", x = "Study") +
  theme_minimal()






