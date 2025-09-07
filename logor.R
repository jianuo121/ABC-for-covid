
############ weighted distance #############
weighted_distance_KL <- function(p, q,N,tau) {
  
  w<-1/(p*(1-p)/N + tau^2)
  w<-w/sum(w)
  
  dist <- sum((p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)))*w)
  
  return(dist)
}

weighted_distance_l1 <- function(m, n, tau, var) {
  w<-1/(var+tau^2)
  normalized_w<-w/sum(w)
  dist <- sum(abs(m-n)*normalized_w)/length(n)
  
  return(dist)
}
##################### distance ########################


distance_abs <- function(m, n) {
  
  dist <- abs(m-n)
  
  return(dist)
}
################# some fixed values setting ##############
# given N1 N2
N1 <- bin_data_anonymised$N1
N2 <- bin_data_anonymised$N2
se_obs=sqrt(lor_data_anonymised$vi)
# number of accepted particles
N <- 1000 

# observed p1 and p2
p1_hat <- bin_data_anonymised$n1 / N1
p2_hat <- bin_data_anonymised$n2 / N2




# define half-normal sampler
rhalfnorm <- function(n, sigma) {
  abs(rnorm(n, 0, sigma))
}

sample_sizes <- c(2000, 5000, 8000)

# prior sets
prior_sets <- list(
  A = list(tau = function() rhalfnorm(1, 0.5),
           gamma_mu = function() rnorm(1, 0, 4),
           g = function(n) rnorm(n, -4, 4)),
  B = list(tau = function() rhalfnorm(1, 0.25),
           gamma_mu = function() rnorm(1, 0, 2),
           g = function(n) rnorm(n, -4, 2)),
  C = list(tau = function() rhalfnorm(1, 0.25),
           gamma_mu = function() rnorm(1, 0, 1),
           g = function(n) rnorm(n, -4, 1))
)

# 存储结果
results_df <- data.frame(
  sample_size = integer(),
  prior_set = character(),
  epsilon = numeric(),
  mean_distance = numeric(),
  sd_distance = numeric(),
  stringsAsFactors = FALSE
)
set.seed(34)
for (n_presample in sample_sizes) {
  for (prior_name in names(prior_sets)) {
    
    cat("Running presample =", n_presample, "with prior =", prior_name, "\n")
    
    distances_presample <- numeric(n_presample)
    prior <- prior_sets[[prior_name]]
    
    for (k in 1:n_presample) {
      # sample parameters
      
      g <- prior$g(nrow(bin_data_anonymised))
      gamma_mu <- prior$gamma_mu()
      tau <- prior$tau()
      
      # simulate
      D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
      
      # calculate p1*, p2*
      p1_star <- D_star$n1_star / N1
      p2_star <- D_star$n2_star / N2
      
      p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
      p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
      
      # distances
      dist_arm1 <- weighted_distance_KL(p1_hat, p1_star, N1, tau)
      dist_arm3 <- weighted_distance_abs(D_star$logOR,
                                         lor_data_anonymised$yi,
                                         tau,
                                         lor_data_anonymised$vi)
      
      distances_presample[k] <- dist_arm1 + dist_arm3
    }
    
    epsilon <- quantile(distances_presample, 0.01)
    cat("Selected epsilon for", prior_name, "with", n_presample, ":", epsilon, "\n")
    
    # 保存到 data frame
    results_df <- rbind(results_df, data.frame(
      sample_size = n_presample,
      prior_set = prior_name,
      epsilon = epsilon,
      mean_distance = mean(distances_presample),
      sd_distance = sd(distances_presample)
    ))
  }
}


print(results_df)



################ run model ###############
run_model_logor <- function(g, gamma_mu, tau, N1, N2, se_obs) {
  

  # binary data
  logit_p2 <- g 
  p2 <- 1 / (1 + exp(-logit_p2))
  n2_star <- rbinom(length(N2), N2, p2)
  
  
  gamma <- rnorm(length(N1)+nrow(lor_data_anonymised), gamma_mu, tau)
  logit_p1 <- g  + gamma[1:length(N1)]
  p1 <- 1 / (1 + exp(-logit_p1)) # inverse logit
  n1_star <- rbinom(length(N1), N1, p1) # simulation
  # logOR data
  
  logOR <- rnorm(nrow(lor_data_anonymised), gamma[(length(N1)+1):(length(N1)+nrow(lor_data_anonymised))], se_obs)
  
  return(list(n1_star = n1_star, n2_star = n2_star, logOR=logOR))
  }






################## ABC Alogrithm ###################

run_abc_logor_1 <- function(N) {
   
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    # simulation data
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    dist_arm1 <- weighted_distance_KL(p1_hat, p1_star,N1,tau)
    dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    #dist_arm1 <- distance_abs(sum(bin_data_anonymised$n1)/sum(N1), sum(D_star$n1)/sum(N1))
    #dist_arm2 <- distance_abs(sum(bin_data_anonymised$n2)/sum(N2), sum(D_star$n1)/sum(N1))
    #dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    distances_presample[k] <- dist_arm1+dist_arm3
  }
  # pick epsilon: 1% quantile of distances_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  accepted <- list(g = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   gamma_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm2=c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    
    dist_arm1 <- weighted_distance_KL(p1_hat, p1_star,N1,tau)
    dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
   
    distance<-dist_arm1+dist_arm3
    
    if (distance <= epsilon) {
      accepted$g <- rbind(accepted$g, g) 
      accepted$gamma_mu <- c(accepted$gamma_mu, gamma_mu)
      accepted$tau <- c(accepted$tau, tau)
      accepted$distance <- c(accepted$distance, distance)
      accepted$dist_arm1<-c(accepted$dist_arm1, dist_arm1)
      accepted$dist_arm3<-c(accepted$dist_arm3, dist_arm3)
      i <- i + 1
    }
    
    j <- j + 1
    acc_rate <- i / j  # acceptance rate
    cat("current acceptance rate = ", round(acc_rate * 100, 5), "%\r")
    
  }
  
  return(accepted)
}

run_abc_logor_2 <- function(N) {
  
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    # simulation data
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    dist_arm1 <- weighted_distance_l1(p1_hat, p1_star,tau,p1_star*(1-p1_star)/N1)
    dist_arm3<- weighted_distance_l1(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    #dist_arm1 <- distance_abs(sum(bin_data_anonymised$n1)/sum(N1), sum(D_star$n1)/sum(N1))
    #dist_arm2 <- distance_abs(sum(bin_data_anonymised$n2)/sum(N2), sum(D_star$n1)/sum(N1))
    #dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    distances_presample[k] <- dist_arm1+dist_arm3
  }
  # pick epsilon: 1% quantile of distances_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  accepted <- list(g = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   gamma_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm3=c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    dist_arm1 <- weighted_distance_l1(p1_hat, p1_star,tau,p1_star*(1-p1_star)/N1)
    dist_arm3<- weighted_distance_l1(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    #dist_arm1 <- distance_abs(sum(bin_data_anonymised$n1)/sum(N1), sum(D_star$n1)/sum(N1))
    #dist_arm2 <- distance_abs(sum(bin_data_anonymised$n2)/sum(N2), sum(D_star$n1)/sum(N1))
    #dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    
    distance<-dist_arm1+dist_arm3
    
    if (distance <= epsilon) {
      accepted$g <- rbind(accepted$g, g) 
      accepted$gamma_mu <- c(accepted$gamma_mu, gamma_mu)
      accepted$tau <- c(accepted$tau, tau)
      accepted$distance <- c(accepted$distance, distance)
      accepted$dist_arm1 <- c(accepted$dist_arm1, dist_arm1)
      accepted$dist_arm3 <- c(accepted$dist_arm3, dist_arm3)
      i <- i + 1
    }
    
    j <- j + 1
    acc_rate <- i / j  # acceptance rate
    cat("current acceptance rate = ", round(acc_rate * 100, 5), "%\r")
    
  }
  
  return(accepted)
}
run_abc_logor_3<- function(N) {
  
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    # simulation data
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    w1<-1/(p1_hat*(1-p1_hat)/N + tau^2)
    w1<-w1/sum(w1)
    w3<-1/(se_obs^2 + tau^2)
    w3<-w3/sum(w3)
    dist_arm1 <- distance_abs(sum(w1*p1_hat), sum(w1*p1_star))
    dist_arm3<- distance_abs(sum(w3*D_star$logOR),sum(w3*lor_data_anonymised$yi))
    distances_presample[k] <- dist_arm1+dist_arm3
  }
  # pick epsilon: 1% quantile of distances_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  accepted <- list(g = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   gamma_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm3<-c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    g <- rnorm(nrow(bin_data_anonymised), -4, 2)
    gamma_mu <- rnorm(1, 0, 2)   
    #tau <- runif(1, 0, 1.5)  
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    #tau<-abs(rcauchy(1, location = 0, scale = 0.5))
    D_star <- run_model_logor(g, gamma_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    w1<-1/(p1_hat*(1-p1_hat)/N + tau^2)
    w1<-w1/sum(w1)
    w3<-1/(se_obs^2 + tau^2)
    w3<-w3/sum(w3)
    dist_arm1 <- distance_abs(sum(w1*p1_hat), sum(w1*p1_star))
    dist_arm3<- distance_abs(sum(w3*D_star$logOR),sum(w3*lor_data_anonymised$yi))
    
    #dist_arm1 <- distance_abs(sum(bin_data_anonymised$n1)/sum(N1), sum(D_star$n1)/sum(N1))
    #dist_arm2 <- distance_abs(sum(bin_data_anonymised$n2)/sum(N2), sum(D_star$n1)/sum(N1))
    #dist_arm3<- weighted_distance_abs(D_star$logOR, lor_data_anonymised$yi,tau,lor_data_anonymised$vi)
    
    
    distance<-dist_arm1+dist_arm3
    
    if (distance <= epsilon) {
      accepted$g <- rbind(accepted$g, g) 
      accepted$gamma_mu <- c(accepted$gamma_mu, gamma_mu)
      accepted$tau <- c(accepted$tau, tau)
      accepted$distance <- c(accepted$distance, distance)
      accepted$dist_arm1 <- c(accepted$dist_arm1, dist_arm1)
      accepted$dist_arm3 <- c(accepted$dist_arm3, dist_arm3)
      i <- i + 1
    }
    
    j <- j + 1
    acc_rate <- i / j  # acceptance rate
    cat("current acceptance rate = ", round(acc_rate * 100, 5), "%\r")
    
  }
  
  return(accepted)
}
################# checking ################   

library(ggplot2)
set.seed(11)
outputs_logor_1<-run_abc_logor_1(N)
outputs_logor_2<-run_abc_logor_2(N)
outputs_logor_3<-run_abc_logor_3(N)
df_means <- data.frame(
  Method = c("ABC1", "ABC2", "ABC3"),
  dist_arm1_mean = c(mean(outputs_logor_1$dist_arm1),
                     mean(outputs_logor_2$dist_arm1),
                     mean(outputs_logor_3$dist_arm1)),
  dist_arm3_mean = c(mean(outputs_logor_1$dist_arm3),
                     mean(outputs_logor_2$dist_arm3),
                     mean(outputs_logor_3$dist_arm3))
)

print(df_means)





########################################## MCMC ############################################



library(nimble)
library(coda)
# N: total number of binary studies
# M: total number of logOR studies

jointModelCode <- nimbleCode({
  # Hyperpriors
  gamma_mu ~ dnorm(0, 1/16)           # Prior for population mean of logOR
  tau ~ T(dnorm(0, 1 / (0.5^2)), 0, )
  
  #tau ~ dunif(0, 1.5)
  
  # Binary data likelihood
  for (i in 1:N) {
    g[i] ~ dnorm(-4,1/16)
    gamma[i] ~ dnorm(gamma_mu, 1/(tau^2))       # Study-specific logOR deviation
    logit(p2[i]) <- g[i]                # Control group log-odds
    logit(p1[i]) <- g[i] + gamma[i]  # Treatment group log-odds
    
    n2[i] ~ dbin(p2[i], N2[i])          # Observed control group events
    n1[i] ~ dbin(p1[i], N1[i])          # Observed treatment group events
  }
  
  # logOR data likelihood
  for (j in 1:M) {
    gamma[N + j] ~ dnorm(gamma_mu, 1/(tau^2))   # Additional delta for logOR studies
    logOR[j] ~ dnorm(gamma[N+j], 1 / (se_obs[j]^2))
  }
})

# Required data: N1, N2, n1, n2 from binary studies
#                logOR, se_obs from logOR studies
#                g: baseline log-odds for control groups in binary studies

# Constants
jointModelConsts <- list(
  N = length(N1),
  M = nrow(lor_data_anonymised),
  N1 = N1,
  N2 = N2,
  se_obs = se_obs
)

# Data
jointModelData <- list(
  n1 = bin_data_anonymised$n1,
  n2 = bin_data_anonymised$n2,
  logOR = lor_data_anonymised$yi
)

# Initial values
jointModelInits <- list(
  gamma_mu = 0,
  tau = 0.1
)

# Parameters to monitor
jointModelParams <- c("gamma_mu","tau","gamma","g")
set.seed(14)
# Run MCMC
samples <- nimbleMCMC(
  code = jointModelCode,
  constants = jointModelConsts,
  data = jointModelData,
  inits = jointModelInits,
  monitors = jointModelParams,
  nchains = 2,
  niter = 100000,
  nburnin = 8000,
  thin = 5
)
coda_samples1 <- mcmc(as.matrix(samples$chain1))
coda_samples2 <- mcmc(as.matrix(samples$chain2))
mcmc_samples <- mcmc.list(coda_samples1, coda_samples2)
trace_samples <- window(mcmc_samples)
plot(trace_samples[, c("gamma_mu","tau","gamma[18]")])
summary(mcmc_samples[, c("gamma_mu","tau","gamma[18]")])


gelman.diag(mcmc_samples)

gelman.plot(mcmc_samples[, c("gamma_mu","tau")])




mcmc_df <- as.data.frame(do.call(rbind, mcmc_samples))  


#Compute the densities of MCMC and ABC
dens_mcmc <- density(mcmc_df$gamma_mu)
dens_abc_1  <- density(outputs_logor_1$gamma_mu)
dens_abc_2<-density(outputs_logor_2$gamma_mu)
dens_abc_3<-density(outputs_logor_3$gamma_mu)

plot(dens_mcmc, col = "blue", lwd = 2, 
     ylim = c(0, 5),
     xlab = expression(theta), ylab = "Density")


lines(dens_abc_1, col = "red", lwd = 2)
lines(dens_abc_2, col = "green", lwd = 2)
lines(dens_abc_3, col = "orange", lwd = 2)

legend("topright", legend = c("MCMC", "ABC1","ABC2","ABC3"), 
       col = c("blue", "red","green","orange"), lwd = 2)
mean(outputs_logor_1$gamma_mu)
mean(outputs_logor_2$gamma_mu)
mean(outputs_logor_3$gamma_mu)
ci_lor_1 <- quantile(outputs_logor_1$gamma_mu, probs = c(0.025, 0.975))
ci_lor_2 <- quantile(outputs_logor_2$gamma_mu, probs = c(0.025, 0.975))
ci_lor_3 <- quantile(outputs_logor_3$gamma_mu, probs = c(0.025, 0.975))
print(ci_lor_1)
print(ci_lor_2)
print(ci_lor_3)

dens_mcmc <- density(mcmc_df$tau)
dens_abc_1  <- density(outputs_logor_1$tau)
dens_abc_2<-density(outputs_logor_2$tau)
dens_abc_3<-density(outputs_logor_3$tau)

plot(dens_mcmc, col = "blue", lwd = 2, 
     ylim = c(0, 5),
     main="The posterior distribution of tau obtained by using MCMC and ABCs",
     xlab = expression(tau), ylab = "Density")

lines(dens_abc_1, col = "red", lwd = 2)
lines(dens_abc_2, col = "green", lwd = 2)
lines(dens_abc_3, col = "orange", lwd = 2)

legend("topright", legend = c("MCMC", "ABC1","ABC2","ABC3"), 
       col = c("blue", "red","green","orange"), lwd = 2)






g_cols <- grep("^g\\[", names(mcmc_df), value = TRUE)

n_post <- 2000  

calc_discrepancy <- function(y) {

  c(mean = mean(y), var = var(y))
  
}

generate_logOR <- function(gamma_mu, tau, lor_data, se_obs) {
  delta <- rnorm(nrow(lor_data), 0, tau)
  gamma_i <- gamma_mu + delta
  logOR_sim <- rnorm(nrow(lor_data), gamma_i, se_obs)
  return(logOR_sim)
}


generate_binary <- function(g, mu, tau) {
  
  gamma<-rnorm(length(N1), mean = mu, sd = tau)
  logit_p2 <- g
  logit_p1 <- g+gamma
  p1 <- 1 / (1 + exp(-logit_p1))
  p2 <- 1 / (1 + exp(-logit_p2))
  n1_star <- rbinom(length(N1), N1, p1)
  n2_star <- rbinom(length(N2), N2, p2)
  return(list(p1_sim = n1_star / N1, p2_sim = n2_star / N2))
}
set.seed(14)
D1_rep <- replicate(n_post, {
  
  i <- sample(1:nrow(mcmc_df), 1)
  g_post <- mcmc_df[i, g_cols]
  g_post<-as.numeric(g_post)
  
  gamma_post<-mcmc_df$gamma_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim<-generate_binary(g_post,gamma_post, tau_post)
  calc_discrepancy(sim$p1_sim)
})



D2_rep <- replicate(n_post, {
 
  i <- sample(1:nrow(mcmc_df), 1)
  g_post <- mcmc_df[i, g_cols]
  g_post<-as.numeric(g_post)
  
  gamma_post<-mcmc_df$gamma_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim<-generate_binary(g_post,gamma_post, tau_post)
  calc_discrepancy(sim$p2_sim)
})



Dlor_rep <- replicate(n_post, {
  i <- sample(1:nrow(mcmc_df), 1)
  gamma_post<-mcmc_df$gamma_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim_lor<-generate_logOR(gamma_post, tau_post, lor_data_anonymised,se_obs)
  calc_discrepancy(sim_lor)
  
})
pval_mean_Dlor <- mean(Dlor_rep["mean", ] >= mean(obs_lor))
pval_var_Dlor <- mean(Dlor_rep["var", ] >= var(obs_lor))
pval_mean_D2 <- mean(D2_rep["mean", ] >=mean(p2_hat) )
pval_var_D2 <- mean(D2_rep["var",] >= var(p2_hat))
pval_mean_D1 <- mean(D1_rep["mean", ] >= mean(p1_hat))
pval_var_D1 <- mean(D1_rep["var", ] >= var(p1_hat))

print(pval_mean_D1)
print(pval_var_D1)
print(pval_mean_D2)
print(pval_var_D2)
print(pval_mean_Dlor)
print(pval_var_Dlor)
obs_lor<-lor_data_anonymised$yi
###########################################
D1_rep <- replicate(n_post, {
  i <- sample(1:nrow(outputs_logor_3$g), 1)
  g_post<-outputs_logor_3$g[i,]
  gamma_post<-outputs_logor_3$gamma_mu[i]
  tau_post<-outputs_logor_3$tau[i]
  sim<-generate_binary(g_post,gamma_post, tau_post)
  calc_discrepancy(sim$p1_sim)
})

D2_rep <- replicate(n_post, {
  
  i <- sample(1:nrow(outputs_logor_3$g), 1)
  g_post<-outputs_logor_3$g[i,]
  gamma_post<-outputs_logor_3$gamma_mu[i]
  tau_post<-outputs_logor_3$tau[i]
  sim<-generate_binary(g_post,gamma_post, tau_post)
  calc_discrepancy(sim$p2_sim)
})


Dlor_rep <- replicate(n_post, {
  i <- sample(1:nrow(outputs_logor_3$g), 1)
  gamma_post<-outputs_logor_3$gamma_mu[i]
  tau_post<-outputs_logor_3$tau[i]
  sim_lor<-generate_logOR(gamma_post, tau_post, lor_data_anonymised,se_obs)
  calc_discrepancy(sim_lor)
  
})


sources <- list(
  ABC1 = outputs_logor_1,
  ABC2 = outputs_logor_2,
  ABC3 = outputs_logor_3
)

n_post <- 2000
results <- data.frame()
set.seed(73)
for (nm in names(sources)) {
  obj <- sources[[nm]]
  
  D1_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$g), 1)
    sim <- generate_binary(as.numeric(obj$g[i, ]), obj$gamma_mu[i], obj$tau[i])
    calc_discrepancy(sim$p1_sim)
  })
  
  D2_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$g), 1)
    sim <- generate_binary(as.numeric(obj$g[i, ]), obj$gamma_mu[i], obj$tau[i])
    calc_discrepancy(sim$p2_sim)
  })
  
  Dlor_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$g), 1)
    sim_lor <- generate_logOR(obj$gamma_mu[i], obj$tau[i], lor_data_anonymised, se_obs)
    calc_discrepancy(sim_lor)
  })
  
  obs_lor <- lor_data_anonymised$yi
  
  results <- rbind(
    results,
    data.frame(
      Source  = nm,
      D1_mean = mean(D1_rep["mean", ] >= mean(p1_hat)),
      D1_var  = mean(D1_rep["var",  ] >= var(p1_hat)),
      D2_mean = mean(D2_rep["mean", ] >= mean(p2_hat)),
      D2_var  = mean(D2_rep["var",  ] >= var(p2_hat)),
      LOR_mean = mean(Dlor_rep["mean",] >= mean(obs_lor)),
      LOR_var  = mean(Dlor_rep["var", ] >= var(obs_lor))
    )
  )
}

print(results)





