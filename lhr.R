# define half-normal sampler
rhalfnorm <- function(n, sigma) {
  abs(rnorm(n, 0, sigma))
}

sample_sizes <- c(2000, 5000, 8000)

# prior sets
prior_sets <- list(
  A = list(tau = function() rhalfnorm(1, 0.5),
           lhr_mu = function() rnorm(1, 0, 4),
           h = function(n) rnorm(n, -4, 4)),
  B = list(tau = function() rhalfnorm(1, 0.25),
           lhr_mu = function() rnorm(1, 0, 2),
           h = function(n) rnorm(n, -4, 2)),
  C = list(tau = function() rhalfnorm(1, 0.25),
           lhr_mu = function() rnorm(1, 0, 1),
           h = function(n) rnorm(n, -4, 1))
)


results_df <- data.frame(
  sample_size = integer(),
  prior_set = character(),
  epsilon = numeric(),
  mean_distance = numeric(),
  sd_distance = numeric(),
  stringsAsFactors = FALSE
)
set.seed(32)
for (n_presample in sample_sizes) {
  for (prior_name in names(prior_sets)) {
    
    cat("Running presample =", n_presample, "with prior =", prior_name, "\n")
    
    distances_presample <- numeric(n_presample)
    prior <- prior_sets[[prior_name]]
    
    for (k in 1:n_presample) {
      # sample parameters
      
      h <- prior$h(nrow(bin_data_anonymised))
      lhr_mu <- prior$lhr_mu()
      tau <- prior$tau()
      
      # simulate
      D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
      
      # calculate p1*, p2*
      p1_star <- D_star$n1_star / N1
      p2_star <- D_star$n2_star / N2
      
      p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
      p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
      
      # distances
      dist_arm1 <- weighted_distance_KL(p1_hat, p1_star, N1, tau)
      dist_arm3 <- weighted_distance_abs(D_star$logHR,
                                         lhr_data_anonymised$yi,
                                         tau,
                                         lhr_data_anonymised$vi)
      
      distances_presample[k] <- dist_arm1 + dist_arm3
    }
    
    epsilon <- quantile(distances_presample, 0.01)
    cat("Selected epsilon for", prior_name, "with", n_presample, ":", epsilon, "\n")
    
    
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
run_model_lhr <- function(h, lhr_mu, tau, N1, N2, se_obs) {
  
  
  # binary data
  lh2 <- h 
  p2 <- 1 - exp(-exp(lh2))
  n2_star <- rbinom(length(N2), N2, p2)
  
  
  delta <- rnorm(length(N1)+nrow(lhr_data_anonymised), 0, tau)
  lh1 <- h  + lhr_mu + delta[1:length(N1)] 
  p1 <- 1 - exp(-exp(lh1)) 
  n1_star <- rbinom(length(N1), N1, p1) # simulation
  # loghr data
  lhr_i <- lhr_mu+delta[-(1:length(N1))]
  logHR <- rnorm(nrow(lhr_data_anonymised), lhr_i, se_obs)
  
  return(list(n1_star = n1_star, n2_star = n2_star, logHR=logHR))
}



################# some fixed values setting ##############
# given N1 N2
N1 <- bin_data_anonymised$N1
N2 <- bin_data_anonymised$N2
se_obs=sqrt(lhr_data_anonymised$vi)
# number of accepted particles
N <- 1000 

# observed p1 and p2
p1_hat <- bin_data_anonymised$n1 / N1
p2_hat <- bin_data_anonymised$n2 / N2




################## ABC Alogrithm ###################

run_abc_lhr_1 <- function(N) {
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    
    # simulation data
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    dist_arm1 <- weighted_distance_KL(p1_hat, p1_star,N1,tau)
    dist_arm3<- weighted_distance_abs(D_star$logHR, lhr_data_anonymised$yi,tau,lhr_data_anonymised$vi)
    distances_presample[k] <- dist_arm1+dist_arm3
    
  }
  # pick epsilon: 1th quantile of distance_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  
  accepted <- list(h = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   lhr_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm3=c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    
    dist_arm1 <- weighted_distance_KL(p1_hat, p1_star,N1,tau)
    dist_arm3<- weighted_distance_abs(D_star$logHR, lhr_data_anonymised$yi,tau,lhr_data_anonymised$vi)
    distance<-dist_arm1+dist_arm3
    
    
    if (distance < epsilon) {
      accepted$h <- rbind(accepted$h, h) 
      accepted$lhr_mu <- c(accepted$lhr_mu, lhr_mu)
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


run_abc_lhr_2 <- function(N) {
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    
    # simulation data
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
    # dsiatance
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    # calculate distance
    dist_arm1 <- weighted_distance_l1(p1_hat, p1_star,tau,p1_star*(1-p1_star)/N1)
    dist_arm3<- weighted_distance_l1(D_star$logHR, lhr_data_anonymised$yi,tau,lhr_data_anonymised$vi)
    
    distances_presample[k] <- dist_arm1+dist_arm3
    
  }
  # pick epsilon: 1th quantile of distance_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  
  accepted <- list(h = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   lhr_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm3=c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    
    dist_arm1 <- weighted_distance_l1(p1_hat, p1_star,tau,p1_star*(1-p1_star)/N1)
    dist_arm3<- weighted_distance_l1(D_star$logHR, lhr_data_anonymised$yi,tau,lhr_data_anonymised$vi)
    
    distance<-dist_arm1+dist_arm3
    
    
    if (distance < epsilon) {
      accepted$h <- rbind(accepted$h, h) 
      accepted$lhr_mu <- c(accepted$lhr_mu, lhr_mu)
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


run_abc_lhr_3 <- function(N) {
  # presample to choose epsilon
  distances_presample <- numeric(5000)
  for (k in 1:5000) {
    
    # Sampling from prior
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    
    # simulation data
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
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
    dist_arm3<- distance_abs(sum(w3*D_star$logHR),sum(w3*lhr_data_anonymised$yi))
    distances_presample[k] <- dist_arm1+dist_arm3
  }
  # pick epsilon: 1th quantile of distance_presample
  epsilon <- quantile(distances_presample, 0.01)
  cat("Selected epsilon: ", epsilon, "\n")
  
  
  accepted <- list(h = matrix(numeric(0), nrow = 0, ncol = nrow(bin_data_anonymised)),
                   lhr_mu = c(), tau = c(), distance = c(), dist_arm1=c(),dist_arm3=c())
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  eps<-1e-10
  
  while(i <= N) {
    
    h <- rnorm(nrow(bin_data_anonymised), -4, 2)
    lhr_mu <- rnorm(1, 0, 2)   
    tau<-abs(rnorm(1, mean = 0, sd = 0.25))
    
    D_star <- run_model_lhr(h, lhr_mu, tau, N1, N2, se_obs)
    
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star/ N2
    
    
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    w1<-1/(p1_hat*(1-p1_hat)/N + tau^2)
    w1<-w1/sum(w1)
    w3<-1/(se_obs^2 + tau^2)
    w3<-w3/sum(w3)
    dist_arm1 <- distance_abs(sum(w1*p1_hat), sum(w1*p1_star))
    dist_arm3<- distance_abs(sum(w3*D_star$logHR),sum(w3*lhr_data_anonymised$yi))
    distance<-dist_arm1+dist_arm3
    
    
    if (distance < epsilon) {
      accepted$h <- rbind(accepted$h, h) 
      accepted$lhr_mu <- c(accepted$lhr_mu, lhr_mu)
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
#################### results ###################   


set.seed(90)
outputs_lhr_1<-run_abc_lhr_1(N)
outputs_lhr_2<-run_abc_lhr_2(N)
outputs_lhr_3<-run_abc_lhr_3(N)
df_means <- data.frame(
  Method = c("ABC1", "ABC2", "ABC3"),
  dist_arm1_mean = c(mean(outputs_lhr_1$dist_arm1),
                     mean(outputs_lhr_2$dist_arm1),
                     mean(outputs_lhr_3$dist_arm1)),
  dist_arm3_mean = c(mean(outputs_lhr_1$dist_arm3),
                     mean(outputs_lhr_2$dist_arm3),
                     mean(outputs_lhr_3$dist_arm3))
)

print(df_means)



########################################## MCMC ############################################


library(nimble)
library(coda)
# N: total number of binary studies
# M: total number of logHR studies

jointModelCode <- nimbleCode({
  # Hyperpriors
  lhr_mu ~ dnorm(0, 1/16 )          # Prior for population mean of logHR
  tau ~ T(dnorm(0, 1 / (0.5^2)), 0, )
  
  # Binary data likelihood
  for (i in 1:N) {
    h[i] ~ dnorm(-4,1/16)
    
    
    delta[i] ~ dnorm(0, 1/(tau^2))      # Study-specific logHR deviation
    p2[i] <- 1 - exp(-exp(h[i]))
    p1[i] <- 1 - exp(-exp(h[i] + lhr_mu + delta[i]))
                   # Control group log-odds
      # Treatment group log-odds
    
    n2[i] ~ dbin(p2[i], N2[i])          # Observed control group events
    n1[i] ~ dbin(p1[i], N1[i])          # Observed treatment group events
  }
  
  # logHR data likelihood
  for (j in 1:M) {
    delta[N + j] ~ dnorm(0, 1/(tau^2))   # Additional delta for logHR studies
    logHR[j] ~ dnorm(lhr_mu + delta[N + j], 1 / (se_obs[j]^2))
  }
})

               

# Constants
jointModelConsts <- list(
  N = length(N1),
  M = nrow(lhr_data_anonymised),
  N1 = N1,
  N2 = N2,
  se_obs = se_obs
)

# Data
jointModelData <- list(
  n1 = bin_data_anonymised$n1,
  n2 = bin_data_anonymised$n2,
  logHR = lhr_data_anonymised$yi
)

# Initial values
jointModelInits <- list(
  lhr_mu = 0,
  tau = 0.1,
  delta = rep(0, length(N1) + nrow(lhr_data_anonymised)),
  h = rep(0, length(N1))
)

# Parameters to monitor
jointModelParams <- c("lhr_mu", "tau", "h")
set.seed(78)
# Run MCMC
samples <- nimbleMCMC(
  code = jointModelCode,
  constants = jointModelConsts,
  data = jointModelData,
  inits = jointModelInits,
  monitors = jointModelParams,
  nchains = 2,
  niter = 300000,
  nburnin = 50000,
  thin = 10
)
coda_samples1 <- mcmc(as.matrix(samples$chain1))
coda_samples2 <- mcmc(as.matrix(samples$chain2))
mcmc_samples <- mcmc.list(coda_samples1, coda_samples2)
trace_samples <- window(mcmc_samples)
plot(trace_samples[, c("lhr_mu", "tau")])

summary(mcmc_samples[, c("lhr_mu", "tau")])
ci_lhr <- quantile(outputs_lhr$lhr_mu, probs = c(0.025, 0.975))
print(ci_lhr)
ci_tau <- quantile(outputs_lhr$tau, probs = c(0.025, 0.975))
print(ci_tau)
mean(outputs_lhr$lhr_mu)
sd(outputs_lhr$lhr_mu)


mcmc_df <- as.data.frame(do.call(rbind, mcmc_samples))  

############################PPC#########################
# Extracting all h[i]
h_cols <- grep("^h\\[", names(mcmc_df), value = TRUE)

generate_logHR <- function(lhr_mu, tau, lhr_data, se_obs) {
  delta <- rnorm(nrow(lhr_data), 0, tau)
  lhr_i <- lhr_mu + delta
  logOR_sim <- rnorm(nrow(lhr_data), lhr_i, se_obs)
}


calc_discrepancy <- function(y) {
  c(mean = mean(y), var = var(y))
}
generate_binary_data <- function(h, mu, tau) {
  lh2 <- h
  lhr <- rnorm(length(N1), mean = mu, sd = tau)
  lh1 <- h  + lhr
  p1 <- 1 - exp(-exp(lh1)) 
  p2 <- 1 - exp(-exp(lh2))
  n1 <- rbinom(length(N1), N1,p1)
  n2<-rbinom(length(N2),N2,p2)
  return(list(p1_sim =n1/N1, p2_sim = n2/N2))
}
set.seed(66)
D1_rep <- replicate(n_post, {
  i <- sample(1:nrow(mcmc_df), 1)
  h_post <- mcmc_df[i, h_cols]
  h_post<-as.numeric(h_post)
  lhr_mu_post<-mcmc_df$lhr_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim<-generate_binary_data(h_post,lhr_mu_post, tau_post)
  calc_discrepancy(sim$p1_sim)

})
D2_rep <- replicate(n_post, {
  i <- sample(1:nrow(mcmc_df), 1)
  h_post <- mcmc_df[i, h_cols]
  h_post<-as.numeric(h_post)
  lhr_mu_post<-mcmc_df$lhr_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim<-generate_binary_data(h_post,lhr_mu_post, tau_post)
  calc_discrepancy(sim$p2_sim)
  
})
Dlhr_rep <- replicate(n_post, {
  i <- sample(1:nrow(mcmc_df), 1)
  lhr_mu_post<-mcmc_df$lhr_mu[i]
  tau_post<-mcmc_df$tau[i]
  sim_lhr<-generate_logHR(lhr_mu_post, tau_post, lhr_data_anonymised,se_obs)
  calc_discrepancy(sim_lhr)
  
})
obs_lhr<-lhr_data_anonymised$yi
pval_mean_Dlhr <- mean(Dlhr_rep["mean", ] >= mean(obs_lhr))
pval_var_Dlhr <- mean(Dlhr_rep["var", ] >= var(obs_lhr))
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

sources <- list(
  ABC1 = outputs_lhr_1,
  ABC2 = outputs_lhr_2,
  ABC3 = outputs_lhr_3
)

n_post <- 2000
results <- data.frame()
set.seed(73)
for (nm in names(sources)) {
  obj <- sources[[nm]]
  
  D1_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$h), 1)
    sim <- generate_binary(as.numeric(obj$h[i, ]), obj$lhr_mu[i], obj$tau[i])
    calc_discrepancy(sim$p1_sim)
  })
  
  D2_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$h), 1)
    sim <- generate_binary(as.numeric(obj$h[i, ]), obj$lhr_mu[i], obj$tau[i])
    calc_discrepancy(sim$p2_sim)
  })
  
  Dlor_rep <- replicate(n_post, {
    i <- sample(1:nrow(obj$h), 1)
    sim_lhr <- generate_logHR(obj$lhr_mu[i], obj$tau[i], lhr_data_anonymised, se_obs)
    calc_discrepancy(sim_lhr)
  })
  
  results <- rbind(
    results,
    data.frame(
      Source  = nm,
      D1_mean = mean(D1_rep["mean", ] >= mean(p1_hat)),
      D1_var  = mean(D1_rep["var",  ] >= var(p1_hat)),
      D2_mean = mean(D2_rep["mean", ] >= mean(p2_hat)),
      D2_var  = mean(D2_rep["var",  ] >= var(p2_hat)),
      LHR_mean = mean(Dlor_rep["mean",] >= mean(obs_lhr)),
      LHR_var  = mean(Dlor_rep["var", ] >= var(obs_lhr))
    )
  )
}

print(results)
#########################plots######################

dens_mcmc <- density(mcmc_df$tau)
dens_abc_1  <- density(outputs_lhr_1$tau)
dens_abc_2<-density(outputs_lhr_2$tau)
dens_abc_3<-density(outputs_lhr_3$tau)

plot(dens_mcmc, col = "blue", lwd = 2, 
     ylim = c(0, 5),
     xlab = expression(tau), ylab = "Density")


lines(dens_abc_1, col = "red", lwd = 2)
lines(dens_abc_2, col = "green", lwd = 2)
lines(dens_abc_3, col = "orange", lwd = 2)

legend("topright", legend = c("MCMC", "ABC1","ABC2","ABC3"), 
       col = c("blue", "red","green","orange"), lwd = 2)


dens_mcmc <- density(mcmc_df$lhr_mu)
dens_abc_1  <- density(outputs_lhr_1$lhr_mu)
dens_abc_2<-density(outputs_lhr_2$lhr_mu)
dens_abc_3<-density(outputs_lhr_3$lhr_mu)

plot(dens_mcmc, col = "blue", lwd = 2, 
     ylim = c(0, 5),
     xlim=c(0,1),
     xlab = expression(theta), ylab = "Density")


lines(dens_abc_1, col = "red", lwd = 2)
lines(dens_abc_2, col = "green", lwd = 2)
lines(dens_abc_3, col = "orange", lwd = 2)
例
legend("topright", legend = c("MCMC", "ABC1","ABC2","ABC3"), 
       col = c("blue", "red","green","orange"), lwd = 2)
mean(outputs_lhr_1$lhr_mu)
mean(outputs_lhr_2$lhr_mu)
mean(outputs_lhr_3$lhr_mu)
ci_lor_1 <- quantile(outputs_lhr_1$lhr_mu, probs = c(0.025, 0.975))
ci_lor_2 <- quantile(outputs_lhr_2$lhr_mu, probs = c(0.025, 0.975))
ci_lor_3 <- quantile(outputs_lhr_3$lhr_mu, probs = c(0.025, 0.975))
print(ci_lor_1)
print(ci_lor_2)
print(ci_lor_3)

