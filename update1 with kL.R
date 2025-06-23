beta_tranform <- function(mu, sd) {
  variance <- sd^2
  
  phi <- (mu * (1 - mu) / variance) - 1
  
  # Calculate alpha and beta
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  
  return(list(alpha = alpha, beta = beta))
}

#
run_model <- function(N1, N2, prob1, prob2) {
  n1_star <- rbinom(length(N1),
                    size = N1,
                    prob = prob1)
  n2_star <- rbinom(length(N2),
                    size = N2,
                    prob = prob2)
  
  return(data.frame(n1_star = n1_star,
                    n2_star = n2_star))
}


distance_KL <- function(p, q) {
sum(p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)))
}

p1_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
p2_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
alpha1<-p1_prior_params$alpha
beta1<-p1_prior_params$beta
alpha2<-p2_prior_params$alpha
beta2<-p2_prior_params$beta
# Given N1 N2
N1 <- bin_data_anonymised$N1
N2 <- bin_data_anonymised$N2
N <- 1000# Number of accepted particles

p1_hat <- bin_data_anonymised$n1 / N1
p2_hat <- bin_data_anonymised$n2 / N2
# presample
distances_presample <- numeric(100)
for (k in 1:100) {
  p1 <- rbeta(1, alpha1, beta1)
  p2 <- rbeta(1, alpha2, beta2)
  D_star <- run_model(N1, N2, p1, p2)
  p1s <- D_star$n1_star / N1
  p2s <- D_star$n2_star / N2
  p1s <- pmin(pmax(p1s, eps), 1 - eps)
  p2s <- pmin(pmax(p2s, eps), 1 - eps)
  d1 <- distance_KL(p1_hat, p1s)
  d2 <- distance_KL(p2_hat, p2s)
  
  distances_presample[k] <- (d1 + d2) / 2
}
epsilon <- quantile(distances_presample, 0.1)
cat("Selected epsilon: ", epsilon, "\n")

run_abc <- function(random_walk = FALSE) {
  
  #epsilon <- 0.4  # distance threshold / tolerance
  
  res <- matrix(NA, nrow = N, ncol = 3)
  colnames(res) <- c("prob1", "prob2", "distance")
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  
  # initialize prob1 and prob2
  # probs_star <- runif(2, 0.05, 0.95)
  probs_star <- c(0.5, 0.7)
  
  logit_probs <- vector(length = 2, mode = "numeric")
  dist_arm <- vector(length = 2, mode = "numeric")
  proposal_sd <- 0.05
  eps <- 1e-10   # small number
  save <- TRUE
  
  
  #### ABC algorithm ####
  
  while(i <= N) {
    
    if (random_walk) {
      # logit transformation
      logit_probs[1] <- log(probs_star[1] / (1 - probs_star[1]))
      logit_probs[2] <- log(probs_star[2] / (1 - probs_star[2]))
      
      logit_proposal <- rnorm(2, mean = logit_probs, sd = proposal_sd)
      
      # transform back to probability
      probs_star <- exp(logit_proposal) / (1 + exp(logit_proposal))
    } else {
      # pure rejection sampling
      probs_star[1] <- rbeta(1, shape1 = p1_prior_params$alpha, shape2 = p1_prior_params$beta)
      probs_star[2] <- rbeta(1, shape1 = p2_prior_params$alpha, shape2 = p2_prior_params$beta)
    }
    
    
    # simulate data
    D_star <- run_model(N1, N2, probs_star[1], probs_star[2])
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star / N2
    
    # nudge away from zero
    p1_star <- pmin(pmax(p1_star, eps), 1 - eps)
    p2_star <- pmin(pmax(p2_star, eps), 1 - eps)
    
    
    dist_arm[1] <- distance_KL(p1_hat, p1_star)
    dist_arm[2] <- distance_KL(p2_hat, p2_star)
  
    distance <- sum(dist_arm) / 2
    
    # Aim for an acceptance rate between 10-40%
    if (distance <= epsilon) {
      # if distance is less than the tolerance
      # Store results
      res[i, ] <- c(probs_star[1], probs_star[2], distance)
      
      i <- i + 1
    }
    
    j <- j + 1
    acc_rate <- i / j  # acceptance rate
    
    # if (acc_rate < 0.2) browser()
    
    cat("current acceptance rate = ", acc_rate, "\r")
  }
  
  if (save) write.csv(res, "ABC_1_res.csv", row.names = FALSE)
  return(list(
    results = res,
    acceptance_rate = acc_rate
  ))
}

dens_prior1 <- density(rbeta(1000, shape1 = alpha1, shape2 = beta1))
dens_prior2 <- density(rbeta(1000, shape1 = alpha2, shape2 = beta2))

colors_prob1 <- adjustcolor("blue", alpha.f = 0.3)
colors_prob2 <- adjustcolor("red", alpha.f = 0.3)
plot(NULL, xlim = c(0, 1),ylim = c(0, 10),
     xlab = "Probability", ylab = "Density", 
     main = "Posterior Densities from 10 Simulations with KL distance")

for (i in 1:10) {
  outs <- run_abc(random_walk = FALSE)
  dens1 <- density(outs$results[,1], from = 0, to = 1)
  dens2 <- density(outs$results[,2], from = 0, to = 1)
  lines(dens1, col = colors_prob1, lwd = 2)
  lines(dens2, col = colors_prob2, lwd = 2)
}

lines(dens_prior1,col="blue",lwd=2,lty=2)

lines(dens_prior2, col = "red", lwd = 2, lty=2)
data_dens1<-density(bin_data_anonymised$n1/N1, from=0, to=1)
lines(data_dens1,col="blue",lty=3)
data_dens2<-density(bin_data_anonymised$n2/N2, from=0, to=1)
lines(data_dens2,col="red",lty=3)
legend("topright", legend = c("Posterior prob1", "Posterior prob2", 
                              "Prior prob1", "Prior prob2","approximately true density of prob1",
                              "approximately true densith of prob2"),
       col = c("blue", "red", "blue", "red","blue","red"),
       lty = c(1, 1, 2, 2,3,3), lwd = 2)

#################################
#RR
distance_KL <- function(p, q) {
sum(p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)))
}

p1_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
p2_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
alpha1<-p1_prior_params$alpha
beta1<-p1_prior_params$beta
alpha2<-p2_prior_params$alpha
beta2<-p2_prior_params$beta
# Given N1 N2
N1 <- bin_data_anonymised$N1
N2 <- bin_data_anonymised$N2
N <- 1000# Number of accepted particles

p1_hat <- bin_data_anonymised$n1 / N1
p2_hat <- bin_data_anonymised$n2 / N2
rr_hat<-p1_hat/p2_hat
eps <- 1e-10   # small number
rr_hat <- pmin(pmax(rr_hat, eps), 1 - eps)
distances_presample <- numeric(100)
for (k in 1:100) {
  p1 <- rbeta(1, alpha1, beta1)
  p2 <- rbeta(1, alpha2, beta2)
  D_star <- run_model(N1, N2, p1, p2)
  p1s <- D_star$n1_star / N1
  p2s <- D_star$n2_star / N2
  rr<-p1s/p2s
  rr <- pmin(pmax(rr, eps), 1 - eps)
  distances_presample[k] <- distance_KL(rr_hat,rr)
}
epsilon <- quantile(distances_presample, 0.1)
cat("Selected epsilon: ", epsilon, "\n")

run_abc <- function(random_walk = FALSE) {
  
  #epsilon <- 0.4  # distance threshold / tolerance
  
  res <- matrix(NA, nrow = N, ncol = 3)
  colnames(res) <- c("prob1", "prob2", "distance")
  
  i <- 1  # counter of accepted
  j <- 1  # counter of proposed
  
  # initialize prob1 and prob2
  # probs_star <- runif(2, 0.05, 0.95)
  probs_star <- c(0.5, 0.7)
  
  logit_probs <- vector(length = 2, mode = "numeric")
  dist_arm <- vector(length = 2, mode = "numeric")
  proposal_sd <- 0.05
  
  save <- TRUE
  
  #### ABC algorithm ####
  
  while(i <= N) {
    
    if (random_walk) {
      # logit transformation
      logit_probs[1] <- log(probs_star[1] / (1 - probs_star[1]))
      logit_probs[2] <- log(probs_star[2] / (1 - probs_star[2]))
      
      logit_proposal <- rnorm(2, mean = logit_probs, sd = proposal_sd)
      
      # transform back to probability
      probs_star <- exp(logit_proposal) / (1 + exp(logit_proposal))
    } else {
      # pure rejection sampling
      probs_star[1] <- rbeta(1, shape1 = p1_prior_params$alpha, shape2 = p1_prior_params$beta)
      probs_star[2] <- rbeta(1, shape1 = p2_prior_params$alpha, shape2 = p2_prior_params$beta)
    }
    
    
    # simulate data
    D_star <- run_model(N1, N2, probs_star[1], probs_star[2])
    
    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star / N2
    
    rr_star <- p2_star / p1_star

    rr_star<-pmin(pmax(rr_star, eps), 1 - eps)
    distance<-distance_KL(rr_hat, rr_star)
    
    # Aim for an acceptance rate between 10-40%
    if (distance <= epsilon) {
      # if distance is less than the tolerance
      # Store results
      res[i, ] <- c(probs_star[1], probs_star[2], distance)
      
      i <- i + 1
    }
    
    j <- j + 1
    acc_rate <- i / j  # acceptance rate
    
    # if (acc_rate < 0.2) browser()
    
    cat("current acceptance rate = ", acc_rate, "\r")
  }
  
  if (save) write.csv(res, "ABC_1_res.csv", row.names = FALSE)
  return(list(
    results = res,
    acceptance_rate = acc_rate
  ))
}

colors_prob_rr <- adjustcolor("green", alpha.f = 0.3)
rr_vec <- (bin_data_anonymised$n1 /N1) /(bin_data_anonymised$n2 /N2)
plot(NULL, xlim = c(0, max(rr_vec)),ylim = c(0, 10),
     xlab = "RR", ylab = "Density", 
     main = "Posterior Densities from 10 Simulations with KL distance")

#Remove outliers using the IQR
q1 <- quantile(rr_vec, 0.25)
q3 <- quantile(rr_vec, 0.75)
iqr <- q3 - q1
rr_vec <- rr_vec[rr_vec > (q1 - 1.5 * iqr) & rr_vec < (q3 + 1.5 * iqr)]

for (i in 1:10) {
  outs <- run_abc(random_walk = FALSE)
  dens_rr <- density(outs$results[,1], from = 0, to = 1)
  lines(dens_rr, col = colors_prob_rr, lwd = 2)
}

data_dens_rr<-density(rr_vec)
lines(data_dens_rr, lty=2, col="red")
legend("topright", legend = c("Posterior of rr",
                             "approximately true density of rr"
                             ),
       col = c("green", "red"),
       lty = c(1,3), lwd = 2)