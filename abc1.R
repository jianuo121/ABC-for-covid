# ABC script

########
# functions

#
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
  n1_star <- rbinom(1,
                    size = N1,
                    prob = prob1)
  n2_star <- rbinom(1,
                    size = N2,
                    prob = prob2)

  return(data.frame(n1_star = n1_star,
                    n2_star = n2_star))
}

#
distance_KL <- function(p, q) {
  sum(p * log(p / q) + (1 - p) * log((1 - p) / (1 - q)))
}

#
distance_euclidian <- function(a, b) {
  (a - b)^2
}

#
distance_jensen_shannon <- function(p, q) {
  # Calculate the average distribution M
  m <- (p + q) / 2

  0.5 * distance_KL(p, m) + 0.5 * distance_KL(q, m)
}

# L1 distance
distance_abs <- function(a, b) {
  abs(a - b)
}

# may way to use log equivalent functions?


run_abc <- function(random_walk = FALSE) {
  #########
  # set-up

  # dummy
  bin_data <- data.frame(N1 = 100, # arm sample sizes
                         N2 = 100,
                         n1 = 20,  # successes
                         n2 = 10)

  N <- 5000 # Number of accepted particles

  # Given N1 N2
  N1 <- bin_data$N1
  N2 <- bin_data$N2

  p1_hat <- bin_data$n1 / N1
  p2_hat <- bin_data$n2 / N2

  epsilon <- 0.2  # distance threshold / tolerance

  res <- matrix(NA, nrow = N, ncol = 3)
  colnames(res) <- c("prob1", "prob2", "distance")

  i <- 1  # counter of accepted
  j <- 1  # counter of proposed

  # initialize prob1 and prob2
  # probs_star <- runif(2, 0.05, 0.95)
  probs_star <- c(0.2, 0.1)

  logit_probs <- vector(length = 2, mode = "numeric")
  dist_arm <- vector(length = 2, mode = "numeric")
  proposal_sd <- 0.05
  eps <- 1e-10   # small number
  save <- TRUE

  p1_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
  p2_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)

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

    # nudge away from zero
    probs_star[1] <- pmin(pmax(probs_star[1], eps), 1 - eps)
    probs_star[2] <- pmin(pmax(probs_star[2], eps), 1 - eps)

    # simulate data
    D_star <- run_model(N1, N2, probs_star[1], probs_star[2])

    p1_star <- D_star$n1_star / N1
    p2_star <- D_star$n2_star / N2

    rr_hat <- p2_hat / p1_hat
    rr_star <- p2_star / p1_star

    # dist_arm[1] <- distance_KL(p1_hat, p1_star)
    # dist_arm[2] <- distance_KL(p2_hat, p2_star)

    dist_arm[1] <- distance_euclidian(p1_hat, p1_star)
    dist_arm[2] <- distance_euclidian(p2_hat, p2_star)

    distance <- sqrt(sum(dist_arm))  # euclidean
    # distance <- sum(dist_arm) / 2

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

    print(cat("current acceptance rate = ", acc_rate, "\r"))
  }

  if (save) write.csv(res, "ABC_1_res.csv", row.names = FALSE)

  res
}

res <- run_abc()


#########
# plots

dens_prob1 <- density(res[,1])
dens_prob2 <- density(res[,2])

p1_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)
p2_prior_params <- beta_tranform(mu = 0.5, sd = 0.2)

dens_prior1 <- density(rbeta(100, shape1 = p1_prior_params$alpha, shape2 = p1_prior_params$beta))
dens_prior2 <- density(rbeta(100, shape1 = p2_prior_params$alpha, shape2 = p2_prior_params$beta))

plot(dens_prob1, col = "blue", lwd = 2,
     xlim = c(0,1), ylim = c(0, 12),
     main = "Density of prob1 and prob2",
     xlab = "Probability", ylab = "Density")
lines(dens_prob2, col = "red", lwd = 2)
legend("topleft", legend = c("prob1", "prob2"),
       col = c("blue", "red"), lwd = 2)

lines(dens_prior1, col = "blue", lty = 2, lwd = 2)
lines(dens_prior2, col = "red", lty = 2, lwd = 2)
