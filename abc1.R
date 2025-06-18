run_model<-function(N1,N2,prob1,prob2){
  
  n1_star<- rbinom(nrow(bin_data_anonymised), size = N1, prob =prob1)
  n2_star<- rbinom(nrow(bin_data_anonymised), size = N2, prob =prob2)  
  return(data.frame(n1_star = n1_star, n2_star = n2_star))
}

calc_distance<- function(p, q) {
  
  kl<- sum(p * log(p / q)+(1-p)*log((1-p)/(1-q)))
  return(kl)
}

#### ABC algorithm ####
N <- 1000 # Number of accepted particles
# Given N1 N2
N1<-bin_data_anonymised$N1
N2<-bin_data_anonymised$N2
epsilon<-2
res <- matrix(NA, nrow = N, ncol = 3)
colnames(res) <- c("prob1", "prob2", "distance")
i <- 1 # Initiate counter of accepted particles
j <- 1 # Initiate counter of proposed particles

# initialize prob1 and prob2
probs_star <- runif(2, 0.05, 0.95)

while(i <= N){ # While the number of accepted particles is less than N_particles
  # do logit transformation
  logit_prob1 <- log(probs_star[1] / (1 - probs_star[1]))
  logit_prob2 <- log(probs_star[2] / (1 - probs_star[2]))
  logit_probs<-c(logit_prob1, logit_prob2)
  logit_proposal <- rnorm(2, mean = logit_probs, sd = 0.1)
  #transform back to probability
  probs_star <- exp(logit_proposal) / (1 + exp(logit_proposal))
  eps<- 1e-10
  probs_star[1]<-pmin(pmax(probs_star[1], eps), 1 - eps)
  probs_star[2]<-pmin(pmax(probs_star[2], eps), 1 - eps)
  # Simulate data set from the model
  D_star <- run_model(N1,N2,probs_star[1], probs_star[2])
  
  # Calculate distance  
  distance<- calc_distance(bin_data_anonymised$n1/N1, D_star$n1_star/N1)+ calc_distance(bin_data_anonymised$n2/N2, D_star$n2_star/N2)
  
  if(distance<=epsilon){ # If the distance is less than the tolerance
    # Store results
    res[i,] <- c(probs_star[1], probs_star[2], distance)
    # Update counter
    i <- i + 1
  }
  j <- j + 1 # Update counter
  acc_rate <- i / j # Calculate the acceptance rate 
  cat("current acceptance rate = ", acc_rate, "\r")
}

# Save data to csv file
write.csv(res, "ABC_1_res.csv", row.names = F)

dens_prob1 <- density(res[,1])
dens_prob2 <- density(res[,2])


plot(dens_prob1, col = "blue", lwd = 2, main = "Density of prob1 and prob2",
     xlab = "Probability", ylab = "Density")


lines(dens_prob2, col = "red", lwd = 2)


legend("topright", legend = c("prob1", "prob2"), col = c("blue", "red"), lwd = 2)