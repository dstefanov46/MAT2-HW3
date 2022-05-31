# set path to wd
setwd("C:\\Users\\User\\Desktop\\Fakultet\\2.Letnik\\mat_2\\part_3\\submission")

# import helper functions
source("helper_funcs.R")

# import libs
library(ggplot2)
library(MASS, lib.loc = "C:/Program Files/R/R-4.2.0/library")
library(mvtnorm)
library(mcmcse)
library(gsubfn)

# import dataset
data <- read.csv("datset.csv")
head(data)

# M-H algo implementation
run_MH <- function(x, Sigma, distr, ...) {
  x_new <- mvrnorm(n = 1, mu = x, Sigma = Sigma)
  p <- compute_p(x, x_new, distr, X, y)
  if (runif(1) > p) {
    x
  }
  else {
    x_new
  }
}

# compute acceptance ratio
compute_p <- function(x, x_new, distr, ...){
  if (distr == 1){
    p <- exp(dmvnorm(x_new, mean = rep(0, 2), sigma = diag(2), log = T) -  
             dmvnorm(x,     mean = rep(0, 2), sigma = diag(2), log = T))
  }
  else if (distr == 2) {
    p <- exp(log(compute_banana_prob(x_new)) - log(compute_banana_prob(x)))
  }
  else if (distr == 3 | distr == 4) {
    p <- exp(compute_log_likelihood(X, y, x_new) - 
             compute_log_likelihood(X, y, x))
  }
  pmin(1, p)
}

# perform single experiment
do_sampling <- function(seed, m, distr, Sigma, start, ...){
  set.seed(seed)
  
  x <- matrix(0, nrow = m + 1, ncol = nrow(Sigma))
  
  # parameters
  x[1,] <- start  # starting position
  
  # main
  count <- 0
  for (i in 1:m) {
    x[i + 1,] <- run_MH(x[i,], Sigma, distr, X, y)
    if (all(x[i + 1,]==x[i,])) {
      count <- count + 1
    }
    if (i %% 100 == 0) {
      print(i)
    }
  }
  
  cat("Rejection ratio: ", count / m, "\n")
  x
}

calc_mcse <- function(x) {
  mcse(x)$se
}

# do experiment
do_experiment <- function(n_seeds, x_lim, m, distr, Sigma, start, mean, ...) {
  betas <- rep(0, ncol(Sigma))
  std_devs <- rep(0, ncol(Sigma))
  sample_list <- list()
  ess_list <- rep(0, ncol(Sigma))
  ess_s_list <- rep(0, ncol(Sigma))
  acf_list <- rep(0, ncol(Sigma))
  for (seed in 1:n_seeds) {
    start_time <- Sys.time()
    x <- do_sampling(100 * seed, m, distr, Sigma, start, X, y)
    end_time <- Sys.time()
    exec_time <- (end_time - start_time)[[1]]
    sample_list[[seed]] <- x
    cat("The betas are", colMeans(x), "\n")
    betas <- rbind(betas, colMeans(x))
    std_devs <- rbind(std_devs, apply(x, 2, function(y) calc_mcse(y)))
    metrics <- list(ess=rep(0, ncol(Sigma)),
                    ess_s=rep(0, ncol(Sigma)),
                    acf=rep(0, ncol(Sigma)))
    for (i in 1:ncol(Sigma)) {
      diag_list <- do_diagnostics(x[,i], m, exec_time)
      metrics$ess[i] <- diag_list$ess
      metrics$ess_s[i] <- diag_list$ess_s
      metrics$acf[i] <- diag_list$acf
    }
    chains1 <- chains1 + do_diagnostics(x[,1], m, exec_time, seed)$g1 + labs(y = "X")
    chains2 <- chains2 + do_diagnostics(x[,2], m, exec_time, seed)$g1 + labs(y = "Y")
    
    ess_list <- rbind(ess_list, metrics$ess)
    ess_s_list <- rbind(ess_s_list, metrics$ess_s)
    acf_list <- rbind(acf_list, metrics$acf)
  }
  list(chains1 = chains1, chains2=chains2, ess_list=ess_list[2:nrow(betas),], 
       ess_s_list=ess_s_list[2:nrow(betas),], acf_list=acf_list[2:nrow(betas),],
       betas = betas[2:nrow(betas),], std_devs = std_devs[2:nrow(std_devs),],
       sample_list=sample_list)
}

### SIM ------------------------------------------------------------------------
chains1 <- ggplot()
chains2 <- ggplot()

# bivariate standard normal
plots <- do_experiment(5, 5, 1000, 1, diag(2), c(0,0), c(0,0))
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# banana function
plots <- do_experiment(5, 20, 1000, 2, diag(2), c(0,0), c(0,0))
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# log reg with 2 features
X <- data[,1:2]
y <- data[,12]
plots <- do_experiment(5, 3, 1000, 3, diag(ncol(X)), c(0,0), c(0,0), X, y)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)
# results from MLE for log reg (3)
log_reg_1 <- glm(formula = y ~ X2,
                 data = data, 
                 family = binomial)
# compare true to sampled coeffs
data.frame(true=log_reg_1$coefficients, sampled=colMeans(plots$betas),
           sd=colMeans(plots$std_devs))


# log reg with 11 features
X <- data[,1:11]
y <- data[,12]
plots <- do_experiment(5, 3, 1000, 4, 0.01 * diag(ncol(X)), rep(0, ncol(X)), 
                       rep(0, ncol(X)), X, y)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)
# results from MLE for log reg (4)
log_reg_2 <- glm(formula = y ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11,
                 data = data, 
                 family = binomial)
# compare true to sampled coeffs
data.frame(true=log_reg_2$coefficients, sampled=colMeans(plots$betas),
           sd=colMeans(plots$std_devs))
