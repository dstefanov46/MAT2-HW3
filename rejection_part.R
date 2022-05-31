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

max_prob_1 <- dmvnorm(c(0, 0), mean = rep(0, 2), sigma = diag(2))
max_prob_2 <- compute_banana_prob(c(0, 5))

compute_prob <- function(distr, x, ...) {
  if (distr == 1) {
    dmvnorm(x, mean = rep(0, 2), sigma = diag(2))
  }
  else if (distr == 2) {
    compute_banana_prob(x)
  }
  else if (distr == 3) {
    compute_log_likelihood(X, y, x)
  }
}

compute_envelope <- function(distr, x) {
  if (distr == 1) {
    # dmvnorm(x, mean = rep(0, 2), sigma = diag(2))
    max_prob_1
  }
  else if (distr == 2) {
    max_prob_2
  }
  else if (distr == 3) {
    -1226.071
  }
}

sample_from_g <- function(distr) {
  if (distr == 1) {
    runif(2, -10, 10)
  }
  else if (distr == 2) {
    c(runif(1, -20, 20), runif(1, -2, 12)) # (-2, 8) > (-10, 10), which slightly better than (-5, 5)
  }                                        # mean for y below 4!!
  else if (distr == 3) {
    c(runif(1, -3, 3), runif(1, -3, 3)) # (-2, 8) > (-10, 10), which slightly better than (-5, 5)
  }
}

do_rejection_sampling <- function(seed, m, M, mu, Sigma, distr, ...) {
  set.seed(seed)
  x <- matrix(0, nrow = m+1, ncol = nrow(Sigma))
  
  count <- 0
  for (i in 1:(m+1)) {
    p <- 1
    frac <- 0
    while(p >= frac) {
      x[i,] <- sample_from_g(distr)
      p <- runif(1)
      if (distr != 3) {
        frac <- compute_prob(distr, x[i,], X, y) / (M * compute_envelope(distr, x[i,]))
      }
      else {
        frac <- exp(compute_prob(distr, x[i,], X, y) - compute_envelope(distr, x[i,])) / M
      }
      count <- count + 1
    }
    if (i %% 100 == 0) {
      print(i)
    }
  }
  list(x=x, count=count)
}

# do experiment
do_experiment <- function(n_seeds, x_lim, m, M, distr, mu, 
                          Sigma, ...) {
  betas <- rep(0, ncol(Sigma))
  std_devs <- rep(0, ncol(Sigma))
  sample_list <- list()
  ess_list <- rep(0, ncol(Sigma))
  ess_s_list <- rep(0, ncol(Sigma))
  acf_list <- rep(0, ncol(Sigma))
  for (seed in 1:n_seeds) {
    start_time <- Sys.time()
    samples <- do_rejection_sampling(100 * seed, m, M, mu, Sigma, distr, X, y)
    cat("Rejections: ", samples$count)
    end_time <- Sys.time()
    exec_time <- (end_time - start_time)[[1]]
    x <- samples$x
    sample_list[[seed]] <- x
    cat("The betas are", colMeans(x), "\n")
    betas <- rbind(betas, colMeans(x))
    std_devs <- rbind(std_devs, apply(x, 2, "sd"))
    metrics <- list(ess=rep(0, ncol(Sigma)),
                    ess_s=rep(0, ncol(Sigma)),
                    acf=rep(0, ncol(Sigma)))
    for (i in 1:ncol(Sigma)) {
      diag_list <- do_diagnostics(x[,i], m, exec_time)
      metrics$ess[i] <- diag_list$ess
      metrics$ess_s[i] <- diag_list$ess_s
      metrics$acf[i] <- diag_list$acf
    }
    ess_list <- rbind(ess_list, metrics$ess)
    ess_s_list <- rbind(ess_s_list, metrics$ess_s)
    acf_list <- rbind(acf_list, metrics$acf)
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
# SIM --------------------------------------------------------------------------
# initial plot settings
chains1 <- ggplot()
chains2 <- ggplot()

# bivariate standard normal
plots <- do_experiment(5, 5, 999, 1, 1, rep(0, 2), diag(2))  # 999, because of how plot function defined
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# banana function
plots <- do_experiment(5, 20, 999, 1, 2, rep(0, 2), diag(2))  # 999, because of how plot function defined
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# log reg with 2 features
X <- data[,1:2]
y <- data[,12]
plots <- do_experiment(5, 3, 999, 1, 3, rep(0, ncol(X)), diag(ncol(X)), X, y)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

temp_samples <- rep(0, 10000)
for (i in 1:10000) {
  temp_sample <- sample_from_g(3)
  temp_samples[i] <- compute_prob(3, temp_sample, X, y)
}
mean(temp_samples)
