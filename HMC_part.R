# set path to wd
setwd("C:\\Users\\User\\Desktop\\Fakultet\\2.Letnik\\mat_2\\part_3\\submission")

# import helper funcs
source("Multiplot.r")
source("helper_funcs.R")
source("HMC.r")

# import libs
library(ggplot2)
library(MASS, lib.loc = "C:/Program Files/R/R-4.2.0/library")
library(mvtnorm)
library(mcmcse)
library(gsubfn)
library(numDeriv)
library(coda)
library(grid)
library(gridExtra)

# import dataset
data <- read.csv("datset.csv")
head(data)

do_sampling <- function(seed, L, epsilon, current_q, m, U, grad_U, ...) {
  set.seed(seed)
  samples <- rep(0, length(current_q))
  count <- 0
  for (i in 1:m) {
    if (i %% 100 == 0) print(i)
    res = HMC(U, grad_U, epsilon, L, current_q)
    samples = rbind(samples, res$next_q)
    current_q = res$next_q
    if (all(res$next_q==samples[i,])) {
      count <- count + 1
    }
    if (i %% 50 == 0) print(m*effectiveSize(samples)/i) 
  }
  
  cat("Rejection ratio: ", count / m, "\n")
  samples
}

calc_mcse <- function(x) {
  mcse(x)$se
}

# do experiment
do_experiment <- function(n_seeds, x_lim, m, L, epsilon, current_q, 
                          U, grad_U, ...) {
  betas <- rep(0, length(current_q))
  std_devs <- rep(0, length(current_q))
  ess_list <- rep(0, length(current_q))
  ess_s_list <- rep(0, length(current_q))
  acf_list <- rep(0, length(current_q))
  sample_list <- list()
  for (seed in 1:n_seeds) {
    start_time <- Sys.time()
    x <- do_sampling(100 * seed, L, epsilon, current_q, m, U, grad_U, X, y)
    end_time <- Sys.time()
    exec_time <- (end_time - start_time)[[1]]
    sample_list[[seed]] <- x
    cat("The betas are", colMeans(x), "\n")
    betas <- rbind(betas, colMeans(x))
    std_devs <- rbind(std_devs, apply(x, 2, function(y) calc_mcse(y)))
    chains1 <- chains1 + do_diagnostics(x[,1], m, exec_time, seed)$g1 + labs(y = "X")
    chains2 <- chains2 + do_diagnostics(x[,2], m, exec_time, seed)$g1 + labs(y = "Y")
    metrics <- list(ess=rep(0, length(current_q)),
                    ess_s=rep(0, length(current_q)),
                    acf=rep(0, length(current_q)))
    for (i in 1:length(current_q)) {
      diag_list <- do_diagnostics(x[,i], m, exec_time)
      metrics$ess[i] <- diag_list$ess
      metrics$ess_s[i] <- diag_list$ess_s
      metrics$acf[i] <- diag_list$acf
    }
    ess_list <- rbind(ess_list, metrics$ess)
    ess_s_list <- rbind(ess_s_list, metrics$ess_s)
    acf_list <- rbind(acf_list, metrics$acf)
  }
  list(ess_list=ess_list[2:nrow(betas),], chains1=chains1, chains2=chains2,
       ess_s_list=ess_s_list[2:nrow(betas),], acf_list=acf_list[2:nrow(betas),],
       betas = betas[2:nrow(betas),], std_devs = std_devs[2:nrow(std_devs),],
       sample_list=sample_list)
}

### SIM ------------------------------------------------------------------------
chains1 <- ggplot()
chains2 <- ggplot()

# bivariate standard normal
plots <- do_experiment(5, 5, 1000, 27, 0.6, c(0, 0),
                       minus_logf_mvn, minus_logf_grad_mvn)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# banana function
plots <- do_experiment(5, 20, 1000, 27, 0.6, c(0, 0),
                       minus_logf, minus_logf_grad)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)

# log reg with 2 features
X <- data[,1:2]
y <- data[,12]
plots <- do_experiment(5, 3, 1000, 5, 0.01, c(0, 0),
                       minus_logf_log_reg, minus_logf_grad_log_reg)
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
plots <- do_experiment(5, 3, 1000, 5, 0.01, rep(0, ncol(X)),
                       minus_logf_log_reg, minus_logf_grad_log_reg)
grid.arrange(plots$chains1, plots$chains2, ncol=2)
calc_mean_and_sd(plots)
# results from MLE for log reg (4)
log_reg_2 <- glm(formula = y ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11,
                 data = data, 
                 family = binomial)
# compare true to sampled coeffs
data.frame(true=log_reg_2$coefficients, sampled=colMeans(plots$betas),
           sd=colMeans(plots$std_devs))
