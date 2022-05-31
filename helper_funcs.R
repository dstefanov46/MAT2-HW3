B <- 0.05

compute_banana_prob <- function(x) {
  exp(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2)
}

minus_logf <- function(x) {
  -(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2 )
}

minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}

minus_logf_mvn <- function(x) {
  log(2*pi) + (x[1]^2 + x[2]^2)/2
}

minus_logf_grad_mvn <- function(x) {
  x
}

minus_logf_log_reg <- function(x) {
  - compute_log_likelihood(X, y, x)
}

minus_logf_grad_log_reg <- function(betas) {
  grad <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    p <- apply(X, 1, function(x, ...) compute_p_2(x, betas))
    grad[j] <- -sum((y - p) * X[,j])
  }
  grad
}

# plot sample density
plot_density <- function(m, x, mean, colour) {
  df <- data.frame(ID = 1:(m+1), x = x)
  g1 <- geom_density(aes(x = x), df, color = colour)
  g1  
}

# do diagnostics
do_diagnostics <- function(x, m, exec_time, seed) {
  # autocorrelation
  calc_acf <- sum(acf(x, lag.max = 100, plot=F, type="correlation")$acf)
  acf(x, lag.max=100)
  cat("ACF: ", calc_acf, "\n")
  
  # traceplot
  df <- data.frame(ID = 1:(m+1), x = x)
  g1 <- geom_line(data = df[1:m,], aes(x = ID, y = x), colour=colours[seed])
  
  # ess, ess/second
  calc_ess <- ess(x)
  calc_ess_s <- ess(x) / exec_time
  cat("ESS: ", calc_ess, "\n")
  cat("ESS/second: ", calc_ess_s, "\n")
  
  list(g1=g1, ess=calc_ess, ess_s=calc_ess_s, acf=calc_acf) 
}

compute_p_2 <- function(x, betas){
  1 / (1 + exp(-sum(betas*x)))
}

# compute first term of log-likelihood sum for only one data point
compute_log_p <- function(x, betas){
  -log(1 + exp(sum(betas*x)))
}

# compute dot product
compute_dot_prod <- function(x, betas) {
  sum(betas*x)
}

# compute log reg likelihood
compute_log_likelihood <- function(X, y, betas){
  term_1 <- apply(X, 1, function(x, ...) compute_log_p(x, betas))
  term_2 <- apply(X, 1, function(x, ...) compute_dot_prod(x, betas))
  sum(term_1 + term_2 * y)
}

save_samples <- function(fname, sample_list) {
  for (i in 1:5) {
    name <- paste0(fname,
                   i, ".csv")
    write.csv(data.frame(sample_list[[i]]), 
              name)
  }
}

save_metrics <- function(plots, fname) {
  write.csv(plots$ess_list, paste0(fname, "ess.csv"))
  write.csv(plots$ess_s_list, paste0(fname, "ess_s.csv"))
  write.csv(plots$acf_list, paste0(fname, "acf.csv"))
  write.csv(plots$betas, paste0(fname, "betas.csv"))
  write.csv(plots$std_devs, paste0(fname, "std_devs.csv"))
}

calc_mean_and_sd <- function(plots) {
  cat("ESS sd: ", apply(plots$ess_list, 2, "sd"), "\n")
  cat("ESS/second sd: ", apply(plots$ess_s_list, 2, "sd"), "\n")
  cat("ACF sd: ", apply(plots$acf_list, 2, "sd"), "\n")
  cat("Betas sd: ", apply(plots$betas, 2, "sd"), "\n")
  cat("Std sd: ", apply(plots$std_devs, 2, "sd"), "\n")
  cat("ESS mean: ", colMeans(plots$ess_list), "\n")
  cat("ESS/second mean: ", colMeans(plots$ess_s_list), "\n")
  cat("ACF mean: ", colMeans(plots$acf_list), "\n")
  cat("Betas mean: ", colMeans(plots$betas), "\n")
  cat("Std mean: ", colMeans(plots$std_devs), "\n")
}

colours <- c("red", "green", "blue", "yellow", "orange")
# colours <- rep("black", 5)