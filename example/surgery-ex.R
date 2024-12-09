library(xactonomial)
library(here)

psi_risk_1a <- function(theta, lower = TRUE) {

  p01_0 <- theta[1]
  p11_0 <- theta[2]
  p01_1 <- theta[4]
  p11_1 <- theta[5]

  if(lower) {
  max(p01_0 + p11_1 - 1,
      -p11_0 + 2 * p11_1 - 1,
      2 * p01_0 - p01_1 - 1)
  } else {

  min(-p01_1 - p11_0 + 1,
      -2 * p11_0 + p11_1 + 1,
      p01_0 - 2 * p01_1 + 1)
  }

}

psi_risk_1b <- function(theta, lower = TRUE) {

  p01_0 <- theta[1]
  p11_0 <- theta[2]
  p01_1 <- theta[4]
  p11_1 <- theta[5]

  if(lower) {
    p01_0 + p11_1 - 1
  } else {
    1 - p01_1 - p11_0
  }



}


data <- list(T0 = c(81, 93 - 81, 8), T1 = c(18, 29 - 18, 4))
data2 <- list(T0 = c(91 - 7, 0, 7), T1 = c(101 - 3, 0, 3))

surgery1a <- list(lower_ci = xactonomial(\(te) psi_risk_1a(te, lower = TRUE), data, alpha = .1,
                                         psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2],
                  upper_ci = xactonomial(\(te) psi_risk_1a(te, lower = FALSE), data, alpha = .1,
                                         psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2])

antibiot1a <- list(lower_ci = xactonomial(\(te) psi_risk_1a(te, lower = TRUE), data2, alpha = .1,
                                         psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2],
                  upper_ci = xactonomial(\(te) psi_risk_1a(te, lower = FALSE), data2, alpha = .1,
                                         psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2])

antibiot1b <- list(lower_ci = xactonomial(\(te) psi_risk_1b(te, lower = TRUE), data2, alpha = .1,
                                          psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2],
                   upper_ci = xactonomial(\(te) psi_risk_1b(te, lower = FALSE), data2, alpha = .1,
                                          psi_limits = c(-1, 1), maxit = 30, chunksize = 200)[1:2])

exres <- list(surgery1a, antibiot1a, antibiot1b)
saveRDS(exres, file = here("example/ex-results.rds"))

exres <- readRDS(here("example/ex-results.rds"))

sapply(exres[3], \(x){

  c(sprintf("bounds estimate: (%.2f to %.2f)", x$lower_ci$estimate, x$upper_ci$estimate),
    sprintf("95%% confidence interval: (%.2f to %.2f)", x$lower_ci$conf.int[1], x$upper_ci$conf.int[2]))
}) |> cat(sep = "\n")


uncondExact2x2(7, 91, 0, 101,
               parmtype = "difference",
               alternative = "less",
               method = "score",
               conf.level = 0.975,
               conf.int = TRUE)



### testing out adaptation of the sample space to fix the missingness
test <- xactonomial(\(te) psi_risk_1b(te, lower = TRUE), data2, alpha = .1,
            psi_limits = c(-1, 1), maxit = 2, chunksize = 20)
p_func <- test$pvalue_function
SSpace <- get("SSpacearr", environment(p_func))
subspace <- SSpace[, 3] <= 7 & SSpace[, 6] <= 3
SSpace <- SSpace[subspace, ]
assign("SSpacearr", SSpace, envir = environment(p_func))
SSpacearr <- SSpace
logC <- get("logC", environment(p_func))
logC <- logC[subspace]
assign("psi_hat", get("psi_hat", environment(p_func))[subspace],
       envir = environment(p_func))
assign("SSpacearr", SSpace, envir = environment(p_func))


flower <- function(x, maxit, chunksize, target, psi_limits, SSpacearr, logC){
  p_func(x, maxit = maxit, chunksize = chunksize,
              lower = TRUE, target = target, psi_limits = psi_limits,
              SSpacearr = SSpacearr, logC = logC) - .1 / 2
}

fupper <- function(x, maxit, chunksize, target, psi_limits, SSpacearr, logC) {
  pvalue_psi0(x, maxit = maxit, chunksize = chunksize,
              lower = FALSE, target = alpha / 2, psi_limits = psi_limits,
              SSpacearr = SSpacearr, logC = logC) - (.1 / 2)
}



lower_limit <- itp_root(flower, -1, 1,
                        fa = -.1 / 2, fb = 1 - .1 / 2, maxiter = 10,
                        maxit = 30, chunksize = 200,
                        target = .1 / 2, psi_limits = c(-1,1),
                        SSpacearr = SSpacearr, logC = logC)
upper_limit <- itp_root(fupper, psi_limits[1], psi_limits[2],
                        fa = 1-alpha / 2, fb = - alpha / 2, maxiter = 10,
                        maxit = maxit, chunksize = chunksize,
                        target = alpha / 2, psi_limits = psi_limits,
                        SSpacearr = SSpacearr, logC = logC)

