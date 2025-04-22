library(xactonomial)


psi_lowerbound <- function(theta, terms = FALSE) {
  p00_0 = theta[,1]
  p00_1 = theta[,5]
  p10_0 = theta[,2]
  p10_1 = theta[,6]
  p01_0 = theta[,3]
  p01_1 = theta[,7]
  p11_0 = theta[,4]
  p11_1 = theta[,8]

  if(terms) {
    rbind(p00_0 - p00_1 - p10_1 - p01_1, p00_0 - p00_1 -
      p10_0 - p10_1 - p01_0, p00_0 - p00_1 + p10_0 - 2 * p10_1 -
      2 * p01_1, -p10_1 - p01_1, -p10_0 - p01_0, -p00_0 + p00_1 -
      2 * p10_0 + p10_1 - 2 * p01_0, -p00_0 + p00_1 - p10_0 -
      p10_1 - p01_1, -p00_0 + p00_1 - p10_0 - p01_0)
  } else {
  pmax(p00_0 - p00_1 - p10_1 - p01_1, p00_0 - p00_1 -
         p10_0 - p10_1 - p01_0, p00_0 - p00_1 + p10_0 - 2 * p10_1 -
         2 * p01_1, -p10_1 - p01_1, -p10_0 - p01_0, -p00_0 + p00_1 -
         2 * p10_0 + p10_1 - 2 * p01_0, -p00_0 + p00_1 - p10_0 -
         p10_1 - p01_1, -p00_0 + p00_1 - p10_0 - p01_0)
  }
}

psi_upperbound <- function(theta, terms = FALSE) {
  p00_0 = theta[,1]
  p00_1 = theta[,5]
  p10_0 = theta[,2]
  p10_1 = theta[,6]
  p01_0 = theta[,3]
  p01_1 = theta[,7]
  p11_0 = theta[,4]
  p11_1 = theta[,8]
  if(terms) {
    rbind(1 - p10_1 - p01_0, 1 + p00_0 + p10_0 - 2 * p10_1 -
             p01_1, 2 - p00_1 - p10_0 - p10_1 - 2 * p01_0, 1 - p10_1 -
             p01_1, 1 - p10_0 - p01_0, 1 + p00_1 - 2 * p10_0 + p10_1 -
             p01_0, 2 - p00_0 - p10_0 - p10_1 - 2 * p01_1, 1 - p10_0 -
             p01_1)
  } else {
  pmin(1 - p10_1 - p01_0, 1 + p00_0 + p10_0 - 2 * p10_1 -
         p01_1, 2 - p00_1 - p10_0 - p10_1 - 2 * p01_0, 1 - p10_1 -
         p01_1, 1 - p10_0 - p01_0, 1 + p00_1 - 2 * p10_0 + p10_1 -
         p01_0, 2 - p00_0 - p10_0 - p10_1 - 2 * p01_1, 1 - p10_0 -
         p01_1)
  }
}

data <- list(Z0 = c(13, 0, 4, 0), Z1 = c(0, 18, 1, 0))
lapply(data, \(x) {

  list(matrix(x, nrow = 2, ncol = 2, byrow = FALSE),
       round(100 * matrix(x / sum(x), nrow = 2, ncol = 2), 1))

})
lapply(data, sum)

xactonomial(data, psi_lowerbound, psi0 = 0, psi_limits = c(-1,1), alternative = "greater",
             maxit = 500, chunksize = 1000)
xactonomial(data, psi_upperbound, psi0 = 0, psi_limits = c(-1,1), alternative = "less",
             maxit = 500, chunksize = 1000)


## bootstrap
round(sapply(list(psi_lowerbound, psi_upperbound),
             \(f) f(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1))), 2)


bsamps <- t(replicate(5000, {
  dstar <- unlist(lapply(data, \(x) {
    newx <- table(factor(sample(rep.int(1:length(x), times = x), sum(x), replace = TRUE),
                         levels = 1:length(x)))

    newx / sum(newx)
  }))
  dstar
}))

psi_lowerbound(bsamps) |> quantile(c(0.025, .975))
psi_upperbound(bsamps) |> quantile(c(0.025, .975))

girlsterms <- cbind(psi_lowerbound(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1), terms = TRUE),
psi_upperbound(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1), terms = TRUE))


######### boys

data <- list(Z0 = c(20, 0, 12, 0), Z1 = c(2, 22, 3, 1))
lapply(data, \(x) {

  list(matrix(x, nrow = 2, ncol = 2, byrow = FALSE),
  round(100 * matrix(x / sum(x), nrow = 2, ncol = 2), 1))

})


xactonomial(data, psi_lowerbound, psi0 = 0, psi_limits = c(-1,1), alternative = "greater",
            maxit = 500, chunksize = 1000)
xactonomial(data, psi_upperbound, psi0 = 0, psi_limits = c(-1,1), alternative = "less",
            maxit = 500, chunksize = 1000)


## bootstrap
round(sapply(list(psi_lowerbound, psi_upperbound),
       \(f) f(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1))), 2)


bsamps <- t(replicate(5000, {
  dstar <- unlist(lapply(data, \(x) {
    newx <- table(factor(sample(rep.int(1:length(x), times = x), sum(x), replace = TRUE),
                         levels = 1:length(x)))

    newx / sum(newx)
  }))
  dstar
}))

psi_lowerbound(bsamps) |> quantile(c(0.025, .975))
psi_upperbound(bsamps) |> quantile(c(0.025, .975))

round(cbind(psi_lowerbound(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1),
                           terms = TRUE),
            psi_upperbound(matrix(unlist(lapply(data, \(x) x / sum(x))), nrow = 1),
                           terms = TRUE),
            girlsterms), 2)

