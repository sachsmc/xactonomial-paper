
sample_data <- function(n, d, theta) {

  gpsz <- length(theta) / d
  cprod <- rep(1, gpsz)
  j <- 1
  res <- vector("list", length = d)
  for(i in 1:d) {
    res[[i]] <- c(rmultinom(1, n[i], theta[j:(j + gpsz - 1)]))
    j <- j + gpsz
  }
  res
}

psi_bc <- function(theta, d) {

  gpsz <- length(theta) / d
  cprod <- rep(1, gpsz)
  j <- 1
  for(i in 1:d) {
    cprod <- cprod * theta[j:(j + gpsz - 1)]
    j <- j + gpsz
  }

  sum((cprod)^(1/d))

}

psi_lb <- function(theta, d) {

  p00_0 <- theta[1]
  p10_0 <- theta[2]
  p01_0 <- theta[3]
  p11_0 <- theta[4]
  p00_1 <- theta[5]
  p10_1 <- theta[6]
  p01_1 <- theta[7]
  p11_1 <- theta[8]


  max(p00_0 - p00_1 - p10_1 - p01_1, p00_0 - p00_1 - p10_0 - p10_1 - p01_0,
       p00_0 - p00_1 + p10_0 - 2 * p10_1 - 2 * p01_1,
       -p10_1 - p01_1,
       -p10_0 - p01_0,
       -p00_0 + p00_1 - 2 * p10_0 + p10_1 - 2 * p01_0,
       -p00_0 + p00_1 - p10_0 -
        p10_1 - p01_1,
       -p00_0 + p00_1 - p10_0 - p01_0)
}

psi_max <- function(theta, d) {

  max(theta)

}

run_one <- function(i, n, d, theta, psi, limits, maxit = 100, csize = 100) {

  true_psi <- psi(theta, d)
  data <- sample_data(rep(n, d), d, theta)

  bsamps <- replicate(500, {
    unlist(lapply(data, \(x) {
      newx <- table(factor(sample(rep.int(1:length(x), times = x), sum(x), replace = TRUE),
                           levels = 1:length(x)))

      newx / sum(newx)
    })) |> psi(d = d)
  })

  cib <- quantile(bsamps, c(.025, .975))

  psi2 <- function(theta) {
    psi(theta, d = d)
  }
  xacto <- xactonomial(psi2, data, alpha = .05, psi_limits = limits, maxit = maxit, chunksize = csize)

  data.frame(cover_boot = true_psi >= cib[1] & true_psi <= cib[2],
             width_boot = cib[2] - cib[1],
             cover_xact = true_psi >= xacto$conf.int[1] & true_psi <= xacto$conf.int[2],
             width_xact = xacto$conf.int[2] - xacto$conf.int[1],
             i = i, n = n, d = d, theta = paste(theta, collapse = ", "),
             true_psi = true_psi)

}

send_mail <- function(...) {
  system('echo "runs complete" | mail -s "runs finished" michael.sachs@sund.ku.dk')
}
