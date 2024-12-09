library(xactonomial)

xactonomial_speedtest <- function(psi, data, psi0 = NULL, alpha = .05, psi_limits,
                        maxit = 50, chunksize = 500) {

  k <- length(data)
  d_k <- sapply(data, length)
  psi_obs <- do.call(psi, lapply(data, \(x) x / sum(x)) |> unlist() |> list())


  SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

  bigdex <- expand_index(sapply(SSpace, nrow))
  psi_hat <- logC <- rep(NA, nrow(bigdex))

  SSpacearr <- array(dim = c(nrow(bigdex), sum(d_k)))

  for(i in 1:nrow(bigdex)) {

    thisS <- lapply(1:length(bigdex[i,]), \(j){
      Sj <- SSpace[[j]][bigdex[i, j],]
      Sj
    })
    SSpacearr[i,] <- unlist(thisS)
    psi_hat[i] <- do.call(psi, lapply(thisS, \(x) x / sum(x)) |> unlist() |> list())
    logC[i] <- sum(sapply(thisS, \(x) log_multinom_coef(x, sum(x))))

  }

  slowtime <- system.time({
    pvalue_psi0_slow(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                lower = TRUE, target = 1,
                SSpacearr = SSpacearr, logC = logC, d_k = d_k)
  })[1]

  fasttime <- system.time({
    xactonomial:::pvalue_psi0(psi0 = psi0, psi = psi, psi_hat = psi_hat,
                     psi_obs = psi_obs, maxit = maxit, chunksize = chunksize,
                     lower = TRUE, target = 1,
                     SSpacearr = SSpacearr, logC = logC, d_k = d_k)
  })[1]

  c(slowtime, fasttime)



}


#' Compute a p value for the test of psi <= psi0
#'
#' @param psi0 The null value
#' @param psi The function of interest
#' @param psi_hat The vector of psi values at each element of the sample space
#' @param psi_obs The observed estimate
#' @param maxit Maximum iterations
#' @param chunksize Chunk size
#' @param lower Do a one sided test of the null that it is less than psi0, otherwise greater.
#' @param target Stop the algorithm if p >= target (for speed)
#' @param SSpacearr The sample space array
#' @param logC The log multinomial coefficient
#' @param d_k The vector of dimensions
#' @returns A p-value
#'

pvalue_psi0_slow <- function(psi0, psi, psi_hat, psi_obs, maxit, chunksize,
                        lower = TRUE, target,
                        SSpacearr, logC, d_k) {

  minus1 <- if(lower) 1 else -1
  II <- if(lower) psi_hat >= psi_obs else psi_hat <= psi_obs

  seqmaxes <- rep(NA, maxit)
  for(i in 1:maxit) {
    theta_cands <- do.call("cbind", lapply(d_k, \(i) get_theta_random(i, chunksize)))

    these_probs <- calc_prob_null(theta_cands, psi, psi0, minus1,
                                   SSpacearr, logC, II)

    if(length(these_probs) == 0) next

    cand <- c(seqmaxes, these_probs)
    if(all(is.na(cand))) seqmaxes[i] <- 1e-12 else {
      seqmaxes[i] <- max(cand, na.rm = TRUE)
    }
    #if(seqmaxes[i] > target + .001) break

  }
  if(all(is.na(seqmaxes))) return(1e-12) else max(seqmaxes, na.rm = TRUE)
}


true_theta <- c(.45, .15, .3, .1, .05, .15, .4, .4)

sample_data <- function(n) {

  T1 <- rmultinom(1, n[1], prob = true_theta[1:4])
  T2 <- rmultinom(1, n[2], prob = true_theta[5:8])

  list(T1 = c(T1), T2 = c(T2))

}

psi <- function(theta) {

  theta1 <- theta[1:4]
  theta2 <- theta[5:8]
  sum(sqrt(theta1 * theta2))

}

psi(true_theta)
set.seed(2024)

nnn <- cbind(seq(5, 30, by = 5), seq(5, 30, by = 5))
res <- NULL

for(i in 5:nrow(nnn)){
data <- sample_data(n = nnn[i,])

time <- xactonomial_speedtest(psi, data, alpha = .05,
                      psi0 = 0.75, psi_limits = c(0,1),
                      maxit = 10, chunksize = 100)


res <- rbind(res, data.frame(totn = rep(sum(nnn[i,]), 2),
                             time = time,
                             method = c("base R", "rust")))
cat(i, "\n")
}

library(ggplot2)
ggplot(res, aes(x = totn, y = time, color = method)) +
  geom_line() + geom_point() + scale_y_continuous("Time (seconds)", trans = "log2") +
  xlab("Total sample size") +
  theme_bw()
ggsave("speedtest.png", width = 5.75, height = 3.5)
