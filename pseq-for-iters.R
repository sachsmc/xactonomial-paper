library(xactonomial)

sample_data <- function(n, true_theta) {

  T1 <- rmultinom(1, n[1], prob = true_theta[1:4])
  T2 <- rmultinom(1, n[2], prob = true_theta[5:8])

  list(T1 = c(T1), T2 = c(T2))

}

psi <- function(theta) {

  theta1 <- theta[1:4]
  theta2 <- theta[5:8]
  sum(sqrt(theta1 * theta2))

}

true_theta <- c(.45, .15, .3, .1, .05, .15, .4, .4)
true_psi <- psi(true_theta)

data <- sample_data(n = c(10, 10), true_theta)

k <- length(data)
d_k <- sapply(data, length)
SSpace <- lapply(data, \(x) sspace_multinom(length(x), sum(x)))

bigdex <- expand_index(sapply(SSpace, nrow))
psi_hat <- logC <- rep(NA, nrow(bigdex))
psi_obs <- do.call(psi, lapply(data, \(x) x / sum(x)) |> unlist() |> list())

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

pvalue_psi0 <- function(psi0, maxit = 5000, chunksize = 50,
                        lower = TRUE, target = .05 / 2, psi_limits = c(0,1),
                        SSpacearr = SSpacearr, logC = logC) {

  minus1 <- if(lower) 1 else -1
  II <- if(lower) psi_hat >= psi_obs else psi_hat <= psi_obs

  seqmaxes <- rep(NA, maxit)
  for(i in 1:maxit) {
    theta_cands <- do.call(cbind, lapply(d_k, \(i) get_theta_random(i, chunksize)))
    #theta_cands <- smart_sample_theta(psi, psi0 = psi0,
    #                                  psi_limits = psi_limits, d_k = d_k,
    #                                  chunksize = chunksize, lower = lower)
    these_probs <- calc_prob_null(theta_cands, psi, psi0, minus1,
                                  SSpacearr, logC, II)
    if(length(these_probs) == 0) next

    cand <- c(seqmaxes, these_probs)
    if(all(is.na(cand))) seqmaxes[i] <- 1e-12 else {
      seqmaxes[i] <- max(cand, na.rm = TRUE)
    }
    #if(seqmaxes[i] > target + .001) break

  }
  #if(all(is.na(seqmaxes))) return(1e-12) else max(seqmaxes, na.rm = TRUE)

  seqmaxes[!is.na(seqmaxes)]
}

library(ggplot2)
library(parallel)

psi0 <- c(.2, .3, .4, .5, .6, .7, .8, .9, .95)

blat <- mclapply(seq_along(psi0)[1:3], \(i) {
  thispseq <- pvalue_psi0(psi0[i], maxit = 5000, SSpacearr = SSpacearr,
                          logC = logC)
  data.frame(
  pseq = thispseq,
  iter = 1:length(thispseq),
  psi0 = psi0[i])


}, mc.cores = 3)


res <- do.call(rbind, blat)
ggplot(res, aes(x = iter, y = pseq, color = factor(psi0))) + geom_line() +
  theme_bw() + scale_x_continuous("Number of iterations") +
  scale_y_continuous("p value") + scale_color_discrete(bquote(psi[0]))
ggsave("stab.pdf", width = 6.75, height = 4.25)


data <- sample_data(n = c(10, 10), true_theta = true_theta)
results <- xactonomial(psi, data, alpha = .05, psi_limits = c(0,1))

results$pvalue_function(.5, 5000, 500, target = 1)
