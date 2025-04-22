library(xactonomial)

source("example/speedtest-functions.R")

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

for(i in 1:nrow(nnn)){
data <- sample_data(n = nnn[i,])

time.fast <- system.time(xactonomial(data, psi, psi0 = 0.75, psi_limits = c(0,1), conf_int = FALSE,
                    ga = FALSE, maxit = 100, chunksize = 100))[3]
time.slow <- system.time(xactonomial_slow(data, psi, psi0 = 0.75, psi_limits = c(0,1), conf_int = FALSE,
                                          ga = FALSE, maxit = 100, chunksize = 100))[3]


res <- rbind(res, data.frame(totn = rep(sum(nnn[i,]), 2),
                             time = c(time.fast, time.slow),
                             method = c("base R", "rust")))
cat(i, "\n")
}

library(ggplot2)
ggplot(res, aes(x = totn, y = time / 60, color = method, linetype = method, shape = method)) +
  geom_line() + geom_point() + scale_y_continuous("Time (minutes)") +
  xlab("Total sample size") +
  theme_bw()
ggsave("speedtest.png", width = 5.75, height = 3.5)
