library(xactonomial)

## comment 1, regularized max


sample_data2 <- function(n) {

  list(rmultinom(1, n, prob = rep(1/4, 4))[, 1])

}

data <- list(c(13, 24, 13))

psi_max <- function(pp) {

  max(pp)

}

pees <- cover <- rep(NA, 1000)
for(i in 1:length(pees)) {
  data <- sample_data2(65)
  runs <- xactonomial(data, psi_max, psi_limits = c(1 / 4, 1), psi0 = 1/ 4,
                      theta_null_points = t(c(1 / 4, 1 / 4, 1 / 4, 1 / 4)),
                      alternative = "two.sided",
                         p_target = .3, conf_int = FALSE)
  pees[i] <- runs$p.value

}
mean(pees < .05)

statistic_max <- function(df) {

  if(is.matrix(df)) {
    pmat <- df + matrix(rep(sqrt(rowSums(df)), ncol(df)), nrow = nrow(df), ncol = ncol(df))
    apply(pmat, 1, \(re) max(re / sum(re)))
  } else {
    pvec <- sqrt(sum(df)) + df
    max(pvec / sum(pvec))
  }

}

sspace <- lapply(data, \(x) matrix(sspace_multinom(length(x), sum(x)), ncol = length(x), byrow = TRUE))

plot(statistic_max(sspace[[1]]), apply(sspace[[1]], 1, \(x) psi_max(x / sum(x))))

xactonomial(data, psi_max, psi_limits = c(1 / 3, 1), psi0 = 1 / 2,
            statistic = statistic_max)

xactonomial(data, psi_max, psi_limits = c(1 / 3, 1), psi0 = 1 / 3, conf_int = TRUE)

## other point about multiple roots
## two sample binomial

data <- list(c(7, 262 - 7), c(30, 494 - 30))
psi_binom <- function(theta) {

  (theta[, 1] - theta[, 3])

}


set.seed(123)
psi.tests <- seq(0.031, .033, length.out = 10)
pees <- rep(NA, length(psi.tests))
for(i in 1:length(psi.tests)) {
   test <- xactonomial(data, psi_binom, psi_limits = c(-1, 1), psi0 = psi.tests[i], alternative = "less",
                        conf_int = FALSE, maxit = 250, chunksize = 100, p_target = .025,
                        ga = TRUE, ga_restart_every = 50)
  pees[i] <- test$p.value
}

png("../nonmonotone.png")
plot(pees ~ psi.tests, ylim = c(0.01, .04), type = "p", pch = 20,
     xlab = "psi0", ylab = "p-value, one-sided vs less")
abline(h = .025, lty = 3)
dev.off()

iter1 <- xactonomial(data, psi_binom, psi_limits = c(-1, 1))
iter1
xactonomial(data, psi_binom, psi_limits = c(iter1$conf.int[2], 1))



xactonomial(psi_binom, data, psi_limits = c(0, 1),
            conf.int = TRUE, maxit = 500, chunksize = 2000, psi_is_vectorized = TRUE)
## point 2 -- hard to sample from problem

psi <- function(theta) {

  theta[, 1]

}

data <- list(c(20, rep(0, 9)))

set.seed(123)
ascent <- xactonomial(data = data, psi, psi0 = .1 ^(1 / 20), conf_int = FALSE, alternative = "greater",
                      psi_limits = c(0, 1), chunksize = 10, maxit = 10000,
                      ga_restart_every = 100,
                      theta_sampler = runif_dk_vects, p_target = .099,
                      ga = TRUE, ga_gfactor = "adapt", ga_lrate = .01)
set.seed(123)
posterior <- xactonomial( data = data, psi, psi0 = .1 ^(1 / 20), conf_int = FALSE, alternative = "greater",
                         psi_limits = c(0, 1), chunksize = 10, maxit = 10000,
                         ga_restart_every = 50,
                         theta_sampler = \(d_k, chunksize) rdirich_dk_vects(chunksize, list(data[[1]] + 1)),
                         p_target = .099,
                         ga = FALSE)
set.seed(123)
combination <- xactonomial(data = data, psi, psi0 = .1 ^(1 / 20), conf_int = FALSE, alternative = "greater",
                      psi_limits = c(0, 1), chunksize = 10, maxit = 10000,
                      ga_restart_every = 100,
                      theta_sampler = \(d_k, chunksize) rdirich_dk_vects(chunksize, list(data[[1]] + 1)),
                      p_target = .099, ga_lrate = .01,
                      ga = TRUE, ga_gfactor = "adapt")

set.seed(123)
uniform <- xactonomial(data = data, psi, psi0 = .1 ^(1 / 20), conf_int = FALSE, alternative = "greater",
                       psi_limits = c(0, 1), chunksize = 10, maxit = 10000,
                       theta_sampler = runif_dk_vects, p_target = .099,
                       ga = FALSE)


library(ggplot2)

pdat <- rbind(data.frame(iteration = seq_along(uniform$p.sequence$p.null),
                 p.value = uniform$p.sequence$p.null, method = "uniform"),
      data.frame(iteration = seq_along(combination$p.sequence$p.null),
                 p.value = combination$p.sequence$p.null, method = "combination"),
      data.frame(iteration = seq_along(posterior$p.sequence$p.null),
                 p.value = posterior$p.sequence$p.null, method = "posterior"),
      data.frame(iteration = seq_along(ascent$p.sequence$p.null),
                 p.value = ascent$p.sequence$p.null, method = "gradient ascent"))
pdat$method <- factor(pdat$method, levels = c("uniform", "posterior", "gradient ascent", "combination"))
ggplot(pdat,
      aes(x = iteration * 10, y = p.value, color = method)) + geom_step() + scale_x_log10() + theme_bw() +
  theme(legend.position = "bottom") + xlab("iteration")+
  geom_hline(yintercept = .1) + scale_color_brewer(type = "qual", palette = "Dark2")

ggsave("../method-compare.png", width = 5.75, height = 4.25)



## speed of calculating combinatinos

k <- 5
n <- 200
choose(k + n - 1, k - 1) / system.time(sspace_multinom(k, n))[3]
