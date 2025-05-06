RNGkind("Marsaglia")
set.seed(1234)

ns <- c(10,50,100,1000)
beta0 <- 3
beta1 <- 5
B <- 10000

beta0_hat <- numeric(B)
beta1_hat <- numeric(B)

simular_estimativas <- function(n) {
  resultados <- replicate(B, {
    X <- runif(n)
    e <- rnorm(n)
    Y <- beta0 + beta1 * X + e
    coef(lm(Y ~ X))
  })
  list(beta0_hat = resultados[1, ], beta1_hat = resultados[2, ])
}

# usando for
# for (n in ns) {
#   resultados_simulacao[[as.character(n)]] <- simular_estimativas(n, B)
# }

# usando lapply
resultados_simulacao <- lapply(ns, simular_estimativas)
names(resultados_simulacao) <- as.character(ns)

# usando for pra plotar os graficos

par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))  

for (n in ns) {
  res <- resultados_simulacao[[as.character(n)]]
  
  # Histograma beta0_hat
  hist(res$beta0_hat, probability = TRUE, col = "lightblue",
       main = bquote(hat(beta)[0] ~ " (n = " ~ .(n) ~ ")"),
       xlab = expression(hat(beta)[0]), breaks = 50)
  curve(dnorm(x, mean = mean(res$beta0_hat), sd = sd(res$beta0_hat)),
        col = "red", lwd = 2, add = TRUE)
  
  # Histograma beta1_hat
  hist(res$beta1_hat, probability = TRUE, col = "lightgreen",
       main = bquote(hat(beta)[1] ~ " (n = " ~ .(n) ~ ")"),
       xlab = expression(hat(beta)[1]), breaks = 50)
  curve(dnorm(x, mean = mean(res$beta1_hat), sd = sd(res$beta1_hat)),
        col = "red", lwd = 2, add = TRUE)
}


library(car)
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))

for (n in ns) {
  res <- resultados_simulacao[[as.character(n)]]
  
  qqPlot(res$beta0_hat, 
         main = bquote("Q-Q de " ~ hat(beta)[0] ~ " (n = " ~ .(n) ~ ")"),
         ylab = "Quantis da amostra", xlab = "Quantis teóricos", col = "black")
  
  qqPlot(res$beta1_hat, 
         main = bquote("Q-Q de " ~ hat(beta)[1] ~ " (n = " ~ .(n) ~ ")"),
         ylab = "Quantis da amostra", xlab = "Quantis teóricos", col = "black")
}
