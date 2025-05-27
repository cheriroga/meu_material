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
    modelo <- lm(Y ~ X)
    r2 <- summary(modelo)$r.squared
    c(coef(modelo), r2)
  })
  list(
    beta0_hat = resultados[1, ],
    beta1_hat = resultados[2, ],
    r2 = resultados[3, ]
  )
}

resultados_simulacao <- lapply(ns, simular_estimativas)
names(resultados_simulacao) <- as.character(ns)

# Plotando a comparação dos Betas estimados com a normal

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

for (n in ns) {
  res <- resultados_simulacao[[as.character(n)]]
  
  qqPlot(res$beta0_hat,
         main = bquote("QQ Plot de " ~ hat(beta)[0] ~ " (n = " ~ .(n) ~ ")"),
         ylab = "Quantis da amostra", xlab = "Quantis teóricos", col = "black")
}

for (n in ns) {
  res <- resultados_simulacao[[as.character(n)]]
  
  qqPlot(res$beta1_hat,
         main = bquote("QQ Plot de " ~ hat(beta)[1] ~ " (n = " ~ .(n) ~ ")"),
         ylab = "Quantis da amostra", xlab = "Quantis teóricos", col = "black")
}

# Distribuição Beta ajustada aos R² simulados

# Função de verossimilhança negativa para a distribuição Beta
log_veross_neg <- function(par, x) {
  a <- par[1]
  b <- par[2]
  if (a <= 0 || b <= 0) return(Inf)
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Ajustando a distribuição Beta aos R² simulados
ajustar_beta_para_r2 <- function(r2) {
  r2 <- r2[r2 > 0 & r2 < 1]  

  mu <- mean(r2)
  sigma2 <- var(r2)
  a0 <- mu * ((mu * (1 - mu)) / sigma2 - 1)
  b0 <- (1 - mu) * ((mu * (1 - mu)) / sigma2 - 1)

  ajuste <- optim(
    par = c(a0, b0),
    fn = log_veross_neg,
    x = r2,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),
    hessian = TRUE
  )

  a_hat <- ajuste$par[1]
  b_hat <- ajuste$par[2]

  hess_inv <- solve(ajuste$hessian)
  se <- sqrt(diag(hess_inv))
  IC_a <- a_hat + c(-1.96, 1.96) * se[1]
  IC_b <- b_hat + c(-1.96, 1.96) * se[2]

  return(list(a = a_hat, b = b_hat, IC_a = IC_a, IC_b = IC_b))
}

# Ajustando a distribuição Beta para cada n
resultados_ajuste <- lapply(as.character(ns), function(n) {
  r2 <- resultados_simulacao[[n]]$r2
  ajustar_beta_para_r2(r2)
})
names(resultados_ajuste) <- as.character(ns)

# Tabela com os parâmetros ajustados e seus intervalos de confiança
tabela_resultados <- do.call(rbind, lapply(names(resultados_ajuste), function(n) {
  ajuste <- resultados_ajuste[[n]]
  data.frame(
    n = as.integer(n),
    a = ajuste$a,
    a_lower = ajuste$IC_a[1],
    a_upper = ajuste$IC_a[2],
    b = ajuste$b,
    b_lower = ajuste$IC_b[1],
    b_upper = ajuste$IC_b[2]
  )
}))

print(tabela_resultados, row.names = FALSE)


# Plotando a comparação dos R² simulados com a distribuição Beta ajustada
library(car)

qqPlot(
  r2_simulado,
  distribution = "beta",
  shape1 = a_hat,
  shape2 = b_hat,
  main = paste("QQ Plot para R² simulado (n =", n_escolhido, ")"),
  xlab = "Quantis Teóricos da Beta",
  ylab = "Quantis Amostrais do R²",
  col = "black",
  pch = 19
)

plot_qq_r2 <- function(n_escolhido, resultados_simulacao, resultados_ajuste) {
  r2_simulado <- resultados_simulacao[[as.character(n_escolhido)]]$r2
  r2_simulado <- r2_simulado[r2_simulado > 0 & r2_simulado < 1]

  parametros_beta <- resultados_ajuste[[as.character(n_escolhido)]]
  a_hat <- parametros_beta$a
  b_hat <- parametros_beta$b

  qqPlot(
    r2_simulado,
    distribution = "beta",
    shape1 = a_hat,
    shape2 = b_hat,
    main = paste("QQ Plot para R² simulado (n =", n_escolhido, ")"),
    xlab = "Quantis Teóricos da Beta",
    ylab = "Quantis Amostrais do R²",
    col = "black",
    pch = 19
  )
}

par(mfrow = c(2, 2))

for (n in ns) {
  plot_qq_r2(n, resultados_simulacao, resultados_ajuste)
}
