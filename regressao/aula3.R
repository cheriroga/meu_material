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

log_veross_neg <- function(par, x) {
  a <- par[1]
  b <- par[2]
  if (a <= 0 || b <= 0) return(Inf)
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Usando o método da matriz H
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

resultados_ajuste <- lapply(as.character(ns), function(n) {
  r2 <- resultados_simulacao[[n]]$r2
  ajustar_beta_para_r2(r2)
})
names(resultados_ajuste) <- as.character(ns)

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



# Usando da função fitdistr() do pacote MASS
library(MASS)

parametros_mle <- lapply(resultados_simulacao, function(res) {
  r2_vals <- res$r2
  r2_vals <- r2_vals[r2_vals > 0 & r2_vals < 1]  
 
  fit <- fitdistr(
    r2_vals,
    dbeta,
    start = list(shape1 = 2, shape2 = 2),
    lower = c(0.01, 0.01)
  )
 
  est <- coef(fit)
  se <- fit$sd

  ic_shape1 <- est["shape1"] + c(-1.96, 1.96) * se["shape1"]
  ic_shape2 <- est["shape2"] + c(-1.96, 1.96) * se["shape2"]

  list(
    shape1 = est["shape1"],
    shape1_ic = ic_shape1,
    shape2 = est["shape2"],
    shape2_ic = ic_shape2
  )
})

tabela_parametros <- do.call(rbind, lapply(names(parametros_mle), function(n) {
  p <- parametros_mle[[n]]
  data.frame(
    n = as.integer(n),
    shape1 = p$shape1,
    shape1_lower = p$shape1_ic[1],
    shape1_upper = p$shape1_ic[2],
    shape2 = p$shape2,
    shape2_lower = p$shape2_ic[1],
    shape2_upper = p$shape2_ic[2]
  )
}))

print(tabela_parametros, row.names = FALSE)