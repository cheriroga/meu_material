# Função para gerar amostra aleatória da distribuição
gerar_amostra <- function(n, theta) {
  # Verifica se os parâmetros estão dentro dos limites especificados
  if (theta < -1 || theta > 1) {
    stop("O parâmetro theta deve estar entre -1 e 1")
  }
  
  # Gera n valores uniformes entre 0 e 1
  u <- runif(n)
  
  # Calcula os valores da amostra usando o método da transformação inversa
  # A função de distribuição acumulada (CDF) é F(x) = (x + 1)/2 + (theta/4)(x^2 - 1)
  # Para inverter, resolvemos a equação quadrática theta*x^2 + 2*x - (1 + theta + 4*u) = 0
  
  amostra <- numeric(n)
  for (i in 1:n) {
    # Coeficientes da equação quadrática
    a <- theta
    b <- 2
    c <- -(1 + theta + 4*u[i])
    
    # Calcula as raízes
    if (theta == 0) {
      # Caso especial quando theta = 0 (distribuição uniforme)
      x <- 2*u[i] - 1
    } else {
      # Caso geral
      discriminant <- b^2 - 4*a*c
      root1 <- (-b + sqrt(discriminant)) / (2*a)
      root2 <- (-b - sqrt(discriminant)) / (2*a)
      
      # Seleciona a raiz que está no intervalo [-1, 1]
      if (root1 >= -1 && root1 <= 1) {
        x <- root1
      } else {
        x <- root2
      }
    }
    amostra[i] <- x
  }
  
  return(amostra)
}

# Exemplo de uso:
# set.seed(123)  # Para reproducibilidade
amostra <- gerar_amostra(n = 1000, theta = -0.8)
hist(amostra, breaks = 30, main = "Histograma da amostra gerada")
