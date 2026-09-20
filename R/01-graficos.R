# Serie
graficar_serie <- function(datos, titulo) {
  library(ggplot2)
  unidad <- attr(datos, "unidad")
  fuente <- attr(datos, "fuente")
  n_obs <- nrow(datos)

  g <- ggplot(datos, aes(x = fecha, y = y)) +
    geom_line(color = "black", linewidth = 0.5) +
    labs(title = titulo, x = "Fecha",y = sprintf("%s", unidad),
      caption = sprintf("Fuente: %s | Número de observaciones: %d", fuente, n_obs)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(hjust = 0, color = "black", size = 9))
  
  return(g)
}

graficar_serie(datos = leer_serie(x = nottem, unidad = "Temperatura", 
                                  fuente = "quien sabe"), titulo = "hola")
not <- leer_serie(x = nottem, unidad = "Temperatura", fuente = "quien sabe")
not$y$n.used
length(not$t)

# Correlogramas
correlograma <- function(datos, m = NULL) {
  library(patchwork)
  library(tidyverse)
  y_vec <- if(is.data.frame(datos)) datos$y else as.numeric(datos)
  n <- length(y_vec)
  
  if(is.null(m)) {
    m <- min(24, floor(n/ 4))}
  
  media <- mean(y_vec)
  acfs <- numeric(m)
  denominador <- sum((y_vec - media)^2)
  
  ci <- 0.95
  climo <- qnorm((1 + ci)/2) / sqrt(n)
  banda_inf <- -climo
  banda_sup <- climo
  
  for (j in 1:m) {
    numerador <- 0
    for (i in (j+1):n) {
      numerador <- numerador + (y_vec[i] - media) * (y_vec[i-j] - media)
    }
    acfs[j] <- numerador / denominador
  }

  pacfs <- pacf(y_vec, lag.max = m, plot = FALSE)
  pacfs <- as.vector(pacfs$acf)
  df <- data.frame(Rezagos = 1:m, acfs = acfs, pacfs = pacfs)
  
  # Correlogramas
  g1 <- ggplot(df, aes(x = Rezagos, y = acfs)) +
    geom_segment(aes(xend = Rezagos, yend = 0), color = "black", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_hline(yintercept = c(banda_inf, banda_sup), color = "red", linetype = "dashed") + 
    labs(title ="Correlograma ACF", x = "Rezago", y = "Autocorrelación") +
    theme_minimal()
  
  g2 <- ggplot(df, aes(x = Rezagos, y = pacfs)) +
    geom_segment(aes(xend = Rezagos, yend = 0), color = "black", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_hline(yintercept = c(banda_inf, banda_sup), color = "red", linetype = "dashed") + 
    labs(title ="Correlograma PACF", x = "Rezago", y = "Autocorrelación") +
    theme_minimal()
  
  attr(acfs, "acfs") <- acfs
  return(g1 / g2)

}

# Bloque de verificación con Nottem
y_vec <- leer_serie(nottem, fuente = "", unidad = "")
y_vec <- y_vec$y 
n <- length(y_vec)
m <- min(24, floor(n / 4))
media <- mean(y_vec)
acfs <- numeric(m)
denominador <- sum((y_vec - media)^2)
ci <- 0.95
climo <- qnorm((1 + ci)/2) / sqrt(n)
banda_inf <- -climo
banda_sup <- climo

for (j in 1:m) {
  numerador <- 0
  for (i in (j+1):n) {
    numerador <- numerador + (y_vec[i] - media) * (y_vec[i-j] - media)
  }
  acfs[j] <- numerador / denominador
}

acf_r <- acf(y_vec, plot = FALSE, lag.max = m)$acf[-1]
max_diferencia <- max(abs(acfs - acf_r))
cat("Máxima diferencia absoluta:", max_diferencia)
# La diferencia es menor a 10^(-12)
