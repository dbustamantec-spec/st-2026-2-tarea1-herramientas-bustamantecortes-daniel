# Prueba de Ljung-Box
ljung_box <- function(r, T_obs, m, p = 0) {
  # Asegurar tipos numéricos para evitar errores de argumentos
  T_obs <- as.numeric(T_obs)
  m <- as.numeric(m)
  p <- as.numeric(p)
  
  h <- 1:m
  r_m <- r[1:m]
  
  # Fórmula exacta del estadístico Q_m
  Q_m <- T_obs * (T_obs + 2) * sum((r_m^2) / (T_obs - h))
  df <- m - p
  
  val_critico <- qchisq(0.95, df)
  p_val <- pchisq(Q_m, df, lower.tail = FALSE)
  
  return(list(
    estadistico = Q_m,
    grados_libertad = df,
    valor_critico_5 = val_critico,
    p_valor = p_val
  ))
}

# Prueba de Jarque-Bera
jarque_bera <- function(e) {
  N <- length(e)
  e_media <- mean(e)

  m2 <- mean((e - e_media)^2)
  m3 <- mean((e - e_media)^3)
  m4 <- mean((e - e_media)^4)
  
  asimetria <- m3 / (m2^(3/2))
  kurtosis <- m4 / (m2^2)

  JB <- (N / 6) * (asimetria^2 + ((kurtosis - 3)^2) / 4)
  gl <- 2
  
  val_critico <- qchisq(0.95, gl)
  p_val <- pchisq(JB, gl, lower.tail = FALSE)
  
  return(list(
    estadistico = JB,
    grados_libertad = gl,
    valor_critico_5 = val_critico,
    p_valor = p_val))
}

# Prueba de Durbin-Watson
durbin_watson <- function(e) {
  num <- sum(diff(e)^2) 
  den <- sum(e^2)
  d <- num / den
  
  # Como los valores críticos de DW dependen de tablas externas (T y variables), 
  # se devuelve NA en las distribuciones según "cuando la distribución lo permite".
  return(list(
    estadistico = d,
    grados_libertad = NA,
    valor_critico_5 = NA,
    p_valor = NA))}
