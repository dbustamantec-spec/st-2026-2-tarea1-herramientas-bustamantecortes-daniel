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

  

medidas <- function(y, yhat, mae_ingenuo = NULL) {
  stopifnot(
    "y y yhat deben tener la misma longitud" = length(y) == length(yhat)
  )
  idx_validos <- !is.na(yhat) & !is.na(y)
  e <- y[idx_validos] - yhat[idx_validos]
  y_val <- y[idx_validos]
  mse <- mean(e^2)
  mad <- mean(abs(e))
  mape <- mean(abs(e / y_val)) * 100
  
  resultados <- c(MSE = mse, MAD = mad, MAPE = mape)
  if (!is.null(mae_ingenuo)) {
    mase <- mad / mae_ingenuo
    resultados <- c(resultados, MASE = mase)
  }
  
  return(resultados)
}

validar_errores <- function(e, p = 0, m = NULL) {
  library(ggplot2)
  library(patchwork)
  e <- na.omit(as.numeric(e))
  N <- length(e)
  
  if (N < 20) {
    warning("Jarque-Bera es una prueba asintótica. Con menos de 20 errores, la decisión se debe tomar con precaución.")
  }
  
  if (is.null(m)) {
    m <- min(24, floor(N / 4))
  }
  
  # Gráfico de errores en el tiempo
  df_e <- data.frame(t = 1:N, error = e)
  g_tiempo <- ggplot(df_e, aes(x = t, y = error)) +
    geom_line(color = "steelblue") +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Errores de pronóstico a un paso", x = "Tiempo", y = "Error")
  
  # Correlograma 
  g_corr <- correlograma(e, m = m)
  
  # Pruebas de hipótesis 
  t_test <- t.test(e, mu = 0)
  jb_test <- jarque_bera(e)
  r_e <- acf(e, lag.max = m, plot = FALSE)$acf[-1, 1, 1]
  lb_test <- ljung_box(r = r_e, T = N, m = m, p = p)
  dw_test <- durbin_watson(e)

  print(g_tiempo / g_corr)
 
  return(list(
    N_errores = N,
    t_test = t_test,
    ljung_box = lb_test,
    jarque_bera = jb_test,
    durbin_watson = dw_test
  ))
}
  return(list(
    estadistico = d,
    grados_libertad = NA,
    valor_critico_5 = NA,
    p_valor = NA))}
