options(stringsAsFactors = FALSE)
set.seed(20260921)

# Cargar las funciones del repositorio
source("R/00-lectura.R")
source("R/01-graficos.R")
source("R/02-metodos.R")
source("R/03-evaluacion.R")

dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# Función para graficar la optimización
graficar_optimizacion <- function(opt, metodo, id) {
  library(ggplot2)
  tab <- opt$rejilla
  optimo <- opt$optimo
  
  if (metodo %in% c("mm", "dmm")) {
    tab <- tab[is.finite(tab$mse), , drop = FALSE]
    g <- ggplot(tab, aes(x = k, y = mse)) + geom_line() + geom_point() + 
      geom_point(data = optimo, color = "red", size = 3) + labs(title = "Optimización de k")
  } else if (metodo == "ses") {
    g <- ggplot(tab, aes(x = alpha, y = mse)) + geom_line() + geom_point() + 
      geom_point(data = optimo, color = "red", size = 3) + labs(title = expression("Optimización de " ~ alpha))
  } else if (metodo == "holt") {
    g <- ggplot(tab, aes(x = alpha, y = beta, fill = mse)) + geom_tile() + 
      geom_point(data = optimo, color = "red", size = 3) + labs(title = "Mapa de MSE (Holt)")
  }
  
  g <- g + theme_minimal()
  ggsave(sprintf("figs/e%d_optimizacion.png", id), g, width = 6, height = 4)
}

ejecutar_ejemplo <- function(id, nombre, x, fuente, unidad, metodo, p_lb, rejilla = NULL) {
  cat("EJEMPLO", id, "-", nombre, "-", metodo, "\n")
  
  # A. Lectura
  datos <- leer_serie(x, fuente, unidad)
  y <- datos$y
  T_total <- length(y)
  
  # B. Gráficos y pruebas de la serie
  g_serie <- graficar_serie(datos, paste("Serie", nombre))
  m_serie <- min(24, floor(T_total / 4))
  g_corr <- correlograma(y, m = m_serie)
  
  ggsave(sprintf("figs/e%d_serie.png", id), g_serie, width = 8, height = 4)
  ggsave(sprintf("figs/e%d_corr_serie.png", id), g_corr, width = 8, height = 6)
  
  cat("\nLjung-Box sobre la serie original (p=0):\n")
  print(ljung_box(acf(y, lag.max=m_serie, plot=FALSE)$acf[-1,1,1], T_total, m_serie, 0))
  
  # Verificación exigida vs R base
  cat("\nVerificación vs Box.test() de R:\n")
  print(Box.test(y, lag = m_serie, type = "Ljung-Box", fitdf = 0))
  
  # C. Partición
  h <- max(1, min(12, floor(0.2 * T_total)))
  n_est <- T_total - h
  y_est <- y[1:n_est]
  y_val <- y[(n_est + 1):T_total]
  
  # D. Ajuste y optimización
  if (metodo %in% c("mm", "ses", "dmm", "holt")) {
    opt <- optimizar(y_est, metodo, rejilla)
    cat("\nParámetros óptimos:\n")
    print(opt$optimo)
    graficar_optimizacion(opt, metodo, id)
  }
  
  fit <- switch(metodo,
                media = ajustar_media(y_est),
                mm = ajustar_mm(y_est, opt$optimo$k),
                ses = ajustar_ses(y_est, opt$optimo$alpha),
                dmm = ajustar_dmm(y_est, opt$optimo$k),
                lineal = ajustar_tendencia(y_est, "lineal"),
                cuadratica = ajustar_tendencia(y_est, "cuadratica"),
                exponencial = ajustar_tendencia(y_est, "exponencial", corregir_sesgo = TRUE),
                holt = ajustar_holt(y_est, opt$optimo$alpha, opt$optimo$beta)
  )
  
  if (metodo %in% c("lineal", "cuadratica", "exponencial")) {
    cat("\nTabla de coeficientes:\n")
    print(fit$parametros$tabla)
  }
  
  # E. Medidas dentro y fuera de muestra
  mae_ingenuo <- mean(abs(diff(y_est)))
  
  cat("\nMedidas en estimación:\n")
  print(medidas(y_est, fit$yhat, mae_ingenuo))
  
  pronostico <- fit$pronosticar(h)
  pronostico_ingenuo <- rep(tail(y_est, 1), h)
  
  cat("\nMedidas extramuestrales (Método):\n")
  m_val <- medidas(y_val, pronostico, mae_ingenuo)
  print(m_val)
  
  cat("\nMedidas extramuestrales (Ingenuo):\n")
  m_ing <- medidas(y_val, pronostico_ingenuo, mae_ingenuo)
  print(m_ing)
  
  # F. Validación de errores
  e_est <- na.omit(y_est - fit$yhat)
  cat("\nValidación de errores a 1 paso:\n")
  val_err <- validar_errores(e_est, p = p_lb, m = min(24, floor(length(e_est)/4)))
  print(val_err[-1]) # Imprime las pruebas ocultando el gráfico bruto
  
  # Guardar gráficos de errores
  # Guardar gráficos de errores forzando el fondo blanco
  ggsave(sprintf("figs/e%d_errores_val.png", id), val_err$grafico, width = 8, height = 6, bg = "white")
  
  # G. Retornar métricas clave para el resumen
  return(data.frame(
    Ejemplo = id, Serie = nombre, Metodo = metodo, 
    MASE_met = m_val["MASE"], MASE_ing = m_ing["MASE"]
  ))
}

# LOS OCHO EJEMPLOS + CONTRAEJEMPLO
resumen <- list()

resumen[[1]] <- ejecutar_ejemplo(1, "Nile", Nile, 
                                 "Durbin and Koopman (2001)", "Nivel", "media", p_lb = 1)

resumen[[2]] <- ejecutar_ejemplo(2, "discoveries", discoveries, 
                                 "World Almanac 1975", "Descubrimientos", "mm", p_lb = 0, rejilla = 2:12)

resumen[[3]] <- ejecutar_ejemplo(3, "Nile", Nile, 
                                 "Durbin and Koopman (2001)", "Nivel", "ses", p_lb = 1, rejilla = seq(0.02, 0.98, by = 0.02))

resumen[[4]] <- ejecutar_ejemplo(4, "airmiles",airmiles, 
                                 "FAA Statistical Handbook", "Millas", "dmm", p_lb = 0, rejilla = 2:12)

resumen[[5]] <- ejecutar_ejemplo(5, "austres", austres, 
                                 "Brockwell and Davis (1996)", "Miles de residentes", "lineal", p_lb = 2)

resumen[[6]] <- ejecutar_ejemplo(6, "austres", austres, 
                                 "Brockwell and Davis (1996)", "Miles de residentes", "cuadratica", p_lb = 3)

resumen[[7]] <- ejecutar_ejemplo(7, "airmiles", airmiles, 
                                 "FAA Statistical Handbook", "Millas", "exponencial", p_lb = 2)

resumen[[8]] <- ejecutar_ejemplo(8, "WWWusage", WWWusage, 
                                 "Durbin J, Koopman SJ (2001", "Usuarios", "holt", p_lb = 2, 
                                 rejilla = expand.grid(alpha = seq(0.05, 0.95, by=0.05), beta = seq(0.05, 0.95, by=0.05)))

# Contraejemplo
resumen[[9]] <- ejecutar_ejemplo(9, "austres", austres, 
                                 "Brockwell and Davis (1996)", "Miles de residentes", "media", p_lb = 1)

# RESUMEN FINAL Y VERIFICACIÓN DE ACF
df_resumen <- do.call(rbind, resumen)
cat("RESUMEN DE MASE (MÉTODO VS INGENUO)\n")
print(df_resumen, row.names = FALSE)
write.csv(df_resumen, "figs/resumen_ejemplos.csv", row.names = FALSE)

cat("VERIFICACIÓN ACF MANUAL vs acf() DE R\n")
series_test <- list(Nile, discoveries, airmiles, austres, WWWusage)
diffs <- sapply(series_test, function(y) {
  m <- min(24, floor(length(y)/4))
  r_r <- acf(y, lag.max = m, plot = FALSE)$acf[-1,1,1]
  
  # Calculo manual
  media <- mean(y)
  den <- sum((y - media)^2)
  r_man <- sapply(1:m, function(h) sum((y[(h+1):length(y)] - media) * (y[1:(length(y)-h)] - media)) / den)
  
  max(abs(r_man - r_r))
})

cat("Diferencias máximas para las 5 series usadas:\n")
print(diffs)
stopifnot("La diferencia de la ACF manual supera 1e-12" = all(diffs < 1e-12))

