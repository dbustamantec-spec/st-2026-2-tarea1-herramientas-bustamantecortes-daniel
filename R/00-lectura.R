leer_serie <- function(x, fuente, unidad) {
  library(tibble)
  # Si es ts
  if (is.ts(x)) {
    freq <- frequency(x)
    año_inicio <- start(x)[1]
    periodo_inicio <- start(x)[2]
    
    # Calcular el mes de inicio 
    mes_inicio <- ifelse(freq == 4, (periodo_inicio - 1) * 3 + 1, periodo_inicio)
    fecha_inicio <- as.Date(paste(año_inicio, mes_inicio, "1", sep = "-"))
    
    # Determinar si la secuencia es de meses, trimestres o años
    paso <- ifelse(freq == 12, "month", ifelse(freq == 4, "quarter", "year"))
    datos <- tibble(
      t = 1:length(x),
      fecha = seq(from = fecha_inicio, by = paso, length.out = length(x)), 
      y = as.numeric(x)
    )
    frecuencia_val <- freq
    
  # Si es csv  
  } else {
    csv <- read.csv(x)
    datos <- tibble(
      t = 1:nrow(csv),
      fecha = as.Date(csv$fecha),
      y = as.numeric(csv$valor)
    )
    frecuencia_val <- round(365 / as.numeric(datos$fecha[2] - datos$fecha[1]))
  }
  
  diferencias <- as.numeric(diff(datos$fecha))
  
  if (any(diferencias <= 0)) {
    stop("Las fechas no son crecientes.")
  }
  if (max(diferencias) - min(diferencias) > 4) {
    stop("Las fechas no están equiespaciadas.")
  }
  
  attr(datos, "frecuencia") <- frecuencia_val
  attr(datos, "fuente") <- fuente
  attr(datos, "unidad") <- unidad
  
  return(datos)
}
