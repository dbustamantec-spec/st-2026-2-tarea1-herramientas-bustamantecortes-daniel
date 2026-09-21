ajustar_media <- function(y) {
  stopifnot(
    "y debe ser un vector numerico, sin NA" = is.numeric(y) && !anyNA(y),
    "y debe tener al menos 2 observaciones"  = length(y) >= 2)
  Tt <- length(y)
  yhat <- rep(NA_real_, Tt)
  m <- y[1]       
  yhat[2] <- m      
  
  for (t in 2:Tt) {
    m <- m + (y[t] - m) / t   
    if (t < Tt) yhat[t + 1] <- m
  }

  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    rep(m, h)
  }
  
  list(yhat = yhat, pronosticar = pronosticar, parametros = list(media_T = m))
}


# 2. Media movil

ajustar_mm <- function(y, k) {
  stopifnot(
    "y debe ser un vector numerico, sin NA"     = is.numeric(y) && !anyNA(y),
    "k debe ser un entero mayor o igual a 2"    = length(k) == 1 && k >= 2 && k == round(k),
    "la ventana k no puede exceder la longitud de la serie" = length(y) > k
  )
  
  Tt <- length(y)
  yhat <- rep(NA_real_, Tt)
  
  suma <- sum(y[1:k])
  mm <- suma / k
  yhat[k + 1] <- mm            
  
  if (Tt > k + 1) {
    for (t in (k + 1):(Tt - 1)) {
      suma <- suma + y[t] - y[t - k]  
      mm <- suma / k
      yhat[t + 1] <- mm
    }
  }
  if (Tt > k) {
    suma <- suma + y[Tt] - y[Tt - k]
    mm <- suma / k
  }
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    rep(mm, h)
  }
  
  list(yhat = yhat, pronosticar = pronosticar, parametros = list(k = k, mm_T = mm))
}


# 3. Suavizamiento exponencial simple

ajustar_ses <- function(y, alpha) {
  stopifnot(
    "y debe ser un vector numerico, sin NA" = is.numeric(y) && !anyNA(y),
    "y debe tener al menos 2 observaciones"  = length(y) >= 2,
    "alpha debe estar en (0, 1)"             = length(alpha) == 1 && alpha > 0 && alpha < 1
  )
  
  Tt <- length(y)
  yhat <- rep(NA_real_, Tt)
  yhat[2] <- y[1] 
  
  if (Tt >= 3) {
    for (t in 2:(Tt - 1)) {
      e_t <- y[t] - yhat[t]
      yhat[t + 1] <- yhat[t] + alpha * e_t
    }
  }
  
  yhat_T1 <- yhat[Tt] + alpha * (y[Tt] - yhat[Tt])   
  yhat_T1_alt <- alpha * y[Tt] + (1 - alpha) * yhat[Tt]
  stopifnot(
    "la forma de correccion de error no coincide con la de promedio ponderado" =
      isTRUE(all.equal(yhat_T1, yhat_T1_alt))
  )
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    rep(yhat_T1, h) 
  }
  
  list(yhat = yhat, pronosticar = pronosticar, parametros = list(alpha = alpha, yhat_T1 = yhat_T1))
}


#4. Doble media movil

ajustar_dmm <- function(y, k) {
  stopifnot(
    "y debe ser un vector numerico, sin NA"  = is.numeric(y) && !anyNA(y),
    "k debe ser un entero mayor o igual a 2" = length(k) == 1 && k >= 2 && k == round(k),
    "la serie debe tener mas de 2k - 1 observaciones" = length(y) > 2 * k - 1
  )
  
  Tt <- length(y)

  MM <- rep(NA_real_, Tt)
  suma <- sum(y[1:k])
  MM[k] <- suma / k
  if (Tt > k) {
    for (t in (k + 1):Tt) {
      suma <- suma + y[t] - y[t - k]
      MM[t] <- suma / k
    }
  }
  
  DMM <- rep(NA_real_, Tt)
  suma2 <- sum(MM[k:(2 * k - 1)])
  DMM[2 * k - 1] <- suma2 / k
  if (Tt > 2 * k - 1) {
    for (t in (2 * k):Tt) {
      suma2 <- suma2 + MM[t] - MM[t - k]
      DMM[t] <- suma2 / k
    }
  }
  
  Et <- 2 * MM - DMM                     
  beta1t <- (2 / (k - 1)) * (MM - DMM)   
  
  yhat <- rep(NA_real_, Tt)
  for (t in (2 * k - 1):(Tt - 1)) {
    yhat[t + 1] <- Et[t] + beta1t[t]
  }
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    Et[Tt] + beta1t[Tt] * seq_len(h)
  }
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = list(k = k, Et = Et, beta1t = beta1t)
  )
}


# 5. Tendencias por minimos cuadrados

ajustar_tendencia <- function(y, tipo, corregir_sesgo = FALSE) {
  stopifnot(
    "y debe ser un vector numerico, sin NA" = is.numeric(y) && !anyNA(y),
    "tipo debe ser 'lineal', 'cuadratica' o 'exponencial'" =
      length(tipo) == 1 && tipo %in% c("lineal", "cuadratica", "exponencial")
  )
  
  Tt <- length(y)
  t_idx <- seq_len(Tt)
  
  if (tipo == "exponencial") {
    stopifnot(
      "la tendencia exponencial requiere y > 0, porque se ajusta sobre log(y)" = all(y > 0)
    )
    y_aj <- log(y)
    X <- cbind(1, t_idx)
    colnames(X) <- c("a", "theta")
  } else if (tipo == "lineal") {
    y_aj <- y
    X <- cbind(1, t_idx)
    colnames(X) <- c("beta0", "beta1")
  } else {
    y_aj <- y
    X <- cbind(1, t_idx, t_idx^2)
    colnames(X) <- c("beta0", "beta1", "beta2")
  }
  
  p <- ncol(X)
  XtX <- crossprod(X)
  beta <- as.numeric(solve(XtX, crossprod(X, y_aj)))
  
  ajustados <- as.numeric(X %*% beta)
  e <- y_aj - ajustados
  sigma2 <- sum(e^2) / (Tt - p)
  
  ee_ols <- sqrt(diag(sigma2 * solve(XtX)))
  
  L <- floor(4 * (Tt / 100)^(2 / 9))
  Xe <- X * e
  meat <- crossprod(Xe)   # rezago 0
  if (L >= 1) {
    for (l in 1:L) {
      w <- 1 - l / (L + 1)
      Gamma_l <- crossprod(Xe[(l + 1):Tt, , drop = FALSE], Xe[1:(Tt - l), , drop = FALSE])
      meat <- meat + w * (Gamma_l + t(Gamma_l))
    }
  }
  XtX_inv <- solve(XtX)
  V_hac <- XtX_inv %*% meat %*% XtX_inv
  ee_robusto <- sqrt(diag(V_hac))
  
  t_val <- beta / ee_robusto
  p_val <- 2 * stats::pt(-abs(t_val), df = Tt - p)
  
  r2 <- 1 - sum(e^2) / sum((y_aj - mean(y_aj))^2)
  dw <- sum(diff(e)^2) / sum(e^2)
  
  tabla <- data.frame(
    coeficiente = colnames(X),
    estimacion  = beta,
    ee_ols      = ee_ols,
    ee_robusto  = ee_robusto,
    t           = t_val,
    p_valor     = p_val,
    row.names = NULL
  )
  
  yhat <- if (tipo == "exponencial") exp(ajustados) else ajustados
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    t_nuevo <- Tt + seq_len(h)
    if (tipo == "lineal") {
      beta[1] + beta[2] * t_nuevo
    } else if (tipo == "cuadratica") {
      beta[1] + beta[2] * t_nuevo + beta[3] * t_nuevo^2
    } else {
      factor_sesgo <- if (corregir_sesgo) exp(sigma2 / 2) else 1
      exp(beta[1] + beta[2] * t_nuevo) * factor_sesgo
    }
  }
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = list(
      tipo = tipo, tabla = tabla, R2 = r2, sigma2 = sigma2, dw = dw,
      L_hac = L, corregir_sesgo = corregir_sesgo
    )
  )
}


# 6. Holt lineal

ajustar_holt <- function(y, alpha, beta) {
  stopifnot(
    "y debe ser un vector numerico, sin NA" = is.numeric(y) && !anyNA(y),
    "y debe tener al menos 2 observaciones"  = length(y) >= 2,
    "alpha debe estar en (0, 1)"             = length(alpha) == 1 && alpha > 0 && alpha < 1,
    "beta debe estar en (0, 1)"              = length(beta)  == 1 && beta  > 0 && beta  < 1
  )
  
  Tt <- length(y)
  yhat <- rep(NA_real_, Tt)
  L <- rep(NA_real_, Tt)
  Tr <- rep(NA_real_, Tt)   
  
  L[1] <- y[1]   
  Tr[1] <- 0
  
  for (t in 2:Tt) {
    yhat[t] <- L[t - 1] + Tr[t - 1]
    L[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
    Tr[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Tr[t - 1]
  }
  
  e <- y[2:Tt] - yhat[2:Tt]
  L_alt  <- L[1:(Tt - 1)]  + Tr[1:(Tt - 1)] + alpha * e
  Tr_alt <- Tr[1:(Tt - 1)] + alpha * beta * e
  stopifnot(
    "las ecuaciones de Holt no coinciden con su forma de correccion de error" =
      isTRUE(all.equal(L[2:Tt], L_alt)) && isTRUE(all.equal(Tr[2:Tt], Tr_alt))
  )
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h > 0 && h == round(h))
    L[Tt] + Tr[Tt] * seq_len(h)
  }
  
  list(
    yhat = yhat,
    pronosticar = pronosticar,
    parametros = list(alpha = alpha, beta = beta, L = L, Tr = Tr)
  )
}

# 7. Optimización de constantes y ventanas
optimizar <- function(y, metodo, rejilla) {
  stopifnot(
    "y debe ser un vector numerico, sin NA" = is.numeric(y) && !anyNA(y),
    "metodo debe ser uno de: mm, dmm, ses, holt" =
      length(metodo) == 1 && metodo %in% c("mm", "dmm", "ses", "holt")
  )
  
  mse_de <- function(ajustar_uno) {
    tryCatch(mean((y - ajustar_uno()$yhat)^2, na.rm = TRUE), error = function(e) NA_real_)
  }
  
  if (metodo %in% c("mm", "dmm")) {
    ajustar <- if (metodo == "mm") ajustar_mm else ajustar_dmm
    mse <- sapply(rejilla, function(k) mse_de(function() ajustar(y, k)))
    tabla <- data.frame(k = rejilla, mse = mse)
    
  } else if (metodo == "ses") {
    mse <- sapply(rejilla, function(a) mse_de(function() ajustar_ses(y, a)))
    tabla <- data.frame(alpha = rejilla, mse = mse)
    
  } else {
    stopifnot(
      "rejilla debe ser un data.frame con columnas alpha y beta para holt" =
        is.data.frame(rejilla) && all(c("alpha", "beta") %in% names(rejilla))
    )
    mse <- mapply(function(a, b) mse_de(function() ajustar_holt(y, a, b)),
                  rejilla$alpha, rejilla$beta)
    tabla <- data.frame(alpha = rejilla$alpha, beta = rejilla$beta, mse = mse)
  }
  
  optimo <- tabla[which.min(tabla$mse), , drop = FALSE]
  list(rejilla = tabla, optimo = optimo)
}
