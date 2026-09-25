https://github.com/dbustamantec-spec/st-2026-2-tarea1-herramientas-bustamantecortes-daniel.git

# Tarea 1: caja de herramientas de pronóstico

Implementación en R de ocho métodos de pronóstico de un paso para series de tiempos
univariados. 

## ¿Qué contiene el repositorio?

La estructura esperada del repositorio es la siguiente:

```text
st-2026-2-tarea1-herramientas-apellidos-nombres/
├── README.md
├── .gitignore
├── R/
│   ├── 00-lectura.R
│   ├── 01-graficos.R
│   ├── 02-metodos.R
│   └── 03-evaluacion.R
├── ejemplos/
│   └── ejemplos.R
├── informe/
│   ├── informe.qmd
│   └── informe.html  (o informe.pdf)
├── figs/
└── sesion-info.txt
```

| Archivo o carpeta | Contenido | Dependencias |
|---|---|---|
| `R/00-lectura.R` | `leer_serie()`: lectura de objetos `ts` o CSV y validación de fechas | R base, `tibble` |
| `R/01-graficos.R` | `graficar_serie()` y `correlograma()` | `ggplot2`, `patchwork` |
| `R/02-metodos.R` | Media simple, media móvil, SES, doble media móvil, tendencias lineal/cuadrática/exponencial, Holt y `optimizar()` | R base |
| `R/03-evaluacion.R` | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()` y `validar_errores()` | R base, `ggplot2`, `patchwork` |
| `ejemplos/ejemplos.R` | Ejecución automática de los ocho ejemplos y el contraejemplo | Archivos de `R/` y paquetes permitidos |
| `informe/informe.qmd` | Descripción, gráficos, pruebas, evaluación y conclusiones | Quarto, Typst y paquetes del proyecto |
| `figs/` | Gráficos y `resumen_ejemplos.csv` generados por `ejemplos.R` | Se crea o actualiza al ejecutar |
| `sesion-info.txt` | Salida de `sessionInfo()` usada para cerrar la entrega | R base |

Los únicos paquetes de R usados por el proyecto son `tibble`, `dplyr`, `tidyr`,
`purrr`, `ggplot2` y `patchwork`, además de las funciones de R base y `stats`
permitidas en la consigna. Los métodos no usan `forecast`, `fable`, `smooth`,
`TTR`, `zoo`, `HoltWinters()`, `filter()`, `lm()`, `decompose()` ni `stl()`.

## ¿Cómo se corre?

Se debe ejecutar el siguiente comando desde la raíz del repositorio:

```r
source("ejemplos/ejemplos.R")
```

El script fija la semilla `20260921`, carga las funciones desde `R/`, ejecuta los
ejemplos en orden y guarda las figuras en `figs/`. Al finalizar también genera:

```text
figs/resumen_ejemplos.csv
```

Para renderizar el informe desde la carpeta `informe/`:

```r
quarto_render("informe/informe.qmd")
```

```text
Tiempo observado de ejemplos.R: 1.2048 minutos
Tiempo observado para renderizar informe.qmd: 49.86 segundos
```

Después de cerrar la entrega, guardar la salida de:

```r
sink("sesion-info.txt")
sessionInfo()
sink()
```

## ¿Cómo se usan las funciones?

Ejemplo mínimo de ajuste, evaluación y pronóstico con suavizamiento exponencial
simple:

```r
source("R/02-metodos.R")
source("R/03-evaluacion.R")

y <- Nile
ajuste <- ajustar_ses(y, alpha = 0.30)

# Pronósticos a los tres períodos siguientes
ajuste$pronosticar(3)

# Errores de un paso disponibles después del calentamiento
errores <- na.omit(y - ajuste$yhat)
medidas(y, ajuste$yhat)
```

Todos los métodos reciben una serie numérica ordenada y sin valores faltantes y
devuelven una lista con `yhat`, `pronosticar` y `parametros`. Las constantes y
ventanas se seleccionan con `optimizar()` cuando corresponde:

```r
source("R/02-metodos.R")

rejilla <- seq(0.02, 0.98, by = 0.02)
optimo <- optimizar(Nile, "ses", rejilla)
fit <- ajustar_ses(Nile, optimo$optimo$alpha)
fit$pronosticar(h = 3)
```

## Convenciones que fijan los números

- La validación extramuestral usa las últimas `h = min(12, floor(0.2*T))`
  observaciones de cada serie.
- El referente ingenuo repite la última observación del tramo de estimación.
- El MASE usa como escala el error absoluto medio del pronóstico ingenuo a un
  paso dentro del tramo de estimación.
- La media simple se actualiza recursivamente: `Yhat[2] = Y[1]`.
- La media móvil y la doble media móvil usan ventanas `k >= 2`; la doble media
  móvil comienza después de `2*k - 1` períodos.
- En SES se inicializa `Yhat[2] = Y[1]` y se usa `0 < alpha < 1`.
- Holt se inicializa con `L[1] = Y[1]`, `T[1] = 0`, `0 < alpha < 1` y
  `0 < beta < 1`.
- Las tendencias se estiman mediante las ecuaciones normales
  `solve(crossprod(X), crossprod(X, y))`, sin `lm()`. La tendencia exponencial
  se ajusta sobre `log(y)` y permite la corrección de sesgo
  `exp(sigma2 / 2)`.
- Los errores estándar robustos de las tendencias usan el núcleo de Bartlett y
  `L = floor(4*(T/100)^(2/9))` rezagos.
- La ACF se calcula manualmente con un único divisor
  `sum((y - mean(y))^2)` para todos los rezagos. La banda de los gráficos es
  `qnorm(0.975) / sqrt(n)`, donde `n` es el número de observaciones de la
  sucesión graficada.
- El Ljung-Box usa `Q = T*(T+2)*sum(r_h^2/(T-h))` y `m-p` grados de libertad.
  Jarque-Bera usa dos grados de libertad. Durbin-Watson se reporta como
  estadístico descriptivo cuando no se dispone de sus cotas tabuladas.

## Resumen de resultados

La siguiente tabla resume los ocho ejemplos del informe. Los valores definitivos
de parámetros y MASE se generan sin intervención al ejecutar
`source("ejemplos/ejemplos.R")` y se guardan en
`figs/resumen_ejemplos.csv`. Deben copiarse aquí desde ese archivo antes de
etiquetar la versión `v1.0`.

| Ejemplo | Serie | Método y parámetros | MASE método | MASE ingenuo |
|---:|---|---|---:|---:|
| 1 | `LakeHuron` | Media simple | 0.8429 | 0.8355 |
| 2 | `discoveries` | Media móvil, `k` óptimo en `9` | 0.5878 | 1.1365 |
| 3 | `Nile` | SES, `alpha` óptimo en `0.24` | 0.8062 | 0.8355 |
| 4 | `airmiles` | Doble media móvil, `k` óptimo en `2` | 1.6179 | 4.4959 |
| 5 | `austres` | Tendencia lineal | 3.9735 | 6.091 |
| 6 | `austres` | Tendencia cuadrática | 1.4549 | 6.091 |
| 7 | `airmiles` | Tendencia exponencial, con corrección de sesgo | 20.8072 | 4.4959 |
| 8 | `LakeHuron` | Holt, `(alpha, beta)` óptimos en la rejilla `(0.95,0.95)` | 1.8903 | 7.6247 |

El contraejemplo adicional usa `austres` con media simple: la media histórica
no representa bien una serie con crecimiento sostenido y debe interpretarse como
una línea base, no como el modelo recomendado.

## Declaración de uso de IA

Se utilizó GitHub Copilot para apoyar la organización y redacción de este 
README a partir de la consigna, los scripts R y el informe del proyecto. Se
recibió una propuesta de estructura, explicación de comandos, convenciones y
tabla de resultados. Se verificó la correspondencia de cada
sección con la rúbrica, se revisó que los nombres de archivos y funciones
coincidieran con la consigna.s
