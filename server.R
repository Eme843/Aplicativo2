library(shiny)
library(kableExtra)
library(data.table)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(quarto)
library(DT)


#============= Emely
# Números Aleatorios

#función que genera números aleatorios bajo el método congruencial multiplicativo
mm <- function(x0, a, m) {
  return((a * x0) %% m)
}

mmsim <- function(x0, a, m, nsim) {
  vec <- numeric(nsim + 1)
  vec[1] <- x0
  for (k in 1:nsim) {
    vec[k + 1] <- mm(vec[k], a, m)
  }
  return((vec[-1]) / m)
}

#función que genera números aleatorios bajo el método congruencial mixto
mi <- function(x0, a, m, c) {
  return((a * x0 + c) %% m)
}

mi_sim <- function(x0, a, m, c, nsim) {
  vec <- numeric(nsim + 1)
  vec[1] <- x0
  for (k in 1:nsim) {
    vec[k + 1] <- mi(vec[k], a, m, c)
  }
  return((vec[-1]) / m)
}

#función que genera números aleatorios bajo el método medios cuadrados

mc <- function(x0, k) {
  return(floor((x0^2 - floor((x0^2) / (10^(2 * k - k / 2))) * 10^(2 * k - k / 2)) / 10^(k / 2)))
}

mc_sim <- function(x0, k, n_sim) {
  vec <- numeric(n_sim + 1)
  vec[1] <- x0
  for (j in 1:n_sim) {
    vec[j + 1] <- mc(vec[j], k)
  }
  return(vec[-1] / (10^k))
}

#función que genera números aleatorios bajo el método de Lehmer

ml <- function(x0, n, c) {
  return(x0 * c - floor((x0 * c) / (10^n)) * 10^n - floor((x0 * c) / (10^n)))
}

ml_sim <- function(x0, n, c, nsim) {
  vec <- numeric(nsim + 1)
  vec[1] <- x0
  for (j in 1:nsim) {
    vec[j + 1] <- ml(x0 = vec[j], n, c)
  }
  return(vec[-1] / (10^n))
}

#-----------------------

# Integrales

generar_aleatorios_integral <- function(n, metodo) {
  if (metodo == "Congruencial Multiplicativo") {
    return(random_cong(a = 7^5, m = 2^31 - 1, x0 = as.numeric(now()), n = n))
  } else if (metodo == "Congruencial Mixto") {
    return(random_mixt(a = 7^5, c = 12345, m =2^31 - 1 , x0 = as.numeric(now()), n = n))
  } else {
    return(NULL)
  }
}

conv_matrix<-function(vector, cols=10){
  res<-rep(NA_real_,ceiling(length(vector)/cols)*cols) #ceiling redondea hacia arriba
  res[1: length(vector)]<- vector
  res<- as.data.frame(matrix(res, nrow=ceiling(length(vector)/cols), ncol=cols, byrow= TRUE))
  colnames(res)<-paste0("Cols", 1:cols) #paste0() concatena texto sin espacios.
  return(res)
}

conv_matrix(c(1,2,3,4,5,6,7,8), cols=5)

as.data.frame(conv_matrix(c(1,2,3,4,5,6,7,8), cols=5))

#-----------------------

#============= Anita
#FUNCIONES DE VARIABLES DISCRETAS
#Binomial
fun_binomial<- function(n,p){
  U<- runif(1)
  pk <- (1-p)^n
  k<- Fk<- 0
  repeat{
    Fk<- Fk+pk
    if(U<Fk){
      return (k)
      break
    }
    
    pk<- ((n-k)/(k+1))*(p/(1-p))*pk
    k<- k+1
  }
}
#Poisson
fun_poisson <- function(lambda) {
  U <- runif(1)
  k <- 0
  L <- exp(-lambda)
  p <- L
  
  while(U > p) {
    k <- k + 1
    L <- L * lambda / k
    p <- p + L
  }
  return(k)
}

# Método de la transformada inversa
fun_transformada_inversa <- function(pmf) {
  U <- runif(1)  # Generamos un número aleatorio uniforme
  S <- 0
  for (j in 1:length(pmf)) {
    S <- S + pmf[j]
    if (U < S) {
      return(j - 1)  # Devuelve el índice correspondiente a la variable aleatoria generada
    }
  }
}

# Método de composición
fun_composicion <- function(pmf1, pmf2, alpha) {
  U <- runif(1)  # Generamos un número aleatorio uniforme
  if (U < alpha) {
    return(sample(1:length(pmf1), 1, prob = pmf1))  # Generar de la primera distribución
  } else {
    return(sample(1:length(pmf2), 1, prob = pmf2))  # Generar de la segunda distribución
  }
}

# Técnica de aceptación y rechazo
fun_aceptacion_rechazo <- function(pmf, qmf) {
  c <- max(pmf / qmf)  # Encontrar la constante c
  repeat {
    Y <- sample(1:length(qmf), 1, prob = qmf, replace = TRUE)  # Generar de la distribución propuesta
    U <- runif(1)  # Generar un número aleatorio uniforme
    if (U < pmf[Y] / (c * qmf[Y])) {  # Aceptar con probabilidad proporcional
      return(Y - 1)  # Retornar el valor generado
    }
  }
}
 
#--------------
#============= Mateo 


# Función general de transformada inversa para continuas
t_inversa_continua <- function(F_inv, nval, ...) {
  U <- runif(nval)
  return(F_inv(U, ...))
}

# Funciones inversas para distribuciones comunes
inversa_exponencial <- function(u, lambda) {
  return(-log(1 - u) / lambda)
}

inversa_uniforme <- function(u, a, b) {
  return(a + (b - a) * u)
}

inversa_pareto <- function(u, alpha, xm) {
  return(xm / (1 - u)^(1/alpha))
}

inversa_weibull <- function(u, lambda, k) {
  return(lambda * (-log(1 - u))^(1/k))
}

# Funcion inversa personalizada
inversa_personalizada <- function(u, funcion_texto) {
  tryCatch({
    x_vals <- seq(0.001, 0.999, length.out = 1000)
    Fx <- sapply(x_vals, function(x) {
      eval(parse(text = funcion_texto))
    })
    
    if(any(diff(Fx) < 0)) {
      return(rep(NA, length(u)))
    }
    
    approx(Fx, x_vals, xout = u, rule = 2)$y
    
  }, error = function(e) {
    message("Error en la función personalizada: ", e$message)
    return(rep(NA, length(u)))
  })
}

# ======= MÉTODO DE RECHAZO - FUNCIONES AUXILIARES =======

# Función general del método de rechazo
metodo_rechazo <- function(nval, f_densidad, g_densidad, g_generadora, c_constante, ...) {
  res <- numeric(nval)
  i <- 0
  iteraciones <- 0
  
  while(i < nval) {
    Y <- g_generadora(1, ...)
    U <- runif(1)
    
    if(U <= f_densidad(Y) / (c_constante * g_densidad(Y, ...))) {
      i <- i + 1
      res[i] <- Y
    }
    iteraciones <- iteraciones + 1
  }
  
  return(list(datos = res, iteraciones = iteraciones, tasa_aceptacion = nval/iteraciones))
}

# Distribuciones predefinidas para el método de rechazo

# 1. Beta(2,4) - usando Uniforme(0,1) como propuesta
f_beta <- function(x) {
  20 * x * (1 - x)^3
}

g_unif <- function(n, a = 0, b = 1) {
  runif(n, a, b)
}

d_unif <- function(x, a = 0, b = 1) {
  dunif(x, a, b)
}

c_beta <- 135/64  # 2.109375

# 2. Normal estándar - usando Exponencial como propuesta para |Z|
f_normal_abs <- function(x) {
  (2/sqrt(2*pi)) * exp(-x^2/2)
}

g_exp <- function(n, rate = 1) {
  rexp(n, rate)
}

d_exp <- function(x, rate = 1) {
  dexp(x, rate)
}

c_normal <- sqrt(2*exp(1)/pi)  # 1.315489

# 3. Gamma(1.5, 1) - usando Exponencial como propuesta
f_gamma <- function(x) {
  (1/gamma(1.5)) * x^(0.5) * exp(-x)
}

c_gamma <- 1.257  # Valor precalculado

# 4. Distribución personalizada para método de rechazo
rechazo_personalizado <- function(nval, f_texto, g_texto, g_gen_texto, c_valor) {
  res <- numeric(nval)
  i <- 0
  iteraciones <- 0
  
  while(i < nval) {
    Y <- eval(parse(text = g_gen_texto))
    f_val <- eval(parse(text = f_texto))
    g_val <- eval(parse(text = g_texto))
    U <- runif(1)
    
    if(U <= f_val / (c_valor * g_val)) {
      i <- i + 1
      res[i] <- Y
    }
    iteraciones <- iteraciones + 1
  }
  
  return(list(datos = res, iteraciones = iteraciones, tasa_aceptacion = nval/iteraciones))
}

#---------------
# =========================
# FUNCIONES GENERADORAS
# =========================

grafico_bin <- function(x) {
  hist(x, breaks = 30, col = "aquamarine",
       main = "Distribución Binomial Simulada",
       xlab = "Valores", ylab = "Frecuencia")
}

grafico_pois <- function(x) {
  hist(x, breaks = 30, col = "orchid",
       main = "Distribución Poisson Simulada",
       xlab = "Valores", ylab = "Frecuencia")
}

grafico_inv <- function(x) {
  hist(x, breaks = 30, col = "lightblue",
       main = "Distribución Transformada Inversa",
       xlab = "Valores", ylab = "Frecuencia")
}

grafico_rechazo <- function(x) {
  hist(x, breaks = 30, col = "lightgreen",
       main = "Distribución Aceptación y Rechazo",
       xlab = "Valores", ylab = "Frecuencia")
}

grafico_composicion <- function(x) {
  hist(x, breaks = 30, col = "lightcoral",
       main = "Distribución Composición",
       xlab = "Valores", ylab = "Frecuencia")
}

#--------
grafico_bin <- function(sim_data, n, p) {
  df_sim <- data.frame(x = sim_data)
  
  # Distribución teórica
  x_teo <- 0:n
  y_teo <- dbinom(x_teo, size = n, prob = p)
  df_teo <- data.frame(x = x_teo, y = y_teo * length(sim_data)) # Escalar a frecuencia
  
  ggplot(df_sim, aes(x = x)) +
    geom_histogram(aes(y = after_stat(count)), bins = n + 1, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_point(data = df_teo, aes(x = x, y = y), color = "red", size = 3) +
    geom_segment(data = df_teo, aes(x = x, xend = x, y = 0, yend = y), color = "red", linewidth = 1, linetype = "dashed") +
    scale_x_continuous(breaks = x_teo) +
    labs(
      title = paste("Simulación Binomial (n=", n, ", p=", p, ")"),
      x = "Número de Éxitos (X)",
      y = "Frecuencia Observada"
    ) +
    theme_minimal()
}

# Función de gráfico base (Poisson)
grafico_pois <- function(sim_data, lambda) {
  df_sim <- data.frame(x = sim_data)
  
  # Rango para la teórica
  max_x <- max(sim_data) + 1
  x_teo <- 0:max_x
  y_teo <- dpois(x_teo, lambda = lambda)
  df_teo <- data.frame(x = x_teo, y = y_teo * length(sim_data)) # Escalar a frecuencia
  
  ggplot(df_sim, aes(x = x)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_point(data = df_teo, aes(x = x, y = y), color = "red", size = 3) +
    geom_segment(data = df_teo, aes(x = x, xend = x, y = 0, yend = y), color = "red", linewidth = 1, linetype = "dashed") +
    scale_x_continuous(breaks = x_teo) +
    labs(
      title = paste("Simulación Poisson (λ=", lambda, ")"),
      x = "Número de Eventos (X)",
      y = "Frecuencia Observada"
    ) +
    theme_minimal()
}

# Gráfico para métodos basados en PMF (Transformada Inversa, AyR, Composición)
grafico_pmf <- function(sim_data, pmf, titulo) {
  df_sim <- data.frame(x = sim_data)
  
  # Rango de valores (asumimos que los valores generados son 0, 1, 2, ...)
  unique_x <- sort(unique(sim_data))
  
  # Teórica (se asume que la pmf está indexada desde 0)
  x_teo <- 0:(length(pmf)-1)
  y_teo <- pmf
  df_teo <- data.frame(x = x_teo, y = y_teo * length(sim_data)) # Escalar a frecuencia
  
  ggplot(df_sim, aes(x = x)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_point(data = df_teo, aes(x = x, y = y), color = "red", size = 3) +
    geom_segment(data = df_teo, aes(x = x, xend = x, y = 0, yend = y), color = "red", linewidth = 1, linetype = "dashed") +
    scale_x_continuous(breaks = unique_x) +
    labs(
      title = paste("Simulación -", titulo),
      x = "Valor de la Variable Aleatoria (X)",
      y = "Frecuencia Observada"
    ) +
    theme_minimal()
}

grafico_inv <- function(sim_data) {
  pmf <- c(0.2, 0.3, 0.5) # Usamos el PMF de ejemplo
  grafico_pmf(sim_data, pmf, "Transformada Inversa")
}

grafico_rechazo <- function(sim_data) {
  pmf <- c(0.2, 0.3, 0.5) # Usamos el PMF de ejemplo
  grafico_pmf(sim_data, pmf, "Aceptación y Rechazo")
}

grafico_composicion <- function(sim_data) {
  # En el método de Composición, si no se definen las PMF en la UI, se usa un ejemplo de mezcla de dos.
  # Aquí asumimos que la distribución resultante se conoce para el plot.
  # Dado que no hay inputs para pmf1 y pmf2 en la UI, este gráfico es solo ilustrativo.
  # Usaremos una PMF de ejemplo simple para el plot:
  pmf_ejemplo <- c(0.1, 0.2, 0.3, 0.4)
  grafico_pmf(sim_data, pmf_ejemplo, "Composición")
}
#-------
#------------------
function(input, output, session) {
  
  #======= Emely
  
  output$parametros_ui <- renderUI({
    
    # Dependiendo del método seleccionado, se mostrarán diferentes controles de entrada
    if (input$metodo == "Congruencial Multiplicativo") {
      tagList(
        numericInput("semilla", "Ingrese un valor inicial:", min = 1, max = 500, value = 30),
        numericInput("divisor", "Ingrese un valor de m:", min = 5, max = 500, value = 37),
        numericInput("constante", "Ingrese un valor de a:", min = 10, max = 500, value = 123),
        numericInput("num", "Cantidad de números a generar:", min = 10, max = 200, value = 30)
      )
    }
    else if (input$metodo == "Congruencial Mixto") {
      tagList(
        numericInput("semilla", "Ingrese un valor inicial:", min = 1, max = 500, value = 30),
        numericInput("divisor", "Ingrese un valor de m:", min = 5, max = 500, value = 37),
        numericInput("constante", "Ingrese un valor de a:", min = 10, max = 500, value = 123),
        numericInput("c", "Ingrese un valor de c para el método mixto:", min = 1, max = 500, value = 200),
        numericInput("num", "Cantidad de números a generar:", min = 10, max = 200, value = 30)
      )
    }
    else if (input$metodo == "Cuadrados Medios") {
      tagList(
        numericInput("x0", "Ingrese un valor inicial x0:", min = 1, max = 10000, value = 9311),
        numericInput("k", "Ingrese el número de dígitos:", min = 1, max = 10, value = 4),
        numericInput("num", "Cantidad de números a generar:", min = 10, max = 200, value = 30)
      )
    }
    else if (input$metodo == "Lehmer") {
      tagList(
        numericInput("x0", "Ingrese un valor inicial x0:", min = 1, max = 10000, value = 4122),
        numericInput("n", "Ingrese el número de dígitos:", min = 1, max = 10, value = 4),
        numericInput("c", "Ingrese un valor de c:", min = 1, max = 100, value = 76),
        numericInput("num", "Cantidad de números a generar:", min = 10, max = 200, value = 30)
      )
    }
  })
  
  # Generación de aleatorios según el método seleccionado
  aleatorios <- eventReactive(input$mostrar, {
    if (input$metodo == "Congruencial Multiplicativo") {
      return(mmsim(x0 = input$semilla, a = input$constante, m = input$divisor, nsim = input$num))
    }
    else if (input$metodo == "Congruencial Mixto") {
      return(mi_sim(x0 = input$semilla, a = input$constante, m = input$divisor, c = input$c, nsim = input$num))
    }
    else if (input$metodo == "Cuadrados Medios") {
      return(mc_sim(x0 = input$x0, k = input$k, n_sim = input$num))
    }
    else if (input$metodo == "Lehmer") {
      return(ml_sim(x0 = input$x0, n = input$n, c = input$c, nsim = input$num))
    }
  })
  
  output$tabla <- DT::renderDataTable({
    res <- aleatorios()
    req(res)
    
    data.frame(Valor = res)
    
  }, options = list(
    pageLength = 10,
    lengthMenu = c(5, 10, 20, 50, 100),
    searching = TRUE,
    paging = TRUE,
    info = TRUE
  ))
  
  
  output$distPlot <- renderPlot({
    x <- aleatorios()
    bins <- seq(min(x), max(x), length.out = input$barras + 1)
    hist(x, breaks = bins, col = "lightcoral", border = 'white', 
         xlab = 'Numeros aleatorios generados', 
         main = paste('Histograma del Método: ', input$metodo),
         ylab = 'Frecuencia')  
  })
  
  # GENERAR EL REPORTE HTML
  output$reporte_aleatorios <- downloadHandler(
    
    filename = function() {
      nombre_base <- paste("Reporte_Numeros_Aleatorios_", gsub(" ", "_", input$metodo), ".html", sep = "")
      return(nombre_base)
    },
    
    content = function(file) {
      
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("El paquete 'rmarkdown' es necesario. Por favor, instálalo.")
      }
      
      req(input$mostrar > 0)
      temp_data <- aleatorios()
      
      params_reporte <- list(
        metodo = input$metodo,
        num_generar = input$num,
        barras = input$barras,
        
        semilla = if (input$metodo %in% c("Congruencial Multiplicativo", "Congruencial Mixto")) input$semilla else input$x0,
        
        constante_a = if (input$metodo %in% c("Congruencial Multiplicativo", "Congruencial Mixto")) input$constante else NULL,
        
        divisor_m = if (input$metodo %in% c("Congruencial Multiplicativo", "Congruencial Mixto")) input$divisor else NULL,
        
        c_val = if (input$metodo == "Congruencial Mixto" || input$metodo == "Lehmer") input$c else NULL,
        
        k_val = if (input$metodo == "Cuadrados Medios") input$k else NULL,
        
        n_val = if (input$metodo == "Lehmer") input$n else NULL,
        
        # Datos Generados
        datos_generados = temp_data
      )
      
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("El paquete 'rmarkdown' es necesario. Por favor, instálalo.")
      }
      
      rmarkdown::render(
        input = "reporte_aleatorios.qmd", # Tu archivo de reporte
        output_file = file,
        output_format = "html_document",  # Especificar el formato HTML
        params = params_reporte,
        envir = new.env(parent = globalenv())
        )
    }
  )
  
  
  
  #------------------
  #Integrales
  
  generar_aleatorios_integral <- function(n, metodo) {
    x0 <- as.numeric(format(Sys.time(), "%OS6")) %% 10000
    
    if (metodo == "Congruencial Multiplicativo") {
      
      return(mmsim(
        x0 = x0,
        a = 7^5,
        m = 2^31 - 1,
        nsim = n
      ))
      
    } else if (metodo == "Congruencial Mixto") {
      
      return(mi_sim(
        x0 = x0,
        a = 7^5,
        m = 2^31 - 1,
        c = 12345,
        nsim = n
      ))
      
    } else if (metodo == "Cuadrados Medios") {
      
      x0_cm <- as.numeric(substr(x0, 1, 4))
      return(mc_sim(
        x0 = x0_cm,
        k = 4,
        n_sim = n
      ))
      
    } else if (metodo == "Lehmer") {
      
      return(ml_sim(
        x0 = x0 %% 10000,
        n = 4,
        c = 16807,
        nsim = n
      ))
    }
    
    return(NULL)
  }
  
  datos_integral <- reactive({
    req(input$calcular)
    
    f <- function(x) eval(parse(text = input$funcion))
    
    a_inf <- ifelse(input$inf_inf, -Inf, input$lim_inf)
    b_inf <- ifelse(input$sup_inf,  Inf,  input$lim_sup)
    
    if (is.infinite(a_inf) && is.infinite(b_inf)) {
      plot_a <- -5
      plot_b <-  5
    } else if (is.infinite(a_inf) && !is.infinite(b_inf)) {
      plot_b <- b_inf
      plot_a <- b_inf - 5
    } else if (!is.infinite(a_inf) && is.infinite(b_inf)) {
      plot_a <- a_inf
      plot_b <- a_inf + 5
    } else {
      plot_a <- a_inf
      plot_b <- b_inf
    }
    
    x_vals <- seq(plot_a, plot_b, length.out = 400)
    y_vals <- sapply(x_vals, f)
    
    data.frame(
      x      = x_vals,
      y      = y_vals,
      a_inf  = a_inf,
      b_inf  = b_inf
    )
  })
  
  # Calcular área bajo la curva (regla trapezoidal)
  area_aprox <- eventReactive(input$calcular, {
    f <- function(x) eval(parse(text = input$funcion))
    a <- ifelse(input$inf_inf, -Inf, input$lim_inf)  
    b <- ifelse(input$sup_inf, Inf, input$lim_sup)
    
    # Caso de la integral de [a, ∞)
    if (a > -Inf && b == Inf) {
      h <- function(y) {
        if (y <= 0 || y >= 1) return(0)
        x <- (1 / y) - 1 + a
        f(x) / y^2
      }
      y_vals <- seq(1e-6, 1 - 1e-6, length.out = 2000)
      dy <- y_vals[2] - y_vals[1]
      area <- (dy/2)*(h(y_vals[1]) +
                        2 * sum(sapply(y_vals[2:(length(y_vals)-1)], h)) +
                        h(y_vals[length(y_vals)]))
      return(area)
    }
    # Caso de la integral de (-∞, b]
    if (a == -Inf && b < Inf) {
      h <- function(y) {
        if (y <= 0 || y >= 1) return(0)
        x <- b - ((1 / y) - 1)
        f(x) / y^2
      }
      y_vals <- seq(1e-6, 1 - 1e-6, length.out = 2000)
      dy <- y_vals[2] - y_vals[1]
      area <- (dy/2)*(h(y_vals[1]) +
                        2 * sum(sapply(y_vals[2:(length(y_vals)-1)], h)) +
                        h(y_vals[length(y_vals)]))
      return(area)
    }
    # Caso de la integral de (-∞, ∞)
    if (a == -Inf && b == Inf) {
      h_pos <- function(y) {
        if (y <= 0 || y >= 1) return(0)
        x <- (1 / y) - 1
        f(x) / y^2
      }
      h_neg <- function(y) {
        if (y <= 0 || y >= 1) return(0)
        x <- -((1 / y) - 1)
        f(x) / y^2
      }
      y_vals <- seq(1e-6, 1 - 1e-6, length.out = 2000)
      dy <- y_vals[2] - y_vals[1]
      
      area_pos <- (dy/2)*(h_pos(y_vals[1]) +
                            2 * sum(sapply(y_vals[2:(length(y_vals)-1)], h_pos)) +
                            h_pos(y_vals[length(y_vals)]))
      area_neg <- (dy/2)*(h_neg(y_vals[1]) +
                            2 * sum(sapply(y_vals[2:(length(y_vals)-1)], h_neg)) +
                            h_neg(y_vals[length(y_vals)]))
      return(area_pos + area_neg)
    }
    
    # Caso de la integral de [a, b]
    datos <- datos_integral()
    a_f <- ifelse(input$inf_inf, -Inf, input$lim_inf)
    b_f <- ifelse(input$sup_inf,  Inf,  input$lim_sup)
    
    delta_x <- (b_f - a_f) / (nrow(datos) - 1)
    return((delta_x/2) * (datos$y[1] +
                            2 * sum(datos$y[2:(nrow(datos)-1)]) +
                            datos$y[nrow(datos)]))
  })
  
  # Gráfica de la función
  output$graf_fun01 <- renderPlot({
    req(input$calcular)
    f <- function(x) eval(parse(text = input$funcion))
    
    datos <- datos_integral()
    area  <- area_aprox()
    
    a <- datos$a_inf[1]
    b <- datos$b_inf[1]
    
    if (any(is.na(datos$y))) {
      plot.new()
      title("Error: la función no se pudo evaluar en algún punto.")
      return()
    }
    
    g <- ggplot(datos, aes(x = x, y = y)) +
      geom_line(color = "blue", linewidth = 1) +
      geom_area(fill = "lightblue", alpha = 0.4) +
      theme_minimal() +
      labs(
        title = paste("Área bajo la curva de f(x) =", input$funcion),
        subtitle = ifelse(is.infinite(a) && is.infinite(b),
                          "Intervalo: (-∞, ∞)",
                          ifelse(is.infinite(a),
                                 paste("Intervalo: (-∞, ", b, "]"),
                                 ifelse(is.infinite(b),
                                        paste("Intervalo: [", a, ", ∞)"),
                                        paste("Intervalo: [", a, ",", b, "]")
                                 ))),
        caption = paste("Área aproximada:", round(area, 6)),
        x = "x", y = "f(x)"
      )
    if (!is.infinite(a)) {
      g <- g + geom_vline(xintercept = a, linetype = "dashed", color = "red")
    }
    if (!is.infinite(b)) {
      g <- g + geom_vline(xintercept = b, linetype = "dashed", color = "blue")
    }
    
  
    g
  })
  
  output$nota_metodo <- renderUI({
    req(input$metodo) 
    
    if (input$metodo == "Cuadrados Medios") {
      return(
        div(
          style = "
          background-color:#fff3cd;
          border-left: 6px solid #ffcc00;
          padding:12px;
          margin-top:10px;
          border-radius:8px;
          color:#856404;
          font-family:'Poppins';
        ",
          strong("Nota sobre el método de Cuadrados Medios:"),
          p("Este método no genera números pseudoaleatorios uniformes de buena calidad. 
           Su ciclo es extremadamente corto y tiende a degenerarse rápidamente, por lo que 
           NO es recomendable para aproximaciones por Monte Carlo.")
        )
      )
    }
    
    if (input$metodo == "Lehmer") {
      return(
        div(
          style = "
          background-color:#fde2e4;
          border-left: 6px solid #ff6b6b;
          padding:12px;
          margin-top:10px;
          border-radius:8px;
          color:#7a0000;
          font-family:'Poppins';
        ",
          strong("Nota sobre el método de Lehmer:"),
          p("En este método, la secuencia generada puede presentar comportamientos particulares 
     que afectan la calidad de los valores pseudoaleatorios. 
     Debido a ello, las aproximaciones por Monte Carlo podrían no ser tan precisas al 
     compararse con otros métodos de generación.")
        )
      )
    }
    
    return(NULL)
  })
  
  output$valores_mc <- renderText({
    req(input$calcular)
    
    datos_mc <- valores_mc_data()
    
    aprox_final <- attr(datos_mc, "aprox_final")
    teorico     <- attr(datos_mc, "teorico")
    
    paste0(
      "Aproximación Monte Carlo: ", round(aprox_final, 6), "\n",
      "Valor teórico: ",
      ifelse(is.na(teorico), "No disponible (la integral puede divergir o no pudo ser calculada)",
             round(teorico, 6))
    )
  })
    
  
  valores_mc_data <- eventReactive(input$calcular, {
    req(input$calcular)
  
    f <- function(x) eval(parse(text = input$funcion))
    a <- ifelse(input$inf_inf, -Inf, input$lim_inf)
    b <- ifelse(input$sup_inf,  Inf,  input$lim_sup)
    
    
    secuencia <- seq(100, 10000, by = 250)
  
    # Cálculo de aproximaciones Monte Carlo
    # -------------------------------------------
    aprox <- sapply(secuencia, function(k) {
      u <- generar_aleatorios_integral(n = k, metodo = input$metodo)
      
      # Casos según tipo de intervalo
      if (a == 0 && b == 1) {
        # [0,1]
        return(mean(sapply(u, f)))
        
      } else if (a > -Inf && b == Inf) {
        # [a, +∞)
        h <- function(y) {
          y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
          x <- (1 / y) - 1 + a
          f(x) / y^2
        }
        return(mean(sapply(u, h)))
        
      } else if (a == -Inf && b < Inf) {
        # (-∞, b]
        h <- function(y) {
          y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
          x <- b - ((1 / y) - 1)
          f(x) / y^2
        }
        return(mean(sapply(u, h)))
        
      } else if (a == -Inf && b == Inf) {
        # (-∞, +∞)
        h_pos <- function(y) {
          y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
          x <- (1 / y) - 1
          f(x) / y^2
        }
        h_neg <- function(y) {
          y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
          x <- -((1 / y) - 1)
          f(x) / y^2
        }
        return(mean(sapply(u, h_pos)) + mean(sapply(u, h_neg)))
        
      } else {
        # [a, b] finito
        x_vals <- a + (b - a) * u
        return(mean((b - a) * sapply(x_vals, f)))
      }
    })
    
    # Cálculo del valor teórico (si existe)
    
    teorico <- tryCatch(
      integrate(f, lower = a, upper = b)$value,
      error = function(e) NA
    )
    
    # Construimos tabla larga para el gráfico
    dt <- data.table(
      Aleatorios = rep(secuencia, times = 2),
      Etiqueta   = rep(c("Aproximación", "Teórico"), each = length(secuencia)),
      Valor      = c(aprox, rep(teorico, length(secuencia)))
    )
    
    attr(dt, "aprox_final") <- tail(aprox, 1)
    attr(dt, "teorico")     <- teorico
    
    dt
  })
  
  # Aproximación por método de Monte Carlo
  output$graf_aprox01 <- renderPlot({
    req(input$calcular) 
    datos_mc <- valores_mc_data()
    
    ggplot(datos_mc, aes(x = Aleatorios, y = Valor, colour = Etiqueta)) +
      geom_line(linewidth = 1) +
      theme_minimal() +
      labs(
        title = "Aproximación por Monte Carlo vs. valor teórico",
        x = "Cantidad de números aleatorios",
        y = "Valor"
      )
  })
      
  observe({
    if (input$inf_inf) {
      disable("lim_inf")
      updateNumericInput(session, "lim_inf", value = NA)
    } else {
      enable("lim_inf")
      if (is.na(input$lim_inf)) {
        updateNumericInput(session, "lim_inf", value = 0)
      }
    }
  })
  
  observe({
    if (input$sup_inf) {
      disable("lim_sup")
      updateNumericInput(session, "lim_sup", value = NA)
    } else {
      enable("lim_sup")
      if (is.na(input$lim_sup)) {
        updateNumericInput(session, "lim_sup", value = 1)
      }
    }
  })
  
  # GENERAR EL REPORTE DE INTEGRALES
  output$reporte_integrales <- downloadHandler(
    
    filename = function() {
      paste("Reporte_Integral_", gsub(" ", "_", input$metodo), ".html", sep = "")
    },
    
    content = function(file) {
      
      req(input$calcular)
    
      data_mc <- valores_mc_data() 
      aprox_final <- attr(data_mc, "aprox_final")
      teorico <- attr(data_mc, "teorico")
      
      area_trapezoidal <- area_aprox()
      
      a_lim <- ifelse(input$inf_inf, -Inf, input$lim_inf)
      b_lim <- ifelse(input$sup_inf, Inf, input$lim_sup)
      
      params_reporte_int <- list(
        funcion = input$funcion,
        lim_inf = a_lim,
        lim_sup = b_lim,
        metodo_aleatorios = input$metodo,
        
        # Resultados
        aprox_mc = aprox_final,
        aprox_trapezoidal = area_trapezoidal,
        valor_teorico = teorico,
        
        # Datos para el Gráfico
        data_mc_plot = data_mc
      )
      
      rmarkdown::render("reporte_integrales.qmd",
                        output_file = file,
                        output_format = "html_document", 
                        params = params_reporte_int,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  #-----------------------------
  
  #======= Anita
  
  #Pestaña Variables Discretas
  # Pestaña Variables Discretas - Versión corregida
  
  # Variable para controlar el estado de la distribución seleccionada
  dist_seleccionada <- reactiveVal(NULL)
  
  # Control de los parámetros mostrados en función de la distribución seleccionada
  observeEvent(input$siguiente, {
    dist_seleccionada(input$distribucion)
  })
  
  output$parametrosUI <- renderUI({
    req(dist_seleccionada())
    
    if(dist_seleccionada() == "Binomial") {
      tagList(
        numericInput("n_bin", "Número de ensayos (n):", value = 10, min = 1),
        numericInput("p_bin", "Probabilidad de éxito (p):", value = 0.5, min = 0, max = 1, step = 0.01),
        numericInput("num_bin", "Número de simulaciones:", value = 1000, min = 100),
        actionButton("calc_bin", "Calcular", class = "btn-success")
      )
    } else if(dist_seleccionada() == "Poisson") {
      tagList(
        numericInput("lambda", "Valor de lambda (λ):", value = 4, min = 0.1, step = 0.1),
        numericInput("num_pois", "Número de simulaciones:", value = 1000, min = 100),
        actionButton("calc_pois", "Calcular", class = "btn-success")
      )
    } else if(dist_seleccionada() == "Transformada Inversa") {
      tagList(
        numericInput("num_inv", "Número de simulaciones:", value = 1000, min = 100),
        actionButton("calc_inv", "Calcular", class = "btn-success")
      )
    } else if(dist_seleccionada() == "Aceptación y Rechazo") {
      tagList(
        numericInput("num_rechazo", "Número de simulaciones:", value = 1000, min = 100),
        actionButton("calc_rechazo", "Calcular", class = "btn-success")
      )
    } else if(dist_seleccionada() == "Composición") {
      tagList(
        numericInput("num_composicion", "Número de simulaciones:", value = 1000, min = 100),
        numericInput("alpha", "Probabilidad de la primera distribución (α):", value = 0.5, min = 0, max = 1, step = 0.01),
        actionButton("calc_composicion", "Calcular", class = "btn-success")
      )
    }
  })
  sim_data_discretas <- reactive({
    req(input$distribucion)
    
    # Define PMFs de ejemplo para TI, AyR, Composición (deberían ser inputs reales en la versión final)
    # Por simplicidad, usamos los valores definidos en las funciones
    pmf_ti <- c(0.2, 0.3, 0.5)
    pmf_ayr <- c(0.2, 0.3, 0.5)
    qmf_ayr <- c(0.4, 0.4, 0.2) # Debe cumplir sum(qmf) = 1
    pmf_comp1 <- c(0.1, 0.9)
    pmf_comp2 <- c(0.5, 0.5)
    
    if (input$distribucion == "Binomial") {
      req(input$calc_bin > 0, input$n_bin, input$p_bin, input$num_bin)
      x <- sim_bin() # Datos simulados
      n <- input$n_bin
      p <- input$p_bin
      
      # Estadísticas Teóricas
      teo_media <- n * p
      teo_varianza <- n * p * (1 - p)
      
      # Plot
      plot_obj <- grafico_bin(x, n, p)
      
      # Parámetros para el reporte
      params_reporte <- paste0("n = ", n, ", p = ", p)
      
    } else if (input$distribucion == "Poisson") {
      req(input$calc_pois > 0, input$lambda, input$num_pois)
      x <- sim_pois()
      lambda <- input$lambda
      
      # Estadísticas Teóricas
      teo_media <- lambda
      teo_varianza <- lambda
      
      # Plot
      plot_obj <- grafico_pois(x, lambda)
      
      # Parámetros para el reporte
      params_reporte <- paste0("λ = ", lambda)
      
    } else if (input$distribucion == "Transformada Inversa") {
      req(input$calc_inv > 0, input$num_inv)
      x <- sim_inv()
      pmf <- pmf_ti
      
      # Estadísticas Teóricas (Media = sum(x * p(x)), Varianza = sum(x^2 * p(x)) - media^2)
      valores_x <- 0:(length(pmf) - 1)
      teo_media <- sum(valores_x * pmf)
      teo_varianza <- sum(valores_x^2 * pmf) - teo_media^2
      
      # Plot
      plot_obj <- grafico_pmf(x, pmf, "Transformada Inversa")
      
      # Parámetros para el reporte
      params_reporte <- paste0("PMF = c(", paste(pmf, collapse = ", "), ")")
      
    } else if (input$distribucion == "Aceptación y Rechazo") {
      req(input$calc_rechazo > 0, input$num_rechazo)
      x <- sim_rechazo()
      pmf <- pmf_ayr
      qmf <- qmf_ayr
      
      # Estadísticas Teóricas (asumimos que la PMF objetivo es la pmf_ayr)
      valores_x <- 0:(length(pmf) - 1)
      teo_media <- sum(valores_x * pmf)
      teo_varianza <- sum(valores_x^2 * pmf) - teo_media^2
      
      # Plot
      plot_obj <- grafico_pmf(x, pmf, "Aceptación y Rechazo")
      
      # Parámetros para el reporte
      params_reporte <- paste0("PMF (Objetivo) = c(", paste(pmf, collapse = ", "), ")",
                               "; QMF (Propuesta) = c(", paste(qmf, collapse = ", "), ")")
      
    } else if (input$distribucion == "Composición") {
      req(input$calc_composicion > 0, input$num_composicion, input$alpha)
      # Nota: La función fun_composicion usa sample(), pero aquí la simularemos con los datos generados
      # x <- sim_composicion() # Esto requiere input$pmf1 y input$pmf2, que no están en la UI
      
      # Para evitar el error por inputs faltantes, usamos una simulación directa aquí si es necesario, 
      # o asumimos que sim_composicion() funciona con los datos de ejemplo (o si los inputs se agregaron)
      # Usaremos la pmf resultante de la mezcla teórica para los stats
      alpha <- input$alpha
      pmf_final <- alpha * pmf_comp1 + (1 - alpha) * pmf_comp2
      x <- replicate(input$num_composicion, fun_composicion(pmf_comp1, pmf_comp2, alpha))
      
      # Estadísticas Teóricas (usando la PMF resultante)
      valores_x <- 1:length(pmf_final) # Asumimos índice 1, 2, ...
      # Se necesita ajustar si fun_composicion devuelve 1, 2, ... o 0, 1, ...
      # sample(1:length(pmf), 1, prob = pmf) devuelve 1, 2, ...
      
      teo_media <- sum(valores_x * pmf_final)
      teo_varianza <- sum(valores_x^2 * pmf_final) - teo_media^2
      
      # Plot
      plot_obj <- grafico_pmf(x, pmf_final, "Composición")
      
      # Parámetros para el reporte
      params_reporte <- paste0("PMF 1 = c(", paste(pmf_comp1, collapse = ", "), ")", 
                               "; PMF 2 = c(", paste(pmf_comp2, collapse = ", "), ")",
                               "; α = ", alpha)
    } else {
      return(NULL)
    }
    
    # Estadísticas Simuladas
    sim_media <- mean(x)
    sim_varianza <- var(x)
    
    list(
      sim_data = x,
      sim_media = sim_media,
      sim_varianza = sim_varianza,
      teo_media = teo_media,
      teo_varianza = teo_varianza,
      metodo = input$distribucion,
      simulaciones = length(x),
      parametros = params_reporte,
      plot = plot_obj
    )
  })
  
  
  # Funciones para manejar las distribuciones discretas
  sim_bin <- eventReactive(input$calc_bin, {
    req(input$n_bin, input$p_bin, input$num_bin)
    replicate(input$num_bin, fun_binomial(input$n_bin, input$p_bin))
  })
  
  output$histBin <- renderPlot({
    req(sim_bin())
    grafico_bin(sim_bin())
  })
  
  output$statsBin <- renderPrint({
    req(sim_bin())
    x <- sim_bin()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
  sim_pois <- eventReactive(input$calc_pois, {
    req(input$lambda, input$num_pois)
    replicate(input$num_pois, fun_poisson(input$lambda))
  })
  
  output$histPois <- renderPlot({
    req(sim_pois())
    grafico_pois(sim_pois())
  })
  
  output$statsPois <- renderPrint({
    req(sim_pois())
    x <- sim_pois()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
  # Resultados de Transformada Inversa
  # Corrección para el método de Transformada Inversa
  sim_inv <- eventReactive(input$calc_inv, {
    req(input$num_inv)  # Asegúrate de que el número de simulaciones no sea vacío
    
    # Reemplazar "input$pmf" por una distribución válida, aquí puedes ajustar según la distribución seleccionada
    pmf <- c(0.2, 0.3, 0.5)  # Ejemplo de PMF
    replicate(input$num_inv, fun_transformada_inversa(pmf))
  })
  
  output$histInv <- renderPlot({
    req(sim_inv())
    grafico_inv(sim_inv())
  })
  
  output$statsInv <- renderPrint({
    req(sim_inv())
    x <- sim_inv()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
  
  # Resultados de Aceptación y Rechazo
  # Corrección para el método de Aceptación y Rechazo
  sim_rechazo <- eventReactive(input$calc_rechazo, {
    req(input$num_rechazo)  # Asegúrate de que el número de simulaciones no sea vacío
    
    # Distribuciones PMF y QMF ejemplo
    pmf <- c(0.2, 0.3, 0.5)  # Ejemplo de PMF
    qmf <- c(0.4, 0.4, 0.2)  # Ejemplo de QMF
    replicate(input$num_rechazo, fun_aceptacion_rechazo(pmf, qmf))
  })
  
  output$histRechazo <- renderPlot({
    req(sim_rechazo())
    grafico_rechazo(sim_rechazo())
  })
  
  output$statsRechazo <- renderPrint({
    req(sim_rechazo())
    x <- sim_rechazo()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
  
  # Resultados de Composición
  sim_composicion <- eventReactive(input$calc_composicion, {
    replicate(input$num_composicion, fun_composicion(input$pmf1, input$pmf2, input$alpha))
  })
  
  output$histComposicion <- renderPlot({
    req(sim_composicion())
    grafico_composicion(sim_composicion())
  })
  
  output$statsComposicion <- renderPrint({
    req(sim_composicion())
    x <- sim_composicion()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
  
   # ================ Mateo
   # ===== TRANSFORMADA INVERSA =====
  
  # UI dinámica para transformada inversa
  output$parametros_continuos_ui <- renderUI({
    req(input$dist_continua)
    
    if (input$dist_continua == "exponencial") {
      tagList(
        numericInput("param1_cont", "λ (tasa):", value = 1, min = 0.1, step = 0.1)
      )
    }
    else if (input$dist_continua == "uniforme") {
      tagList(
        numericInput("param1_cont", "a (límite inferior):", value = 0),
        numericInput("param2_cont", "b (límite superior):", value = 1)
      )
    }
    else if (input$dist_continua == "pareto") {
      tagList(
        numericInput("param1_cont", "α (forma):", value = 2, min = 0.1, step = 0.1),
        numericInput("param2_cont", "xm (escala):", value = 1, min = 0.1, step = 0.1)
      )
    }
    else if (input$dist_continua == "weibull") {
      tagList(
        numericInput("param1_cont", "λ (escala):", value = 1, min = 0.1, step = 0.1),
        numericInput("param2_cont", "k (forma):", value = 1, min = 0.1, step = 0.1)
      )
    }
    else if (input$dist_continua == "personalizada") {
      tagList(
        textInput("funcion_personalizada_cont", "Función de distribución F(x):",
                  placeholder = "Ej: x^2 para 0<x<1"),
        helpText("Ingrese F(x) en términos de x. Ejemplos:"),
        helpText("• x^2 (F(x)=x²)"),
        helpText("• (x^2 + x)/2"),
        helpText("• 1 - exp(-x) (exponencial)")
      )
    }
  })
  
  # Datos para transformada inversa
  datos_continuos <- eventReactive(input$generar_continuos, {
    req(input$nval_continua, input$dist_continua)
    n <- as.numeric(input$nval_continua)
    
    if (input$dist_continua == "exponencial") {
      req(input$param1_cont)
      lambda <- as.numeric(input$param1_cont)
      datos_gen <- t_inversa_continua(inversa_exponencial, n, lambda = lambda)
      nombre <- paste("Exponencial(λ =", lambda, ")")
      color <- "lightblue"
    }
    else if (input$dist_continua == "uniforme") {
      req(input$param1_cont, input$param2_cont)
      a <- as.numeric(input$param1_cont)
      b <- as.numeric(input$param2_cont)
      datos_gen <- t_inversa_continua(inversa_uniforme, n, a = a, b = b)
      nombre <- paste("Uniforme(a =", a, ", b =", b, ")")
      color <- "lightgreen"
    }
    else if (input$dist_continua == "pareto") {
      req(input$param1_cont, input$param2_cont)
      alpha <- as.numeric(input$param1_cont)
      xm <- as.numeric(input$param2_cont)
      datos_gen <- t_inversa_continua(inversa_pareto, n, alpha = alpha, xm = xm)
      nombre <- paste("Pareto(α =", alpha, ", xm =", xm, ")")
      color <- "lightcoral"
    }
    else if (input$dist_continua == "weibull") {
      req(input$param1_cont, input$param2_cont)
      lambda <- as.numeric(input$param1_cont)
      k <- as.numeric(input$param2_cont)
      datos_gen <- t_inversa_continua(inversa_weibull, n, lambda = lambda, k = k)
      nombre <- paste("Weibull(λ =", lambda, ", k =", k, ")")
      color <- "lightgoldenrod"
    }
    else if (input$dist_continua == "personalizada") {
      req(input$funcion_personalizada_cont)
      datos_gen <- inversa_personalizada(runif(n), input$funcion_personalizada_cont)
      nombre <- "Distribución Personalizada"
      color <- "lightsteelblue"
    }
    
    list(datos = datos_gen, nombre = nombre, color = color)
  })
  
  # Histograma para transformada inversa
  output$hist_continuo <- renderPlot({
    datos <- datos_continuos()
    req(datos$datos)
    
    if (all(is.na(datos$datos))) {
      plot(1, 1, type = "n", xlab = "", ylab = "", 
           main = "Error en la función personalizada")
      text(1, 1, "Revise la sintaxis de la función", col = "red")
    } else {
      ggplot(data.frame(x = datos$datos), aes(x = x)) +
        geom_histogram(aes(y = ..density..), bins = 30, 
                       fill = datos$color, color = "black", alpha = 0.7) +
        geom_density(color = "red", size = 1) +
        labs(title = paste("Distribución de", datos$nombre),
             x = "Valores generados", y = "Densidad") +
        theme_minimal()
    }
  })
  
  # Estadísticas para transformada inversa
  output$estadisticas_continuas <- renderTable({
    datos <- datos_continuos()$datos
    req(datos)
    
    if (!all(is.na(datos))) {
      stats <- data.frame(
        Estadística = c("Media", "Desviación Estándar", "Mínimo", "Máximo", "Cantidad"),
        Valor = c(
          round(mean(datos), 4),
          round(sd(datos), 4),
          round(min(datos), 4),
          round(max(datos), 4),
          length(datos)
        )
      )
      stats
    }
  })
  
  # Tabla para transformada inversa
  output$tabla_continuos <- renderDataTable({
    datos <- datos_continuos()$datos
    req(datos)
    
    if (!all(is.na(datos))) {
      data.frame(Valor = datos)
    }
  }, options = list(pageLength = 10))
  
  # ===== MÉTODO DE RECHAZO =====
  
  # UI dinámica para método de rechazo
  output$parametros_rechazo_ui <- renderUI({
    req(input$dist_rechazo)
    
    if (input$dist_rechazo == "personalizada_rechazo") {
      tagList(
        textInput("f_densidad", "Función de densidad f(x):", 
                  placeholder = "20*x*(1-x)^3"),
        textInput("g_densidad", "Densidad de la distribución propuesta g(x):",
                  placeholder = "1"),
        textInput("g_generadora", "Generadora de g (ej: runif(1)):",
                  placeholder = "runif(1)"),
        numericInput("c_rechazo", "Constante c:", value = 2.0, min = 1, step = 0.1),
        helpText("Ejemplo para Beta(2,4):"),
        helpText("f(x): 20*x*(1-x)^3"),
        helpText("g(x): 1"),
        helpText("Generadora g: runif(1)"),
        helpText("c: 2.109375")
      )
    } else {
      helpText("Parámetros predefinidos para la distribución seleccionada")
    }
  })
  
  # Datos para método de rechazo
  datos_rechazo <- eventReactive(input$generar_rechazo, {
    req(input$dist_rechazo, input$nval_rechazo)
    
    n <- as.numeric(input$nval_rechazo)
    
    if (input$dist_rechazo == "beta_24") {
      resultado <- metodo_rechazo(n, f_beta, d_unif, g_unif, c_beta)
      nombre <- "Beta(2,4)"
      color <- "lightcoral"
    }
    else if (input$dist_rechazo == "normal") {
      resultado_abs <- metodo_rechazo(n, f_normal_abs, d_exp, g_exp, c_normal)
      datos_con_signo <- resultado_abs$datos * sample(c(-1, 1), n, replace = TRUE)
      
      resultado <- list(
        datos = datos_con_signo,
        iteraciones = resultado_abs$iteraciones,
        tasa_aceptacion = resultado_abs$tasa_aceptacion
      )
      nombre <- "Normal Estándar"
      color <- "lightblue"
    }
    else if (input$dist_rechazo == "gamma_15") {
      resultado <- metodo_rechazo(n, f_gamma, d_exp, g_exp, c_gamma)
      nombre <- "Gamma(1.5, 1)"
      color <- "lightgreen"
    }
    else if (input$dist_rechazo == "personalizada_rechazo") {
      req(input$f_densidad, input$g_densidad, input$g_generadora, input$c_rechazo)
      
      resultado <- rechazo_personalizado(
        n, 
        input$f_densidad, 
        input$g_densidad, 
        input$g_generadora, 
        as.numeric(input$c_rechazo)
      )
      nombre <- "Distribución Personalizada"
      color <- "lightgoldenrod"
    }
    
    list(
      datos = resultado$datos,
      nombre = nombre,
      color = color,
      iteraciones = resultado$iteraciones,
      tasa_aceptacion = resultado$tasa_aceptacion
    )
  })
  
  # Histograma para método de rechazo
  output$hist_rechazo <- renderPlot({
    datos <- datos_rechazo()
    req(datos$datos)
    
    ggplot(data.frame(x = datos$datos), aes(x = x)) +
      geom_histogram(aes(y = ..density..), bins = 30, 
                     fill = datos$color, color = "black", alpha = 0.7) +
      geom_density(color = "red", size = 1) +
      labs(title = paste("Distribución de", datos$nombre, "- Método de Rechazo"),
           x = "Valores generados", y = "Densidad") +
      theme_minimal()
  })
  
  # Estadísticas para método de rechazo
  output$estadisticas_rechazo <- renderTable({
    datos <- datos_rechazo()
    req(datos$datos)
    
    stats <- data.frame(
      Estadística = c("Media", "Desviación Estándar", "Mínimo", "Máximo", 
                      "Iteraciones totales", "Tasa de aceptación", "Cantidad"),
      Valor = c(
        round(mean(datos$datos), 4),
        round(sd(datos$datos), 4),
        round(min(datos$datos), 4),
        round(max(datos$datos), 4),
        datos$iteraciones,
        paste0(round(datos$tasa_aceptacion * 100, 2), "%"),
        length(datos$datos)
      )
    )
    stats
  })
  
  # Tabla para método de rechazo
  output$tabla_rechazo <- renderDataTable({
    datos <- datos_rechazo()$datos
    req(datos)
    
    data.frame(Valor = datos)
  }, options = list(pageLength = 10)) 
  
  
}

