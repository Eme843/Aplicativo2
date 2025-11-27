library(shiny)
library(kableExtra)
library(data.table)
library(ggplot2)
library(lubridate)
library(tidyverse)


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
  
  output$tabla <- renderUI({
    res <- aleatorios()
    req(res)
    
    
    num_cols <- ceiling(sqrt(length(res)))  
    n_rows <- ceiling(length(res) / num_cols)
    
    res <- c(res, rep(NA, n_rows * num_cols - length(res)))
    
    matrix_data <- matrix(res, nrow = n_rows, ncol = num_cols, byrow = TRUE)
    colnames(matrix_data) <- paste("Col", 1:num_cols)  
    
    table_html <- kbl(matrix_data, escape = FALSE) %>%
      kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "bordered"), 
                    position = "left", font_size = 14) %>%
      row_spec(0, background = "#40E0D0", color = "white", bold = TRUE)
    
    # Retornar el HTML como un objeto renderizado en UI
    HTML(table_html)
  })
  
  output$distPlot <- renderPlot({
    x <- aleatorios()
    bins <- seq(min(x), max(x), length.out = input$barras + 1)
    hist(x, breaks = bins, col = "#00CD00", border = 'white', 
         xlab = 'Numeros aleatorios generados', 
         main = paste('Histograma del Método: ', input$metodo),
         ylab = 'Frecuencia')  
  })
  
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
  
  # Funciones para manejar las distribuciones discretas
  sim_bin <- eventReactive(input$calc_bin, {
    req(input$n_bin, input$p_bin, input$num_bin)
    replicate(input$num_bin, fun_binomial(input$n_bin, input$p_bin))
  })
  
  output$histBin <- renderPlot({
    req(sim_bin())
    hist(sim_bin(), breaks = 30, col = "aquamarine", 
         main = "Distribución Binomial Simulada",
         xlab = "Valores", ylab = "Frecuencia")
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
    hist(sim_pois(), breaks = 30, col = "orchid", 
         main = "Distribución Poisson Simulada",
         xlab = "Valores", ylab = "Frecuencia")
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
    req(sim_inv())  # Asegura que la simulación se haya ejecutado correctamente
    hist(sim_inv(), breaks = 30, col = "lightblue", 
         main = "Distribución Transformada Inversa",
         xlab = "Valores", ylab = "Frecuencia")
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
    req(sim_rechazo())  # Asegura que la simulación se haya ejecutado correctamente
    hist(sim_rechazo(), breaks = 30, col = "lightgreen", 
         main = "Distribución Aceptación y Rechazo",
         xlab = "Valores", ylab = "Frecuencia")
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
    hist(sim_composicion(), breaks = 30, col = "lightcoral", 
         main = "Distribución Composición",
         xlab = "Valores", ylab = "Frecuencia")
  })
  
  output$statsComposicion <- renderPrint({
    req(sim_composicion())
    x <- sim_composicion()
    cat("Media:", mean(x), "\nVarianza:", var(x))
  })
  
}
