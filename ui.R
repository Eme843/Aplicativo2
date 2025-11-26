library(shiny)
library(kableExtra)
library(data.table)
library(ggplot2)
library(lubridate)
library(tidyverse)

# Define UI for application that draws a histogram
fluidPage(
  
  # Application title
  titlePanel(HTML("<h1 style='font-family: Poppins; font-weight:bold; color: #1F4E79;'>Números Aleatorios</h1>")),
  
  tabsetPanel(
    tabPanel("Números Aleatorios",
             
             # Sidebar with a slider input for number of bins
             sidebarLayout(
               sidebarPanel(
                 sliderInput("semilla",
                             "Ingrese un valor inicial:",
                             min=1,
                             max=500,
                             value = 30),
                 sliderInput("divisor",
                             "Ingrese un valor de m:",
                             min=5,
                             max=500,
                             value = 37),
                 sliderInput("constante",
                             "Ingrese un valor de a:",
                             min=10,
                             max=500,
                             value = 123),
                 sliderInput("c",
                             "Ingrese un valor de c para el método mixto:",
                             min=1,
                             max=500,
                             value = 200),
                 sliderInput("num",
                             "Cantidad de números a generar:",
                             min=10,
                             max=200,
                             value = 30),
                 actionButton("mostrar", "Mostrar resultados",class = "btn-danger",style = "color : #F0FFF0; background-color: #EE0000; border-color: #F0FFF0", icon("eye"))
               ),
               # 
               # Show a plot of the generated distribution
               mainPanel(
                 conditionalPanel(condition = "input.mostrar!=0",
                                  h4("Tabla de resultados - Metodo congruencial multiplicativo:"),
                                  br(), #generar una linea en blanco
                                  tableOutput("tabla"),
                                  br(),
                                  h4("Tabla de resultados - Metodo congruencial mixto:"),
                                  br(),
                                  tableOutput("tabla2"),
                                  br(),
                                  h4("Gráfico:"),
                                  fluidRow(
                                    column(width=2,
                                           numericInput("barras", "Número de barras:", value=10, min=2, max=50),
                                    ),
                                    column(width=5,
                                           h5("Congruencial Multiplicativo"),
                                           plotOutput("distPlot")
                                           
                                    ),
                                    column(width=5,
                                           h5("Congruencial Mixto"),
                                           plotOutput("HistBin")
                                    )
                                  )
                 )
               )
             )
    ),
    tabPanel("Integrales",
             sidebarLayout(
               sidebarPanel(
                 textInput("funcion", "Ingrese la funcion a integrar:", value = "1-x"),
                 numericInput("lim_inf", "Límite inferior del intervalo:", value = 0),
                 numericInput("lim_sup", "Límite superior del intervalo:", value = 1),
                 radioButtons("metodo", "Seleccione el metodo para generar los numeros aleatorios:",
                              c("Congruencial Multiplicativo", "Congruencial Mixto")),
                 actionButton("calcular", "Calcular Área",class = "btn-lg btn-success", style = "color : #F0FFF0; background-color: #9ACD32; border-color: #8B8B83", icon("chart-area"))
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 conditionalPanel( condition = "input.calcular!=0", #Boton tiene por defecto valor de 0
                                   h4("Gráfica de la funcion a integrar:"),
                                   plotOutput("graf_fun01"),
                                   h4("Aproximaciòn:"),
                                   plotOutput("graf_aprox01")
                 )
               )
             )
    ),
    tabPanel("Variables Aleatorias Discretas",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("distribucion", "Seleccione el método:",
                              c("Binomial", "Poisson", "Transformada Inversa", "Aceptación y Rechazo", "Composición")),
                 actionButton("siguiente", "Seleccionar", class = "btn-primary"),
                 
                 # Panel de parámetros dinámico
                 uiOutput("parametrosUI")
               ),
               mainPanel(
                 # Resultados para Binomial
                 conditionalPanel(
                   condition = "input.distribucion == 'Binomial' & input.siguiente > 0 & input.calc_bin > 0",
                   h4("Resultados - Binomial"),
                   plotOutput("histBin"),
                   verbatimTextOutput("statsBin")
                 ),
                 
                 # Resultados para Poisson
                 conditionalPanel(
                   condition = "input.distribucion == 'Poisson' & input.siguiente > 0 & input.calc_pois > 0",
                   h4("Resultados - Poisson"),
                   plotOutput("histPois"),
                   verbatimTextOutput("statsPois")
                 ),
                 
                 # Resultados para Transformada Inversa
                 conditionalPanel(
                   condition = "input.distribucion == 'Transformada Inversa' & input.siguiente > 0 & input.calc_inv > 0",
                   h4("Resultados - Transformada Inversa"),
                   plotOutput("histInv"),
                   verbatimTextOutput("statsInv")
                 ),
                 
                 # Resultados para Aceptación y Rechazo
                 conditionalPanel(
                   condition = "input.distribucion == 'Aceptación y Rechazo' & input.siguiente > 0 & input.calc_rechazo > 0",
                   h4("Resultados - Aceptación y Rechazo"),
                   plotOutput("histRechazo"),
                   verbatimTextOutput("statsRechazo")
                 ),
                 
                 # Resultados para Composición
                 conditionalPanel(
                   condition = "input.distribucion == 'Composición' & input.siguiente > 0 & input.calc_composicion > 0",
                   h4("Resultados - Composición"),
                   plotOutput("histComposicion"),
                   verbatimTextOutput("statsComposicion")
                 ),
                 
                 # Mensaje cuando no hay resultados
                 conditionalPanel(
                   condition = "(input.distribucion == 'Poisson' & input.siguiente > 0 & input.calc_pois == 0) | 
                      (input.distribucion == 'Binomial' & input.siguiente > 0 & input.calc_bin == 0) |
                      (input.distribucion == 'Transformada Inversa' & input.siguiente > 0 & input.calc_inv == 0) |
                      (input.distribucion == 'Aceptación y Rechazo' & input.siguiente > 0 & input.calc_rechazo == 0) |
                      (input.distribucion == 'Composición' & input.siguiente > 0 & input.calc_composicion == 0)",
                   wellPanel(
                     h4("Instrucciones:"),
                     p("Por favor ingrese los parámetros y presione 'Calcular' para ver los resultados")
                   )
                 )
               )
             )
    )
  ),
)
