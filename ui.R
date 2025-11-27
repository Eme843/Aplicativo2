library(shiny)
library(kableExtra)
library(data.table)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(shinyjs)

# Define UI for application that draws a histogram
fluidPage(
  
  # Application title
  titlePanel(HTML("<h1 style='font-family: Poppins; font-weight:bold; color: #1F4E79;'>Números Aleatorios</h1>")),
  
  tabsetPanel(
    tabPanel("Números Aleatorios",
             
             # Sidebar with a slider input for number of bins
             sidebarLayout(
               sidebarPanel(
                 selectInput("metodo", "Seleccione el Método", 
                             choices = c("Congruencial Multiplicativo", "Congruencial Mixto", 
                                         "Cuadrados Medios", "Lehmer"),
                             selected = "Congruencial Multiplicativo"),
                 
                 actionButton("seleccionar", "Seleccionar", class = "btn-primary", 
                              style = "color : #F0FFF0; background-color: #1E90FF; border-color: #1E90FF", 
                              icon("check")),
                 
                 conditionalPanel(
                   condition = "input.seleccionar > 0",
                   uiOutput("parametros_ui"),
                 
    
                   actionButton("mostrar", "Mostrar resultados", class = "btn-danger", 
                                style = "color : #F0FFF0; background-color: #EE0000; border-color: #F0FFF0", 
                                icon("eye"))
                 )
                ),
               # Show a plot of the generated distribution
               mainPanel(
                 h4(style = "font-weight: bold;", "Números aleatorios generados:"),
                 br(),
                 conditionalPanel(condition = "input.mostrar!=0",
                                  div(style = "max-height: 300px; overflow-y: scroll;",
                                      uiOutput("tabla")
                                  ),
                                  br(),
                                  h4(style = "font-weight: bold;", "Gráfico:"),
                                  fluidRow(
                                    column(width = 8,
                                           plotOutput("distPlot")
                                           ),
                                    column(width = 4,
                                           sliderInput("barras", "Número de barras para el histograma:", 
                                                       min = 5, max = 50, value = 10)
                                    )
                                    )
                                   )
               )
             )
    ),
    tabPanel("Integrales",
             useShinyjs(),
             sidebarLayout(
               sidebarPanel(
                 textInput("funcion", "Ingrese la funcion a integrar:", value = "1-x"),
                 
                 fluidRow(
                   column(7, numericInput("lim_inf", "Límite inferior del intervalo:", value = 0)),
                   column(5, checkboxInput("inf_inf", "Límite inferior: -∞", value = FALSE))
                 ),
                 
                 fluidRow(
                   column(6, numericInput("lim_sup", "Límite superior del intervalo:", value = 1)),
                   column(6, checkboxInput("sup_inf", "Límite superior: +∞", value = FALSE))
                 ),
                 
                 radioButtons("metodo", "Seleccione el método de generación de números aleatorios:",
                              choices = c("Congruencial Multiplicativo",
                                          "Congruencial Mixto",
                                          "Cuadrados Medios",
                                          "Lehmer")),
                 
                 actionButton("calcular", "Calcular Área",
                              class = "btn-lg btn-success",
                              style = "color : #F0FFF0; background-color: #9ACD32; border-color: #8B8B83",
                              icon("chart-area"))
               ),
               
               mainPanel(
                 conditionalPanel(
                   condition = "input.calcular!=0",
                   
                   fluidRow(
                     column(
                       width = 8, offset = 2,
                       h4("Gráfica de la función a integrar:"),
                       div(
                         style = "max-width: 100%; height: 350px; 
                   overflow: hidden; 
                   border: 1px solid #ddd; 
                   padding: 10px; 
                   margin-bottom: 25px;",
                         plotOutput("graf_fun01", height = "320px")
                       )
                     )
                   ),
                   fluidRow(
                     column(
                       width = 8, offset = 2,
                       h4("Aproximación:"),
                       verbatimTextOutput("valores_mc"),
                       uiOutput("nota_metodo"),
                       
                       br(),
                       div(
                         style = "max-width: 100%; height: 350px; 
                   overflow: hidden; 
                   border: 1px solid #ddd; 
                   padding: 10px;",
                         plotOutput("graf_aprox01", height = "320px")
                       )
                     )
                   )
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
    ),
    tabPanel("Variables Aleatorias Continuas")
  )
)
