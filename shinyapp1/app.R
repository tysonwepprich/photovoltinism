# Shiny App for simulating photothermographs and diapause decisions
# Compare two scenarios with different species parameters/years/sites

# Setup ----
library(shiny)
library(ggplot2)
library(viridis)
library(dplyr)
library(purrr)
library(daymetr)
library(ggridges)
library(shinycssloaders)
source('helpers.R')

# prefix is used to define the two scenarios to simulate
renderInputs <- function(prefix) {
  wellPanel(
    fluidRow(
      column(6,
             sliderInput(paste0(prefix, "_", "nsim"), "Number of individuals to simulate:", min = 1, max = 500, value = 100, step = 1),
             sliderInput(paste0(prefix, "_", "year_range"), "Range of observations (in Years):", min = 1980, max = 2016, value = c(2005, 2010), step = 1, sep = ""),
             numericInput(paste0(prefix, "_", "lat"), "Site latitude (in decimal units, limited to North America):", min = 6.6, max = 83.3, value = 44.56, step = .1),
             numericInput(paste0(prefix, "_", "lon"), "Site longitude (in decimal units, limited to North America):", min = -178.2, max = -49, value = -123.26, step = .1),
             sliderInput(paste0(prefix, "_", "dev_temp"), "Temperature range for development (in Celsius):", min = -10, max = 45.0, value = c(10.0, 37.0), step = 0.1),
             sliderInput(paste0(prefix, "_", "cdl_mu"), "Critical photoperiod (hours of light):", min = 8, max = 24, value = 15.5, step = 0.1),
             sliderInput(paste0(prefix, "_", "cdl_sd"), "Std. dev. critical photoperiod:", min = 0, max = 3, value = 0.5, step = 0.1),
             sliderInput(paste0(prefix, "_", "lambda"), "Population growth rate between generations:", min = .5, max = 5, value = 2, step = 0.1)
             
      ),
      column(6,
             sliderInput(paste0(prefix, "_", "emerg_dd"), "Mean overwinter emergence (in degree-days):", min = 0.0, max = 2000, value = 100, step = 1),
             sliderInput(paste0(prefix, "_", "emerg_dd_sd"), "Std. dev. overwinter emergence:", min = 0.0, max = 50, value = 15, step = 1),
             sliderInput(paste0(prefix, "_", "gen_dd"), "Mean development time egg to teneral adult (in degree-days):", min = 100, max = 2000, value = 367, step = 1),
             sliderInput(paste0(prefix, "_", "gen_dd_sd"), "Std. dev. development time egg to teneral adult:", min = 0.0, max = 50, value = 30, step = 1),
             sliderInput(paste0(prefix, "_", "povip_dd"), "Mean pre-oviposition time (in degree-days):", min = 0, max = 500, value = 126, step = 1),
             sliderInput(paste0(prefix, "_", "povip_dd_sd"), "Std. dev. pre-oviposition time:", min = 0.0, max = 50, value = 10, step = 1),
             sliderInput(paste0(prefix, "_", "ovip_dd"), "Mean time until oviposition (in degree-days, Poisson distribution):", min = 0, max = 200, value = 30, step = 1)
      )
    ),
    p(actionButton(paste0(prefix, "_", "recalc"),
                   "Run simulation", icon("refresh"))
    )
    
    
    # TODO: reactive options for plots once simulation runs
    # see conditionalPanel function
  )
} # close renderInputs


# Define UI ----
ui <- fluidPage(
  
  # Define UI for application that plots random distributions
  fluidPage(theme="simplex.min.css",
            tags$style(type="text/css",
                       "label {font-size: 12px;}",
                       ".recalculating {opacity: 1.0;}"
            ),
            
            # Application title
            tags$h2("Development, voltinism and photoperiod-based diapause decisions simulated for insects"),
            p("An adaptation of the model from ",
              tags$a(href="https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/14-2071.1", "Grevstad & Coop 2015"),
              "with individual variation in critical photoperiod and development requirements"),
            hr(),
            
            # TODO: Add more instructions/explanations above steps for user input
            
            fluidRow(
              column(6, tags$h3("Scenario A")),
              column(6, tags$h3("Scenario B"))
            ),
            fluidRow(
              column(6, renderInputs("a")),
              column(6, renderInputs("b"))
            ),
            # fluidRow(
            #   column(6, uiOutput("a_year_plot")
            #   ),
            #   column(6, uiOutput("b_year_plot")
            #   )
            # ),
            fluidRow(
              column(6,
                     withSpinner(plotOutput("a_distPlot", height = "600px"))
              ),
              column(6,
                     withSpinner(plotOutput("b_distPlot", height = "600px"))
              )
            ),
            br(),
            fluidRow(
              column(6,
                     withSpinner(plotOutput("a_conseqPlot", height = "600px"))
              ),
              column(6,
                     withSpinner(plotOutput("b_conseqPlot", height = "600px"))
              )
            )
  )
  
  
) # close ui fluidPage


# TODO: 
# If plot options modified (after simulations run), modify plot output only



temp_param_names <- c("year_range", "lat", "lon")
gdd_param_names <- c("dev_temp")
sim_param_names <- c("nsim", "cdl_mu", "cdl_sd", "lambda", "emerg_dd", 
                     "emerg_dd_sd", "gen_dd", "gen_dd_sd", "povip_dd", 
                     "povip_dd_sd", "ovip_dd")


# Define server logic ----

# Define server logic required to generate and plot a random distribution
#
# Idea and original code by Pierre Chretien
# Small updates by Michael Kapler
#
server <- function(input, output, session) {
  
  
  get_temp_params <- function(prefix) {
    params <- lapply(temp_param_names, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    names(params) <- temp_param_names
    params
  }
  
  get_gdd_params <- function(prefix) {
    params <- lapply(gdd_param_names, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    names(params) <- gdd_param_names
    params
  }
  
  get_sim_params <- function(prefix) {
    params <- lapply(sim_param_names, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    names(params) <- sim_param_names
    params
  }
  
  # Delay simulation until action button pressed
  # If accum_gdd or nsim or traits modified AND action button = run simulations and output plots
  
  # output$a_distPlot <- renderPlot({
  #   if (input$a_recalc == 0){
  #     return()
  #   }
  #   isolate({ 
  #     temp_a <- get_temp_data(userinput = get_temp_params("a"))
  #     gdd_a <- get_gdd_data(temp_data = temp_a, userinput = get_gdd_params("a"))
  #     results_a <- run_sims(pars = get_sim_params("a"), gdd = gdd_a)
  #     
  #     # # debug without reactive
  #     # temp_a <- get_temp_data(userinput = list("year_range" = c(2005, 2010), "lat" = 45, "lon" = -100))
  #     # gdd_a <- get_gdd_data(temp_data = temp_a, userinput = list("dev_temp" = c(10, 30)))
  #     # results_a <- run_sims(pars = as.list(pars), gdd = gdd_a)
  # 
  #     })
  #   plot_sim(results = results_a, gdd = gdd_a) 
  #   
  # })
  
  # Scenario A
  
  temp_a <- eventReactive(input$a_recalc, {
    result <- get_temp_data(userinput = get_temp_params("a"))
    return(result)
  })
  
  gdd_a <- eventReactive(input$a_recalc, {
    result <- get_gdd_data(temp_data = temp_a(), userinput = get_gdd_params("a"))
    return(result)
  })
  
  results_a <- eventReactive(input$a_recalc, {
    result <- run_sims(pars = get_sim_params("a"), gdd = gdd_a())
    return(result)
  })
  
  output$a_distPlot <- renderPlot({
    if (input$a_recalc == 0){
      return()
    }
    
    # on hold: selecting plot years after simulation
    # a_plot_yrs <- reactive({ifelse(exists("input$a_plot_yrs") == FALSE, 
    #                      unique(results_a()$year),
    #                      input$a_plot_yrs)})
    # isolate({ 
      plot_sim(results = results_a(), gdd = gdd_a()) 
    # })
  })
  # # on hold: selecting plot years after simulation
  # # Issues with selecting years for plots
  # # Below, none selected upon start
  # output$a_year_plot <- renderUI({
  #   yrs <- unique(results_a()$year)
  #   selectizeInput("a_plot_yrs", "Years to plot", 
  #                  choices = yrs, multiple = TRUE)
  # })
  output$a_conseqPlot <- renderPlot({
    if (input$a_recalc == 0){
      return()
    }
    plot_conseq(results = results_a(), gdd = gdd_a()) 
  })
  
  # Scenario B
  
  temp_b <- eventReactive(input$b_recalc, {
    result <- get_temp_data(userinput = get_temp_params("b"))
    return(result)
  })
  
  gdd_b <- eventReactive(input$b_recalc, {
    result <- get_gdd_data(temp_data = temp_b(), userinput = get_gdd_params("b"))
    return(result)
  })
  
  results_b <- eventReactive(input$b_recalc, {
    result <- run_sims(pars = get_sim_params("b"), gdd = gdd_b())
    return(result)
  })
  
  output$b_distPlot <- renderPlot({
    if (input$b_recalc == 0){
      return()
    }
    
    plot_sim(results = results_b(), gdd = gdd_b()) 
  })

  output$b_conseqPlot <- renderPlot({
    if (input$b_recalc == 0){
      return()
    }
    plot_conseq(results = results_b(), gdd = gdd_b()) 
  })
}

#TODO: adjust ggridge axis to match gdd on fig1
# get rid of y axis labels, redundant with legend

# Run the app ----
shinyApp(ui = ui, server = server)