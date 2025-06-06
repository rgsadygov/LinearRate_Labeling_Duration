library(ggplot2)
library(plotly)
library(bslib)
library(shiny)
library(shinythemes)

library(shinyjs)

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  useShinyjs(),
  #cosmo
  navbarPage(
    "Turnover rate estimation",
    #==================================================
    #==================== labeling duraiton============
    tabPanel(
      "Labeling duration estimation",
      sidebarLayout(
        sidebarPanel(
          textInput('page1_peptide', label = "Peptide", value = "VGAFTVVCK"),
          numericInput(
            "page1_expected_k",
            label = "Expected turnover rate",
            min = 0.00,
            max = 0.9,
            value = 0.05,
            step = 0.01
          ),
          numericInput(
            "page1_bwe",
            label = "Body water Enrichment",
            min = 0.01,
            max = 0.9,
            value = 0.046,
            step = 0.01
          ),
          numericInput(
            "page1_delta_I0",
            label = (HTML(paste0("∆I", tags$sub(
              "0"
            )))),
            min = 0.01,
            max = 0.9,
            value = 0.031,
            step = 0.01
          ),
          numericInput(
            "page1_max_error",
            label = "Maximum k relative error (%)",
            min = 0,
            max = 100,
            value = 25,
            step = 0.5
          ),
          sliderInput(
            "page1_maxTime",
            "Max. Labeling Duration",
            min = 1,
            max = 100,
            value = 24
          ),
          
          actionButton("page1_submit", "Submit", class = "btn  btn-danger")
        ),
        mainPanel(
          verbatimTextOutput("page1_status"),
          # Output: A tabset that combines three panels ----
          navset_tab(
            # title = "Visualizations",
            # Panel with plot ----
            nav_panel("Plot", plotlyOutput("page1_plot")),
            
            # Panel with summary ----
            nav_panel("Summary", tableOutput("page1_summaryTable")),
             
          ),
          
          tags$p(
            "Remark: The labeling duration should be set so that I0(t) is between the theoretical upper bound (red line) and the lower bound (black line)."
          ),
          verbatimTextOutput("page1_recommendation"),
          
        )
      )
    ),
    
    
    #==================================================
    #==================== labeling  Range============
    tabPanel(
      "Labeling duration estimation [Range]",
      sidebarLayout(
        sidebarPanel(
          textInput('page2_peptide', label = "Peptide", value = "AQTAHIVLEDGTK"),
          numericInput(
            "page2_expected_k_low",
            label = "Expected slowest turnover rate",
            min = 0.00,
            max = 0.9,
            value = 0.06,
            step = 0.01
          ),
          numericInput(
            "page2_expected_k_high",
            label = "Expected highest turnover rate",
            min = 0.00,
            max = 0.9,
            value = 0.12,
            step = 0.01
          ),
          numericInput(
            "page2_bwe",
            label = "Body water Enrichment",
            min = 0.01,
            max = 0.9,
            value = 0.046,
            step = 0.01
          ),
          numericInput(
            "page2_delta_I0",
            label = (HTML(paste0("∆I", tags$sub(
              "0"
            )))),
            min = 0.01,
            max = 0.9,
            value = 0.031,
            step = 0.01
          ),
          numericInput(
            "page2_max_error",
            label = "Maximum relative error (%)",
            min = 0,
            max = 100,
            value = 25,
            step = 0.5
          ),
          sliderInput(
            "page2_maxTime",
            "Max. Labeling Duration",
            min = 1,
            max = 100,
            value = 24
          ),
          
          actionButton("page2_submit", "Submit", class = "btn  btn-info")
        ),
        mainPanel(
          verbatimTextOutput("simstatus_3"),
          # plotlyOutput("simplot"),
          # h3(""),
          # tableOutput("simdatatable")
          
          # Output: A tabset that combines three panels ----
          navset_tab(
            # title = "Visualizations",
            # Panel with plot ----
            nav_panel("Plot", plotlyOutput("simplot_3")),
            
            # Panel with summary ----
            nav_panel("Summary", tableOutput("simdatatable_3")),
            
            # # Panel with table ----
            # nav_panel("Table", tableOutput("table"))
          ),
          
          tags$p(
            "Remark: The labeling duration should be set so that I0(t) is between the theoretical upper bound (red line) and the lower bound (black line)."
          ),
          verbatimTextOutput("simstatus_recommendation_3"),
          
        )
      )
    ),
    
    #==================================================
    
  )
)
