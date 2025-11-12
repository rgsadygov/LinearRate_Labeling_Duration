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
    #================== labeling  Range (Page 2) ======
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
          verbatimTextOutput("page2_status"),
          
          # Output: A tabset that combines three panels ----
          navset_tab(
            # title = "Visualizations",
            # Panel with plot ----
            nav_panel("Plot", plotlyOutput("page2_plot")),
            
            # Panel with summary ----
            nav_panel("Summary", tableOutput("page2_summaryTable")),
            
            # # Panel with table ----
            # nav_panel("Table", tableOutput("table"))
          ),
          
          tags$p(
            "Remark: The labeling duration should be set so that I0(t) is between the theoretical upper bound (red line) and the lower bound (black line)."
          ),
          verbatimTextOutput("page2_recommendation"),
          
        )
      )
    ),
    
    
    #==================================================
    #================= labeling duraiton (page 1) =====
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
    
    #==================================================
    #================= neh (page 3) =====
    tabPanel("NEH", mainPanel(
      fluidRow(
        column(
          6,
          numericInput(
            "page3_neh_Alanine",
            label = "Alanine (A)",
            min = 0.00,
            max = 15,
            value = 4.00,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Cysteine",
            label = "Cysteine (C)",
            min = 0.00,
            max = 15,
            value = 1.62,
            step = 0.01
          ),
          numericInput(
            "page3_neh_AsparticAcid",
            label = "Aspartic Acid (D)",
            min = 0.00,
            max = 15,
            value = 1.89,
            step = 0.01
          ),
          numericInput(
            "page3_neh_GlutamicAcid",
            label = "Glutamic Acid (E)",
            min = 0.00,
            max = 15,
            value = 3.95,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Phenylalanine",
            label = "Phenylalanine (F)",
            min = 0.00,
            max = 15,
            value = 0.32,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Glycine",
            label = "Glycine (G)",
            min = 0.00,
            max = 15,
            value = 2.06,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Histidine",
            label = "Histidine (H)",
            min = 0.00,
            max = 15,
            value = 2.88,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Isoleucine",
            label = "Isoleucine (I)",
            min = 0.00,
            max = 15,
            value = 1.00,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Lysine",
            label = "Lysine (K)",
            min = 0.00,
            max = 15,
            value = 0.54,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Leucine",
            label = "Leucine (L)",
            min = 0.00,
            max = 15,
            value = 0.60,
            step = 0.01
          )
        ),
        column(
          6,
          numericInput(
            "page3_neh_Methionine",
            label = "Methionine (M)",
            min = 0.00,
            max = 15,
            value = 1.12,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Asparagine",
            label = "Asparagine (N)",
            min = 0.00,
            max = 15,
            value = 1.89,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Proline",
            label = "Proline (P)",
            min = 0.00,
            max = 15,
            value = 2.59,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Glutamine",
            label = "Glutamine (Q)",
            min = 0.00,
            max = 15,
            value = 3.95,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Arginine",
            label = "Arginine (R)",
            min = 0.00,
            max = 15,
            value = 3.43,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Serine",
            label = "Serine (S)",
            min = 0.00,
            max = 15,
            value = 2.61,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Threonine",
            label = "Threonine (T)",
            min = 0.00,
            max = 15,
            value = 0.20,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Valine",
            label = "Valine (V)",
            min = 0.00,
            max = 15,
            value = 0.56,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Tryptophan",
            label = "Tryptophan (W)",
            min = 0.00,
            max = 15,
            value = 0.08,
            step = 0.01
          ),
          numericInput(
            "page3_neh_Tyrosine",
            label = "Tyrosine (Y)",
            min = 0.00,
            max = 15,
            value = 0.42,
            step = 0.01
          )
        )
      ),
      br(),
      div(
        style = "text-align:center;",
        actionButton("page3_save_neh", "Save", class = "btn btn-danger btn-lg")
      ),
      br(),
    ))
    
    
    #==================================================
  )
)
