
source("utils.R")
source("Labeling _duration_bounds.R")


library(munsell)
library(ggplot2)
library(plotly)
library(dplyr)
library(shiny)
library(tidyr)
library(shinyjs)
source("utils.R") 
source("Labeling _duration_bounds.R")


# Define server logic 
server <- function(input, output, session) {
 
   observe({
    # Code to execute on page load
    cat("App has loaded!\n")
    runjs('$(document).ready(function() {$("#page1_submit").click();});')
    runjs('$(document).ready(function() {$("#page2_submit").click();});')
    
  })
  
  #====================================================
  #page 1
  #====================================================
  
  page1_inputdata <- reactive({
    dict <- data.frame(
      key = c("bwe", "time", "peptide", "neh", "k", "delta_I0", "kt"),
      value = c(
        input$page1_bwe,
        input$page1_maxTime,
        input$page1_peptide,
        getNeh(input$page1_peptide),
        input$page1_expected_k,
        input$page1_delta_I0,
        solve_for_kt_vals(input$page1_max_error / 100)
      )
    )
    output <- dict
  })
  
  
  output$page1_plot <- renderText(if (input$page1_submit > 0)
  {
    temp = isolate(page1_inputdata())
    
    if (temp$value[temp$key == "kt"] == Inf) {
      return (
        "Unable to determine the labeling durations for the given inputs. Please adjust the input parameters."
      )
    }
    else{
      return (paste(c(
        paste(c(
          "Peptide: ", temp$value[temp$key == "peptide"], "(NEH", (temp$value[temp$key == "neh"]), ")"
        )), paste(c(" Max. Days: ", temp$value[temp$key == "time"])), paste(c(
          " kt: ", sprintf("%.4f", as.numeric(temp$value[temp$key == "kt"]))
        )), paste(c(" BWE: ", temp$value[temp$key == "bwe"]))
      )
      , sep = ","))
    }
  }
  else{
    return ("Server is ready for calculation.")
  })
  
  
  
  simulateGetLowerAndUpperLimits <- reactive({
    temp_input = (page1_inputdata())
    temp_res = getLowerAndUpperLimits(
      as.numeric(temp_input$value[temp_input$key == "neh"]),
      as.numeric(temp_input$value[temp_input$key == "bwe"]),
      as.numeric(temp_input$value[temp_input$key == "time"]),
      as.numeric(temp_input$value[temp_input$key == "k"]),
      (temp_input$value[temp_input$key == "peptide"]),
      as.numeric(temp_input$value[temp_input$key == "delta_I0"]),
      as.numeric(temp_input$value[temp_input$key == "kt"])
    )
    
    output <- temp_res
  })
  
  
  output$page1_summaryTable <- renderTable ({
    if (input$page1_submit > 0)
    {
      temp_res = isolate(simulateGetLowerAndUpperLimits())
      
      temp_res %>%
        rename(
          "Labeling Duration" =    "time",
          "Lower Limit (theo.)" =  "theo_lb_values",
          "Upper Limit (theo.)" =   "theo_ub_values",
          "Lower Limit (exp.)" =   "lb",
          "Upper Limit (exp.)" =   "ub",
          "In Range" =              "InRange"
        )
      
      
      
      
    }
  })
  
  
  
  output$page1_plot <- renderPlotly(if (input$page1_submit > 0)
  {
    temp_res = isolate(simulateGetLowerAndUpperLimits()) 
    
    temp_plot = getPossibleRangePlot_plotly(temp_res) 
    
    # ggplotly(temp_plot)
  })
  
  output$page1_recommendation <- renderText(if (input$page1_submit > 0)
  {
    temp_res = isolate(simulateGetLowerAndUpperLimits())
    possibleRange <- temp_res %>% filter(InRange == TRUE)
    min_time <- min(possibleRange$time)
    max_time <- max(possibleRange$time)
    
    return (paste(
      c(
        "Minimum Labeling Duration",
        min_time,
        "\nMaximum Labeling Duration",
        max_time
      )
    ))
  }
  else{
    return ("")
  })
  
  
  #====================================================
  #====================================================
  #page 2
  
  
  page2_inputdata <- reactive({
    dict <- data.frame(
      key = c(
        "bwe",
        "time",
        "peptide",
        "neh",
        "k_low",
        "k_high",
        "delta_I0",
        "kt"
      ),
      value = c(
        input$page2_bwe,
        input$page2_maxTime,
        input$page2_peptide,
        getNeh(input$page2_peptide),
        input$page2_expected_k_low,
        input$page2_expected_k_high,
        input$page2_delta_I0,
        solve_for_kt_vals(input$page2_max_error / 100)
      )
    )
    output <- dict
  })
  
  
  output$page2_status <- renderText(if (input$page2_submit > 0)
  {
    temp = isolate(page2_inputdata())
    
    
    return (paste(c(paste(
      c("Peptide: ", temp$value[temp$key == "peptide"], "(NEH", (temp$value[temp$key == "neh"]), ")")
    ), paste(
      c(" Max. Days: ", temp$value[temp$key == "time"])
    ), paste(c(
      " kt: ", sprintf("%.4f", as.numeric(temp$value[temp$key == "kt"]))
    )), paste(
      c(" BWE: ", temp$value[temp$key == "bwe"])
    ))
    , sep = ","))
  }
  else{
    return ("Server is ready for calculation.")
  })
  
  
  
  simulateGetLowerAndUpperLimits_Range <- reactive({
    temp_input = (page2_inputdata())
    temp_res = getLowerAndUpperLimits_Range(
      as.numeric(temp_input$value[temp_input$key == "neh"]),
      as.numeric(temp_input$value[temp_input$key == "bwe"]),
      as.numeric(temp_input$value[temp_input$key == "time"]),
      as.numeric(temp_input$value[temp_input$key == "k_low"]),
      as.numeric(temp_input$value[temp_input$key == "k_high"]),
      (temp_input$value[temp_input$key == "peptide"]),
      as.numeric(temp_input$value[temp_input$key == "delta_I0"]),
      as.numeric(temp_input$value[temp_input$key == "kt"])
    )
    
    # print((temp_res))
    
    output <- temp_res
  })
  
  
  output$page2_summaryTable <- renderTable ({
    if (input$page2_submit > 0)
    {
      temp_input = isolate(page2_inputdata())
      temp_res = isolate(simulateGetLowerAndUpperLimits_Range())
      
      temp_res %>%
        mutate(across(where(is.numeric), ~ format(
          .x, digits = 5, nsmall = 5
        )))
      
      
      
      colname_ub_low = paste("Upper Limit (exp.), k=", as.numeric(temp_input$value[temp_input$key == "I0_t_low_k"]))
      colname_exp_lb_high = paste("Lower Limit (exp.), k=", as.numeric(temp_input$value[temp_input$key == "I0_t_high_k"]))
      colname_InRange_low = paste("In Range, k=", as.numeric(temp_input$value[temp_input$key == "k_low"]))
      colname_InRange_upper =      paste("In Range, k=", as.numeric(temp_input$value[temp_input$key == "k_high"]))
      
      
      temp_res <- temp_res %>%
        rename_with(
          ~ case_when(
            . == "time" ~ "Labeling Duration",
            . == "theo_lb_values" ~ "Lower Limit (theo.)",
            . == "theo_ub_values" ~ "Upper Limit (theo.)", 
            . == "ub_low" ~ colname_ub_low, 
            . == "exp_ub_high" ~ colname_exp_lb_high,
            . == "InRange_low" ~ colname_InRange_low,
            . == "InRange_upper" ~ colname_InRange_upper,
            TRUE ~ .
          )
        )
      
      
    }
  })
  
  
  
  output$page2_plot <- renderPlotly(if (input$page2_submit > 0)
  {
    temp_res = isolate(simulateGetLowerAndUpperLimits_Range())
    
    
    temp_input = isolate(page2_inputdata())
    
    temp_plot = getPossibleRangePlot_plotly_range(temp_res,
                                                  as.numeric(temp_input$value[temp_input$key == "k_low"]),
                                                  as.numeric(temp_input$value[temp_input$key == "k_high"]))
    
    
    # ggplotly(temp_plot)
  })
  
  
  
  output$page2_recommendation <- renderText(if (input$page2_submit > 0)
  {
    temp_res = isolate(simulateGetLowerAndUpperLimits_Range())
    possibleRange <- temp_res %>% filter(InRange_low == TRUE &
                                           InRange_upper == TRUE)
    min_time <- min(possibleRange$time)
    max_time <- max(possibleRange$time)
    
    return (paste(
      c(
        "Minimum Labeling Duration",
        min_time,
        "\nMaximum Labeling Duration",
        max_time
      )
    ))
  }
  else{
    return ("")
  })
  
  
  
  
  
  
  
  
}
