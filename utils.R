library(plotly)
library(dplyr)
library(tidyr)


# Constants
pw <- 0.046
ph <- 1.5574e-4

# Function to compute I0_t
# inputs theoretical relative monoisotpe abundance: I0_0,
# I0 value at the plateau of labeling:  I0_asymp, turnover rate:k, and labeling
# duration:t
get_I0_t <- function(I0_0, I0_asymp, k, t) {
  return(I0_asymp + (I0_0 - I0_asymp) * exp(-k * t))
}

# Function to compute  relative monoisotpe abundance value at the plateau of
# labeling inputs theoretical relative monoisotpe abundance: I0_0
# and number of exchangeable hydrogens: neh
get_I0_asymp <- function(I0_0, neh) {
  return(I0_0 * ((1 - pw / (1 - ph))^neh))
}

# funtion to compute the numberof exchangable hydrogens for a given
# peptide sequence
getNeh <- function(peptide) {
  aa_neh = list(
    "A" = 4   ,
    "C" = 1.62,
    "D" = 1.89,
    "E" = 3.95,
    "F" = 0.32,
    "G" = 2.06,
    "H" = 2.88,
    "I" = 1   ,
    "K" = 0.54,
    "L" = 0.6 ,
    "M" = 1.12,
    "N" = 1.89,
    "P" = 2.59,
    "Q" = 3.95,
    "R" = 3.43,
    "S" = 2.61,
    "T" = 0.2 ,
    "V" = 0.56,
    "W" = 0.08,
    "Y" = 0.42
  )
  
  
  res = 0
  for (aa in strsplit(peptide, "")[[1]]) {
    res = res + as.numeric(aa_neh[as.character(aa)])
  }
  
  return (res)
  
}

# Function to create a plot showing experimental and theoretical I0 bounds 
# over time and recommended labeling duration
getPossibleRangePlot_plotly <- function(ans) {
  
  # Convert columns data type to numeric values
  ans$time <- as.numeric(ans$time)
  ans$lb <- as.numeric(ans$lb)
  ans$ub <- as.numeric(ans$ub)
  ans$theo_lb_values <- as.numeric(ans$theo_lb_values)
  ans$theo_ub_values <- as.numeric(ans$theo_ub_values)
  
  # Reshape data for easier plotting
  ans_long <- ans %>%
    pivot_longer(
      cols = c(lb, ub, theo_lb_values, theo_ub_values),
      names_to = "type",
      values_to = "value"
    )
  
  # Identify region where theoretical values fall within experimental bounds
  possibleRange <- ans %>% filter(InRange == TRUE)
  min_time <- min(possibleRange$time) - 0.5
  max_time <- max(possibleRange$time)
  max_y <- max(possibleRange$theo_ub_values)
  min_y <- min(possibleRange$theo_lb_values)
  
  # Define rectangle for visualizing valid I0 range area
  rect_x <- c(min_time,
              max_time + 2.5E-1,
              max_time + 2.5E-1,
              min_time,
              min_time)
  rect_y <- c(min_y, min_y, max_y, max_y, min_y)
  
  
  #plot
  p <- plot_ly(data = ans_long,
               x = ~ time,
               y = ~ value) %>%
    
    # green rectangle showing valid I0 region
    add_trace(
      x = rect_x,
      y = rect_y,
      fill = "toself",
      fillcolor = "green",
      opacity = 0.1,
      line = list(color = "rgba(0,0,0,0)"),
      type = "scatter",
      mode = "lines",
      showlegend = FALSE,
      hoverinfo = "none",
      inherit = FALSE
    ) %>%
    
    # experimental I0(t) as black markers
    add_trace(
      data = subset(ans_long, type == "lb"),
      type = "scatter",
      mode = "markers",
      marker = list(size = 8, color = "black"),
      name = " ",
      # "Lower bound (exp.)",
      showlegend = FALSE
    ) %>%
    
    #theoretical lower bound as black dotted line
    add_trace(
      data = subset(ans_long, type == "theo_lb_values"),
      type = "scatter",
      mode = "lines",
      line = list(dash = "dot", color = "black"),
      name = "Lower bound (theo.)",
      showlegend = TRUE
    ) %>%
    
    # theoretical lower bound as red dotted line
    add_trace(
      data = subset(ans_long, type == "theo_ub_values"),
      type = "scatter",
      mode = "lines",
      line = list(dash = "dot", color = "red"),
      name = "Upper bound (theo.)",
      showlegend = TRUE
    ) %>%
    
    # Layout styling
    layout(
      title = "",
      #Lower and Upper Bound Limits",
      xaxis = list(
        title = "Labeling Duration",
        gridcolor = "white",
        zeroline = FALSE,
        linecolor = "black",
        showline = TRUE
      ),
      yaxis = list(
        title = "I<sub>0</sub>(t)",
        range = c(
          min(ans$lb, ans$theo_lb_values) - 0.01,
          max(ans$lb, ans$theo_ub_values) + 0.01
        ),
        gridcolor = "white",
        zeroline = FALSE,
        linecolor = "black",
        showline = TRUE
      ),
      # legend = list(title = list(text = "Legend")),
      legend = list(
        x = 0.5,
        # Moves legend to the top-left
        y = 1,
        # Aligns legend at the top
        xanchor = "left",
        yanchor = "top",
        orientation = "h",
        # Places legend items side by side
        title = NULL  # Removes legend title
      ),
      plot_bgcolor = "white",
      font = list(size = 14),
      title = list(font = list(size = 14))
    ) %>% # Updated titlefont syntax
    
    config(displayModeBar = FALSE)
  
  return(p)
}


# Function to create a plot showing experimental and theoretical I0 bounds 
# over time and recommended range of labeling durations
getPossibleRangePlot_plotly_range <- function(ans, klow, khigh) {
  
  # Convert columns data type to numeric values
  ans$time <- as.numeric(ans$time)
  ans$theo_lb_values <- as.numeric(ans$theo_lb_values)
  ans$theo_ub_values <- as.numeric(ans$theo_ub_values)
  ans$I0_t_low_k <- as.numeric(ans$I0_t_low_k)
  ans$I0_t_high_k <- as.numeric(ans$I0_t_high_k)
  
  # Reshape data for easier plotting
  ans_long <- ans %>%
    pivot_longer(
      cols = c(I0_t_low_k, I0_t_high_k, theo_lb_values, theo_ub_values),
      names_to = "type",
      values_to = "value"
    )
  
  # Identify region where theoretical values fall within experimental bounds
  possibleRange <- ans %>% filter(InRange_low == TRUE &
                                    InRange_upper == TRUE)
  min_time <- min(possibleRange$time) - 0.5
  max_time <- max(possibleRange$time)
  max_y <- max(possibleRange$theo_ub_values)
  min_y <- min(possibleRange$theo_lb_values)
   
  
  # Define rectangle for visualizing valid I0 range area
  rect_x <- c(min_time,
              max_time + 2.5E-1,
              max_time + 2.5E-1,
              min_time,
              min_time)
  rect_y <- c(min_y, min_y, max_y, max_y, min_y)
  
  #plot
  p <- plot_ly(data = ans_long,
               x = ~ time,
               y = ~ value) %>%
    
    # green rectangle showing valid I0 region
    add_trace(
      x = rect_x,
      y = rect_y,
      fill = "toself",
      fillcolor = "green",
      opacity = 0.1,
      line = list(color = "rgba(0,0,0,0)"),
      type = "scatter",
      mode = "lines",
      showlegend = FALSE,
      hoverinfo = "none",
      inherit = FALSE
    ) %>%
    
    # experimental I0(t) for small k as black markers
    add_trace(
      data = subset(ans_long, type == "I0_t_low_k"),
      type = "scatter",
      mode = "markers",
      marker = list(size = 8, color = "black"),
      name = paste("k=", klow),
      showlegend = TRUE
    ) %>%
    
    # experimental I0(t) for large k as black markers
    add_trace(
      data = subset(ans_long, type == "I0_t_high_k"),
      type = "scatter",
      mode = "markers",
      marker = list(size = 8, color = "blue"),
      name = paste("k=", khigh),
      showlegend = TRUE
    ) %>%
    
    #theoretical lower bound as black dotted line
    add_trace(
      data = subset(ans_long, type == "theo_lb_values"),
      type = "scatter",
      mode = "lines",
      line = list(dash = "dot", color = "black"),
      name = "Lower bound (theo.)",
      showlegend = TRUE
    ) %>%
    
    # theoretical lower bound as red dotted line
    add_trace(
      data = subset(ans_long, type == "theo_ub_values"),
      type = "scatter",
      mode = "lines",
      line = list(dash = "dot", color = "red"),
      name = "Upper bound (theo.)",
      showlegend = TRUE
    ) %>%
    
    # Layout styling
    layout(
      title = "",
      #Lower and Upper Bound Limits",
      xaxis = list(
        title = "Labeling Duration",
        gridcolor = "white",
        zeroline = FALSE,
        linecolor = "black",
        showline = TRUE
      ),
      yaxis = list(
        title = "I<sub>0</sub>(t)",
        range = c(
          min(ans$I0_t_low_k, ans$I0_t_high_k) - 0.015,
          max(ans$I0_t_low_k, ans$exp_lb_high) + 0.02
        ),
        gridcolor = "white",
        zeroline = FALSE,
        linecolor = "black",
        showline = TRUE
      ),
      # legend = list(title = list(text = "Legend")),
      legend = list(
        x = 0.5,
        # Moves legend to the top-left
        y = 1,
        # Aligns legend at the top
        xanchor = "left",
        yanchor = "top",
        orientation = "h",
        # Places legend items side by side
        title = NULL  # Removes legend title
      ),
      plot_bgcolor = "white",
      font = list(size = 14),
      title = list(font = list(size = 14))
    ) %>% # Updated titlefont syntax
    
    config(displayModeBar = FALSE)
  
  return(p)
}

solve_for_kt_vals <- function(rd) {
  relative_error <- function(x) {
    if (x <= 0 || x >= 1)
      return(Inf)  # Keep x,kt in (0, 1)
    (-log(1 - x) - x) / x - rd
  }
  # Use uniroot within a safe interval
  result <- tryCatch({
    uniroot(relative_error,
            interval = c(1e-6, 0.999),
            tol = 1e-10)$root
  }, error = function(e) {
    NA  # Return NA if no root found
  })

  return(result)
}
