source("Isotopes.R")

ph <- 1.5574E-4
getLowerAndUpperLimits <- function(neh,
                                   bwe,
                                   maxLabelingDuration,
                                   k,
                                   peptideSequence,
                                   delta_I0,
                                   kt) {
  # Initialize sequences for x and y
  
  I0_0 = Isotope_Distribution(peptideSequence)[1]
  
  ph <- 1.5574E-4
  time <- seq(1, maxLabelingDuration, 1)
  I0_asymp = I0_0 * (1 - bwe / (1 - ph))^neh
  theo_ub = I0_0 - delta_I0
  theo_lb = as.numeric((1 - kt) * I0_0 + kt * I0_asymp)
  
  
  
  # Initialize an empty dataframe to store the results with the correct column names
  result_df <- data.frame(
    time = integer(0),
    theo_lb_values = numeric(0),
    theo_ub_values = numeric(0),
    lb = numeric(0),
    ub = numeric(0),
    InRange = logical(),
    stringsAsFactors = FALSE
  )
  
  numberOfDigits = 4
  # Iterate over each value of x and y, computing their product and storing it in the dataframe
  for (t_i in time) {
    exp_lb = I0_asymp + (I0_0 - I0_asymp) * exp(-k * t_i)
    exp_ub = I0_asymp + (I0_0 - I0_asymp) * exp(-k * t_i)
    
    
    result_df <- rbind(
      result_df,
      data.frame(
        time = as.integer(t_i),
        theo_lb_values = format(theo_lb, digits = numberOfDigits) ,
        theo_ub_values = format(theo_ub, digits = numberOfDigits) ,
        lb = format(exp_lb, digits = numberOfDigits) ,
        ub = format(exp_ub, digits = numberOfDigits) ,
        InRange = (exp_lb >= theo_lb & exp_ub <= theo_ub)
      )
    )
  }
  
  return (result_df)
}


getLowerAndUpperLimits_Range <- function(neh,
                                         bwe,
                                         maxLabelingDuration,
                                         k_low,
                                         k_high,
                                         peptideSequence,
                                         delta_I0,
                                         kt) {
  # Initialize sequences for x and y
  
  I0_0 = Isotope_Distribution(peptideSequence)[1]
  ph <- 1.5574E-4
  time <- seq(1, maxLabelingDuration, 1)
  I0_asymp = I0_0 * (1 - bwe / (1 - ph))^neh
  theo_ub = I0_0 - delta_I0
  theo_lb = as.numeric((1 - kt) * I0_0 + kt * I0_asymp)
  
  # Initialize an empty dataframe to store the results with the correct column names
  result_df <- data.frame(
    time = integer(0),
    theo_lb_values = numeric(0),
    theo_ub_values = numeric(0),
    lb_low = numeric(0),
    ub_low = numeric(0),
    lb_upper = numeric(0),
    ub_upper = numeric(0),
    InRange_low = logical(),
    InRange_upper = logical(),
    stringsAsFactors = FALSE
  )
  
  numberOfDigits = 4
  # Iterate over each value of x and y, computing their product and storing it in the dataframe
  for (t_i in time) {
    I0_t_low_k = I0_asymp + (I0_0 - I0_asymp) * exp(-k_low * t_i)
    I0_t_high_k = I0_asymp + (I0_0 - I0_asymp) * exp(-k_high * t_i)
    
    
    result_df <- rbind(
      result_df,
      data.frame(
        time = as.integer(t_i),
        theo_lb_values = format(theo_lb, digits = numberOfDigits) ,
        theo_ub_values = format(theo_ub, digits = numberOfDigits) ,
        
        
        I0_t_low_k = format(I0_t_low_k, digits = numberOfDigits) ,
        I0_t_high_k = format(I0_t_high_k, digits = numberOfDigits) ,
        
        
        
        InRange_low = (I0_t_low_k >= theo_lb &
                         I0_t_low_k <= theo_ub),
        InRange_upper = (I0_t_high_k >= theo_lb &
                           I0_t_high_k <= theo_ub)
        
      )
    )
  }
  
  return (result_df)
}
