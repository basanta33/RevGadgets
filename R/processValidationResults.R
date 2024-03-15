#' Process Validation Results 
#' 
#' Reads in and processes output for simulation based calibration.
#' calculates HPD intervals and coverage probabilities 
#'
#'
#' @param analysis_name (character string; no default) Name of the analysis
#' @param n_reps (single numeric value; default: 1000)
#'  number of replicates of the validation 
#' @param num_bins (single numeric value; default: 50) number of bins for 
#'  coverage probability
#' 
#' @return A data frame with the coverage probabilities
#' 
#' @export

processValidation <- function(analysis_name, n_reps = 1000, n_bins = 50) {
  if (missing(analysis_name)) {
    stop("Missing required argument: Analysis Name")
  }
  
  # Define directories
  results_dir = paste0("results_", analysis_name)
  output_dir = paste0("output_", analysis_name)
  
  # create directories
  dir.create(results_dir, showWarnings = FALSE)
  
  # function to read data from the file 
  read_data <- function(file_path) {
    if (file.exists(file_path)) {
      return(read.table(file_path, 
                        sep = "\t", 
                        header = TRUE, 
                        skip = 0, 
                        check.names = FALSE))
    } else {
      return(NULL)
    }
  }
  
  # get the list of parameters 
  parameters <- colnames(read_data(paste0(output_dir, 
                                          "/Validation_Sim_0/posterior_samples.var")))
  parameters <- parameters[-1]   # remove the first column ("Iteration")
  
  # exclude "branch_rates" if present 
  parameters <- parameters[parameters != "branch_rates"]

  # Print list of parameters 
  cat("parameters:\n", parameters, "\n\n")
  
  # Iterate over each parameter
  for (param in parameters) {
    # initialize variables 
    coverage_probs <- data.frame(total_count = numeric(0), 
                                 in_count = numeric(0), 
                                 hpd_width = numeric(0), 
                                 stringsAsFactors = FALSE)
    hpd_width <- seq(from = 0.0, to = 1.0, by = 1/n_bins)
    
    for (i in 1:(n_bins+1)) {
      coverage_probs[i,] = c(total_count = 0, in_count = 0, hpd_width = hpd_width[i])
    }
    
    # initialize progress bar 
    pb <- txtProgressBar(min = 0, max = n_reps, char = "*", style = 3)
    
    # iterate over each replication 
    for (i in 1:n_reps) {
      setTxtProgressBar(pb, i)
      
      # read in the data 
      data <- read_data(paste0(output_dir, 
                         "/Validation_Sim_", 
                         i-1, 
                         "/posterior_samples.var"))
      if (is.null(data)) next
      
      # extract samples 
      num_samples = length(data[,1])
      
      x <- as.mcmc(data[round(0.25*num_samples):num_samples, param])
      
      true_val_ext <- ifelse(param == "branch_rates", ".out", ".txt")
      true_val <- read.table(file=paste0(output_dir, 
                                         "/Validation_Sim_", 
                                         i-1, 
                                         "/", 
                                         param, 
                                         true_val_ext))[1,1]
      
      # calculate coverage probabilities 
      for (k in 1:(n_bins + 1)) {
        hpd <- HPDinterval(x, prob = hpd_width[k])
        if (true_val >= hpd[1,1] && true_val <= hpd[1,2]){
          coverage_probs[k, "in_count"] <- coverage_probs[k, "in_count"] + 1
        }
        coverage_probs[k, "total_count"] <- coverage_probs[k, "total_count"] + 1
      }      
    }
    close(pb)
    
    # calculate frequency of coverage
    coverage_probs$freq = coverage_probs$in_count / coverage_probs$total_count
    
    # save coverage probabilities 
    saveRDS(coverage_probs, file = paste0(results_dir, "/", param, ".rds"))
    
    # print results to the screen 
    cat(param,"\n")
    cat("HPD-width:\t\t",hpd_width,"\n")
    cat("Coverage-freq:\t\t",coverage_probs$freq,"\n")
  }
  
  return(coverage_probs)
}


alpha_tachy <- processValidation("tachy")
