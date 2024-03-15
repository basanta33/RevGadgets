#' Plot Coverage Probabilities
#' 
#' Plots the output of simulation based calibration validation results
#' 
#' The returned plot is a ggplot object 
#' 
#' @param df (data frame; no default) with coverage probabilities
#' @param analysis_name (character string; no default) name of the analysis
#' @param figures_dir (character string; no default) directory for the figures 
#' @import coda
#' @import ggplot2
#' @export


generate_coverage_plots <- function(analysis_name, 
                                    num_reps = 10000, 
                                    num_bins = 50, 
                                    results_dir, 
                                    figs_dir) {
  library(coda)
  library(ggplot2)
  
  # Create directories if they don't exist
  dir.create(figs_dir, showWarnings = FALSE)
  dir.create(results_dir, showWarnings = FALSE)
  
  # Read data
  in_file <- paste0("output_", analysis_name, "/Validation_Sim_1/posterior_samples.var")
  data <- read.table(in_file, sep="\t", header=TRUE, skip=0, check.names=FALSE)
  parameters <- colnames(data)[-1]
  
  # Iterate over each parameter
  for (param in parameters) {
    # Read coverage probabilities
    coverage_probs <- readRDS(file = paste0(results_dir, "/", param, ".rds"))
    
    # Generate plot
    p <- ggplot(coverage_probs) +
      geom_bar(stat="identity", aes(x=hpd_width, y=freq), colour="lightgray", fill="lightgray") +
      theme_classic() +
      xlab("HPD width") + ylab("coverage probability") + ggtitle(param) +
      geom_segment(aes(x=0, y=0, xend=1, yend=1), linetype="dashed", size=1.5, show.legend=FALSE) +
      theme(legend.position="none", plot.title = element_text(hjust = 0.5))
    
    # Save plot
    ggsave(paste0(figs_dir, "/hpd_width_vs_coverage_", param, ".pdf"), plot=p, width=10, height=10, units="cm")
  }
}

generate_coverage_plots("tachy", results_dir = "results_tachy", figs_dir = "figures_tachy")
