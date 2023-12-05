library(ggplot2)
library(tidyverse)
library(readr)
library(lessR)
library(psych)
library(dplyr)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(car)

## function to filter uorfs by reads, start codon, and overlap
filter_uorfs <- function(finput, foutput) {
  output <- read_csv(finput, col_names = TRUE)
  output <- na.omit(output)
  exp = colnames(output)[startsWith(colnames(output), "Experiment")]# NEED TO CHANGE
  
  output_filt <- output %>%
    filter(rowSums(across(all_of(exp), .fns = ~. >= 10, .names = "count_")) > 0)
  
  output_filt$Start <- as.numeric(as.character(output_filt$Start))
  output_filt$'Actual Stop' <- as.numeric(as.character(output_filt$'Actual Stop'))
  
  output_filt$length <- output_filt$'Actual Stop' - output_filt$Start
  
  output_filt_2 <- data.frame() 
  output_filt_filt <- data.frame()
  genes = unique(output_filt$Gene)
  
  for (gene in genes) {
    subset_gene = output_filt[output_filt$Gene == gene,]
    stops = unique(subset_gene$'Actual Stop')
    for (stop in stops) {
      subset_frame = subset_gene[subset_gene$'Actual Stop' == stop, ]
      # for(frame in unique(subset_stop$in_frame)) {
      #   subset_frame = subset_stop[subset_stop$in_frame == frame, ]
        filtered_df <- subset_frame %>%
          group_by(across(all_of(exp))) %>%
          filter(Start == max(Start))
        
        atg_df <- filtered_df[filtered_df$`Start Codon` == "ATG",]
        if(nrow(atg_df > 0)) {
          filtered_df <- atg_df
        }
        
        condition <- ""
        for (experiment_name in exp) {
          condition <- paste0(condition, 
                              "filtered_df$", experiment_name, 
                              " >= max(filtered_df$", experiment_name, ") * 0.8 & ")
        }
        
        condition <- substr(condition, 1, nchar(condition) - 3) 
        
        new_df <- filtered_df[eval(parse(text = condition)), ]
        
        if (nrow(new_df) > 0) {
          output_filt_filt <- rbind(output_filt_filt, new_df)
        } else {
          output_filt_filt <- rbind(output_filt_filt, filtered_df)
        }
        
        output_filt_2 <- rbind(output_filt_2, filtered_df)
    }
  }
  
  write.csv(output_filt_filt[complete.cases(output_filt_filt), ], file = foutput, row.names = FALSE)
}

# TO RUN:
filter_uorfs("INPUT.csv", "OUTPUT.csv")