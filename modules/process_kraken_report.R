#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# process_kraken_report.R
# Usage:
# Keiler Collier, 12 May 2024

# This is a script designed to get the taxids of contaminants in a readset found by Kraken2
# We do this by searching for "Bacteria, and then assuming everything aftr is noneukaryotic and a contaminant.
# It returns a text file of all non-eukaryotic taxids found in the dataset.
  # This text file can then be passed to 'fasta_name_filt_OVH' to remove the tagged data.

# Usage notes:
  # Unix commands related to filtering kraken reports:
  # cat $raw_report | grep -v "0.0[0-9]" > $filt_report
  # cat Burhinus_test.fq | grep --color=auto 'kraken:taxid|2$'

###################################################
### PARSE ARGUMENTS ###
###################################################
ARGS <- commandArgs(trailingOnly = TRUE)
# If no argument is supplied, return error:
if (length(ARGS)==0)
{
  stop("Usage: ./process_kraken_report.R <kraken_report.txt> <output_prefix>", call. = FALSE)
} else if(length(ARGS)==1)
{
  ARGS[2] = "output"
}

###################################################
### LOAD PACKAGES AND DEFINE VARIABLES ###
###################################################
# FIXME - implement some kind of checking to make sure these packages are installed, and install if not
library(readr)
library(stringr)
library(dplyr)

# Set variables from ARGS
message("Parameters interpreted as:")
file_string <- ARGS[1]
message(paste(" Kraken2 report:", file_string))
output_prefix <- ARGS[2]
message(paste(" Output prefix:", output_prefix))
message("")

###################################################
### HELPER FUNCTIONS ###
###################################################
# None currently, actually. Pretty simple script.

###################################################
### MAIN FUNCTION ###
###################################################
# Read in file
message("Reading in Kraken2 report: ")
report <- read_tsv(file_string, 
                   col_names = FALSE)
# Add column names
colnames_6 = c("percent_of_reads", "num_reads_clade_nested", "num_reads_clade_direct", "taxon_rank", "taxid", "taxon") #if you used Kraken on reads only
colnames_8 = c("percent_of_reads", "num_reads_clade_nested", "num_reads_clade_direct", "num_mins_clade_nested", "num_mins_clade_direct", "taxon_rank", "taxid", "taxon") #if you included minimizer data
ifelse(length(report) == 6, names(report) <- colnames_6, names(report) <- colnames_8)

# Filter out taxa with no assigned reads:
report <- report %>% filter(percent_of_reads > 0.00)

# FIXME : Cast taxon_rank as factor and provide levels
report <- report %>%
  mutate(taxon_rank = as_factor(taxon_rank))

# Basic order is U, R, D, K, P, C, O, F, G, S.
# figure out a way to account for taxonomic sublevels. Probably use unique(report$taxon_rank) and begins_with()
levels(report$taxon_rank) 


# strategy:
  #Get taxids for all entries under D, until you hit the next D.

contaminant_taxids <- report %>% slice(which(report$taxon == "Bacteria"):nrow(report)) %>% # get taxa for first domain after Eukaryota (Bacteria). taxon_rank will equal "D"
  select(taxid) %>% 
  mutate(taxid = str_c("kraken:taxid|",taxid))

# FIXME : write the data frame we just made of all the contaminants to a text file that C++ can deal with.
filename <- str_c(output_prefix, ".csv")
write_csv(contaminant_taxids, file = filename, col_names = FALSE)

message(str_c("Output written to", getwd()))
message("Program finishing.")







