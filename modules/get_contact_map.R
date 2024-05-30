#!/usr/bin/env Rscript

###################################################
### INFORMATION ###
###################################################
# get_contact_map.R
# KCollier 28 May 2024

# This is a command-line-usable script that outputs a synteny/contact map for two aligned genomes.
# It is generally assumed the output comes from RagTag.
# tidyverse is a dependency

###################################################
### PARSE ARGUMENTS ###
###################################################
default_args <- list(help = FALSE,
                     paf = NULL,
                     query_lab = "Query",
                     target_lab = "Target",
                     raw = FALSE,
                     min_len = 0,
                     contig_plots = FALSE,
                     prefix = "output")

args <- R.utils::commandArgs(trailingOnly = TRUE,
                             asValues = TRUE,
                             defaults = default_args)

print_help <- function() {
  sink(stdout(), type = "message")
  message("\nget_contact_map.R, Keiler Collier, 28 May 2023, v.0.2.0)
  
  Suggested usage: get_contact_map.R --paf <input_paf> --prefix <output_prefix>
          
  Options: 
          --raw: If set, all alignments are reported in the dotplot and quality filtering is not performed. Any micro/macro bounding remains.
                 It is not recommended to set this, as secondary, short and low-mapping-quality alignments can confound interpretation.
                 If --raw is unset (default behavior), secondary mappings, alignments under 500kbp and alignments with a mapq > 40 are removed.
         
          --target_lab <'Target'>: The label of the genome you are mapping the query sequences to. Defaults to 'Target'. 
                 More generally, this is the name of your 'truth set'. Reads in Target that don't have an alignment with Query will not be included in the dotplot.
          
          --query_lab <'Query'>: The label of the genome that you are attempting to find matches in. Defaults to 'Query'.
                 This is your genome-of-interest. All reads in Query will be displayed, regardless of matches with Target.
          
          --min_len <0>: Takes an integer value and sets the minimum alignment length for an alignment to be mapped. If unset, no micro/macro analyses will be done.
                Defaults to 0. Suggested boundary between avian macro/microchromosomes is 10 000 000bp. (10000000)
          
          --contig_plots: If set, 'n' number of contig-wise alignment plots will be written, where 'n' is the number of high-quality plots in 'Query'.
                IT IS NOT RECOMMENDED TO RUN THIS WITH --raw.
          ")
  sink(NULL, type = "message")
  quit()
}

# Sanity check for required inputs:
if (is.null(args$paf)) {
  message("Required inputs missing. See --help for suggested usage.")
  stop("get_contact_map.R --paf <input_paf> --prefix <output_prefix>
       OR
       get_contact_map.R --target <input_fasta> --query <input_fasta> --prefix <output_prefix>", call. = FALSE)
}

###################################################
### CHECKS FOR ALL NECESSARY PACKAGES ###
###################################################
### CHECK AND LOAD LIBS
check_installed_libs <- function(libs) {
  message("Checking for uninstalled libraries...")
  new_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
  if(length(new_libs)) {
    message("Missing packages found. Installing.")
    install.packages(new_libs)
  }
  message("")
}
lib_load <- function(libs) {
  # Function to smooth out library loading process
  message("Loading libraries:")
  for (i in 1:length(libs)) { 
    suppressPackageStartupMessages(library(libs[i], character.only = TRUE))
    message(paste("\t",toString(libs[i]),"successfully loaded"))
  }
  message("All packages loaded\n")
}

libs <- c("pafr", "ggplot2", "magrittr", "dplyr", "stringr", "scales") # list of libs
# Load libraries present in system
check_installed_libs(libs)
lib_load(libs)

print_params <- function(args) {
  message("User-defined parameters:")
  
  for (i in 1:length(args)) {
    name <- names(args[i])
    if (!name == "") {
      message(stringr::str_c("\t",name, ": ", toString(args[i])))
    }
  }
  message("")
}
print_params(args)


###################################################
### MAIN FUNCTION ###
###################################################
# Read in PAF
message("Reading in Pairwise mApping Format (PAF): ")
paf <- pafr::read_paf(args$paf)

# Filter PAF, if 'raw' is unset (default behavior):
# removes secondary alignments, alignments under 500 000 bp, and alignments with mapq under 40 are removed.
if (args$raw == FALSE) {
  min_clean_aln_len = 5e5 #500000 bp. 
  # Reasonable for most chromosomal-level arrangementsl; may exclude a couple of micro
  # We can optimize later.
  paf <- paf %>%
    dplyr::filter(tp == "P" & alen >= min_clean_aln_len & mapq > 40)
}
paf <- dplyr::arrange(paf, desc(qlen)) # sort by size

###################################################
### PLOT ALIGNMENTS ###
###################################################
save_plot <- function(plot, filename, filetype, out_dir = getwd()) {
  filetype = filetype
  ggsave(
    filename = file.path(out_dir,filename),
    plot = plot,
    device = filetype,
    width = 500,
    height = 250,
    units = "mm",
    dpi = 300
  )
}

message("Plotting alignments:")
message("")
prefix_plotname <- stringr::str_replace(args$prefix, "_", " ")
prefix_filename <- stringr::str_replace(args$prefix, " ", "_")

# Define theme
assembly_theme <- theme(
  title = element_text(size = 18, face = "bold"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  legend.position = "none")

# Scatterplot of alignment length
scatterplot <- paf %>%
  ggplot(aes(x = alen, y = dv, fill = tp)) +
  geom_point(shape = 21, color = "black", cex = 2) +
  labs(title = stringr::str_c(prefix_plotname," alignment length scatterplot")) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_percent()) +
  xlab("Alignment length (bp)") +
  ylab("Alignment divergence") +
  theme_bw() +
  assembly_theme

# Dotplot of all alignments:
dotplot_all <- pafr::dotplot(paf, label_seqs = TRUE, order_by = "qstart") + 
  labs(title = stringr::str_c(prefix_plotname, " dotplot all")) +
  theme_bw() +
  assembly_theme

# If 'min_len' has been set, make dotplots for micro and macro chromosomes
if (args$min_len != 0) {
  # Macrochromosomes
  paf_macro <- paf %>%
    filter(alen >= args$min_len)
  dotplot_macro <- pafr::dotplot(paf_macro, label_seqs = TRUE, order_by = "qstart") +
    labs(title = stringr::str_c(prefix_plotname, " dotplot macro")) +
    theme_bw() +
    assembly_theme
  
  #Microchromosomes
  paf_micro <- paf %>%
    filter(alen < args$min_len) 
  dotplot_micro <- pafr::dotplot(paf_micro, label_seqs = TRUE, order_by="qstart") +
    labs(title = stringr::str_c(prefix_plotname, " dotplot micro")) +
    theme_bw() +
    assembly_theme
}

# get individual plots of all macro+micro, if we set min_len
# FIXME - integrate this with micro/macro information if min_len is set
if (args$contig_plots) {
  out_dir <- stringr::str_c(getwd(),stringr::str_c(args$prefix,"alignment_plots",sep = "_"),sep="/")
  dir.create(out_dir)
  
  if (args$min_len != 0 & FALSE) {
    # FIXME - have some way of tagging which are micro and which are macro
    # Code is effectively commented out right now
  } else {
    for (i in 1:length(paf$alen)) {
      message(stringr::str_c("\tPlotting alignment ", i, "/", length(paf$alen), ": "))
      message(stringr::str_c("\t\t", paf$tname[i], " and ", paf$qname[i]))
      
      # plot the damn thing
      plot_title <- stringr::str_c(args$prefix, "alignment", paf$tname[i], paf$qname[i], sep = " ")
      
      aln_plot <- pafr::plot_synteny(paf, q_chrom = paf$qname[i], t_chrom = paf$tname[i]) +
        labs(title = plot_title) +
        theme_bw() +
        assembly_theme
      
      # get the plot title and filename title
      file_prefix <- stringr::str_c(args$prefix,"aligned",i,paf$tname[i],paf$qname[i],sep = "_")
      file_name <- stringr::str_c(file_prefix, "png", sep = '.')
      
      # save the plot
      save_plot(plot = aln_plot, 
                filename = stringr::str_c(file_prefix, "png", sep = '.'), 
                filetype = "png", 
                out_dir = out_dir)
    }
  }
}

###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
# Save scatterplot of alignments
save_plot(plot = scatterplot, filename = stringr::str_c(prefix_filename, "_alnlen_scatterplot.png"), filetype = "png")

# Contact map of all alignments
save_plot(plot = dotplot_all, filename = stringr::str_c(prefix_filename, "_all_alns.png"), filetype = "png")
# Contact map of micro and macro chromosomes, if applicable
if (args$min_len != 0) {
  #macro
  save_plot(plot = dotplot_macro, filename = stringr::str_c(prefix_filename, "_macro_alns.png"), filetype = "png")
  #micro
  save_plot(plot = dotplot_micro, filename = stringr::str_c(prefix_filename, "_micro_alns.png"), filetype = "png")
}

# Save list of individual contigs:
message(stringr::str_c("Output written to ", getwd()))
message("Program terminating")
