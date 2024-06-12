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
                     q_lab = "Query",
                     t_lab = "Target",
                     raw = FALSE,
                     min_macro = 0, # FIXME - this could be set to just, the greater 50% distribution of the contigs
                     contig_plots = FALSE,
                     output_paf = TRUE,
                     prefix = NULL)

args <- R.utils::commandArgs(trailingOnly = TRUE,
                             asValues = TRUE,
                             defaults = default_args)

if (args$help) {
  print_help <- function() {
    sink(stdout(), type = "message")
    message("\nget_contact_map.R, Keiler Collier, 28 May 2023, v.0.2.0)
  
  Suggested usage: get_contact_map.R --paf <input_paf> --prefix <output_prefix>
  
  This is a script designed to visualize whole-genome alignments, filtering out low-quality or spurious alignments as necessary. 
  By default, it outputs a scatterplot of alignment length and percent divergence (<prefix>_alnlen_scatterplot.png) and a contact map of
  all query sequences matched to any target sequences (<prefix>_all_alns.png). Target sequences that are not aligned with any query 
  sequences are not reported- if these are of interest, it is recommended to run and visualize two separate alignments.
  By default, a filtered .paf file is also output.
  If the --min_macro option is set to a nonzero value, two additional contact maps are produced- <prefix>_macro_alns.png, for alignments 
  smaller than the value of --min_macro, and <prefix>_macro_alns.png, for all alignments greater than it.
  
  All outputs are written to your working directory. See options below for further detail.

  
  Options: 
          --raw: If set, all alignments are reported in the dotplot and quality filtering is not performed. Any micro/macro bounding remains.
                 It is not recommended to set this, as secondary, short and low-mapping-quality alignments can confound interpretation.
                 If --raw is unset (default behavior), secondary mappings, alignments under 500kbp and alignments with a mapq > 40 are removed.
         
          --min_macro <0>: Takes an integer value and sets the minimum alignment length for an alignment to be mapped. If unset, no micro/macro analyses will be done.
                Defaults to 0. Suggested boundary between avian macro/microchromosomes is 10 000 000bp. Also accepts scientific notation (eg., 1e7)
          
          --contig_plots: If set, 'n' number of contig-wise alignment plots will be written, where 'n' is the number of high-quality plots in 'Query'.
                Running with --raw set to 'TRUE' will suppress this argument.
                
          --output_paf: This option outputs a csv of all filtered/plotted alignments. Turned on by default; set to 'FALSE' to suppress.
                          
          --prefix: The output prefix used for any plots and figures. Defaults to the prefix of your input.paf file.
          
          --t_lab: The label used for the query aligment of the contact maps. Defaults to \"Target\".

          --q_lab: The label used for the query aligment of the contact maps. Defaults to \"Query\".
          ")
    sink(NULL, type = "message")
    quit()
  }
  # Actually invoke the function
  print_help()
}

# Sanity check for required inputs:
if (is.null(args$paf)) {
  message("Required inputs missing. See --help for suggested usage.")
  stop("get_contact_map.R --paf <input_paf> --prefix <output_prefix>", call. = FALSE)
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

libs <- c("pafr", "ggplot2", "magrittr", "dplyr", "stringr", "scales", "readr") # list of libs
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

# Sanity check for missing arguments:
if (is.null(args$paf)) { 
  message("PAF file not found. Supply with --paf")
  quit()
}
# silently set prefix to name of PAF is prefix isn't set
  # FIXME - this doesn't take into account filestrings as input, causing problems in labeling.
if (is.null(args$prefix)) { 
  args$prefix = unlist(stringr::str_split(tail(unlist(stringr::str_split(args$paf, '/')), 1), fixed('.')))[1]
  }
print_params(args)


###################################################
### MAIN FUNCTION ###
###################################################
SAM_check <- function(paf_obj, SAM_tag) {
  # this is a function to sanity-check whether or not an optional SAM-formatted tag is present in your PAF
  #SAM_tags = c("am","as","bc","bq","cc","cm","co","cp","cq","cs","ct","e2","fi","fs","fz","lb","h0","h1","h2","hi","ih","mc","md","mq","nh","nm","oq","op","oc","pg","pq","pt","pu","qt","q2","r2","rg","rt","sa","sm","tc","u2","uq")
  return(SAM_tag %in% colnames(paf_obj))
}

# Read in PAF
message("Reading in Pairwise mApping Format (PAF): ")
paf <- pafr::read_paf(args$paf)
# Filter PAF, if 'raw' is unset (default behavior):
# removes secondary alignments, alignments under 500 000 bp, and alignments with mapq under 40 are removed.
if (args$raw == FALSE) {
  min_clean_aln_len = 5e5 #500000 bp. 
  # Reasonable for most chromosomal-level arrangements; may exclude a couple of micro
  # We can optimize later.
  paf <- paf %>% 
    dplyr::filter(alen >= min_clean_aln_len & mapq > 40)
  
  # if SAM flags are present, filter them here
  if (SAM_check(paf, "tp")) {
    paf <- paf %>%
      dplyr::filter(tp == "P")
  }
}
paf <- dplyr::arrange(paf, desc(qlen)) # sort by size

# Write CSV of filtered alignments
if (args$output_paf) {
  write_paf <- function(paf_df = paf) {
    # Have to relabel the data frame's columns to it works downstream
    if (SAM_check(paf, "tp")) {
      paf_formatted <- paf_df %>% 
        dplyr::mutate(tp = stringr::str_c("tp:A:",tp),
                      cm = stringr::str_c("cm:i:",cm),
                      s2 = stringr::str_c("s2:i:",s2),
                      dv = stringr::str_c("dv:f:",dv),
                      rl = stringr::str_c("rl:i:",rl))
    } else { paf_formatted = paf_df }
    
    readr::write_delim(x = dplyr::as_tibble(paf_formatted), 
                       file = stringr::str_c(args$prefix,"alignments.paf", sep = "_"), 
                       delim = "\t",
                       append = FALSE,
                       col_names = FALSE)
  }
  write_paf(paf) # We've defined write_paf as a function so we're not storing an extra data frame anywhere in the main function
}

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
if (SAM_check(paf, "tp")) {
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
}

# Dotplot of all alignments:
# FIXME - deconstruct 'pafr::dotplot' so you can make it look better and see how it works
dotplot_all <- pafr::dotplot(paf, label_seqs = TRUE, order_by = "qstart", xlab = args$q_lab, ylab = args$t_lab) + 
  labs(title = stringr::str_c(prefix_plotname, " dotplot all")) +
  theme_bw() +
  assembly_theme

# If 'min_macro' has been set, make dotplots for micro and macro chromosomes
if (args$min_macro != 0) {
  # Macrochromosomes
  paf_macro <- paf %>%
    filter(alen >= args$min_macro)
  if (nrow(paf_macro) == 0) {
    message(stringr::str_c("No alignments greater than", args$min_macro, "found", sep = " "))
  } else {
    dotplot_macro <- pafr::dotplot(paf_macro, label_seqs = TRUE, order_by = "qstart", xlab = args$q_lab, ylab = args$t_lab) +
      labs(title = stringr::str_c(prefix_plotname, " dotplot macro")) +
      theme_bw() + 
      assembly_theme
  }

  #Microchromosomes
  paf_micro <- paf %>%
    filter(alen < args$min_macro) 
  if (nrow(paf_micro) == 0) {
    message(stringr::str_c("No alignments smaller than", args$min_macro, "found", sep = " "))
  } else {
    dotplot_micro <- pafr::dotplot(paf_micro, label_seqs = TRUE, order_by="qstart", xlab = args$q_lab, ylab = args$t_lab) +
      labs(title = stringr::str_c(prefix_plotname, " dotplot micro")) +
      theme_bw() +
      assembly_theme
  }
}

# get individual plots of all macro+micro, if we set min_macro
if (args$contig_plots) {
  message("Producing individual contig plots: \n")
  if (args$raw == TRUE) {
    message(
    "ERROR: You have enabled --raw, which is incompatible with contig plots.\nThis is because it will produce potentially thousands of images, most of which are of low-quality alignments.\nTo get contig plots, turn off --raw.\n")

  } else {
    message("blah")
    out_dir <- stringr::str_c(getwd(),stringr::str_c(args$prefix, "alignment_plots", sep = "_"),sep="/")
    dir.create(out_dir)
    
    if (args$min_macro != 0 & FALSE) {
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
}

###################################################
### SAVE PLOTS TO LOCAL FILESYSTEM ###
###################################################
# Save scatterplot of alignments
if (SAM_check(paf, "tp")) {
  save_plot(plot = scatterplot, filename = stringr::str_c(prefix_filename, "_alnlen_scatterplot.png"), filetype = "png")
}

# Contact map of all alignments
save_plot(plot = dotplot_all, filename = stringr::str_c(prefix_filename, "_all_alns.png"), filetype = "png")
# Contact map of micro and macro chromosomes, if applicable
if (args$min_macro != 0) {
  #macro
  if (exists("dotplot_macro")) {
    save_plot(plot = dotplot_macro, filename = stringr::str_c(prefix_filename, "_macro_alns.png"), filetype = "png")
  }
  #micro
  if (exists("dotplot_micro")) {
    save_plot(plot = dotplot_micro, filename = stringr::str_c(prefix_filename, "_micro_alns.png"), filetype = "png")
  }
}

# Save list of individual contigs:
message(stringr::str_c("Output written to '", getwd(), args$prefix),"*'")
message("Program terminating")
