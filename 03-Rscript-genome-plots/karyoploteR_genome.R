#!/usr/bin/env Rscript

# Function to install missing packages
install_if_missing <- function(packages) {
  installed <- installed.packages()
  repos <- getOption("repos")
  if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org/"))
  }
  for (pkg in packages) {
    if (!pkg %in% rownames(installed)) {
      if (pkg %in% c("BiocManager", "karyoploteR", "regioneR", "zoo", "rtracklayer", "AnnotationDbi")) {
        if (!"BiocManager" %in% rownames(installed)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# List of required packages
required_packages <- c("karyoploteR", "regioneR", "zoo", "BiocManager", "rtracklayer", "AnnotationDbi")

# Install missing packages
install_if_missing(required_packages)

# Load required libraries
library(karyoploteR)
library(regioneR)
library(zoo)
library(rtracklayer)

# Function to display help information
print_help <- function() {
  cat("
Usage: Rscript plot_karyotype.R <custom_genome_name> <custom_genome_end> <gff_file_path>

Description:
  This script generates karyotype plots using the karyoploteR library. It requires a custom genome name, the end position of the genome, and a GFF file path containing gene annotations.

Arguments:
  <custom_genome_name>  : The name of the assembly, must match what is specified in the gff  (e.g., 'cy0333_M3_MHC_v2-genes').
  <custom_genome_end>   : The end position of the assembly (e.g., 5227476).
  <gff_file_path>       : The file path to the GFF file containing gene annotations.

Example:
  Rscript plot_karyotype.R \"cy0333_M3_MHC_v2-genes\" 5227476 \"/path/to/your/file.gff\"

")
}

# Function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3 || any(args %in% c("-h", "--help"))) {
    print_help()
    quit(status = 1)
  }
  list(
    custom_genome_name = args[1],
    custom_genome_end = as.numeric(args[2]),
    gff_file_path = args[3]
  )
}

# Parse command line arguments
args <- parse_args()

# STEP 2 SET CUSTOM GENOME
custom.genome <- toGRanges(data.frame(chr=c(args$custom_genome_name), start=c(1), end=c(args$custom_genome_end)))

# STEP 3 LOAD GFF
features <- import(args$gff_file_path)
genes <- features[features$type == "gene"]

# STEP 4 PLOT PARAMETERS AND BEGIN PLOT
pp <- getDefaultPlotParams(plot.type=2)
pp$topmargin=20
pp$bottomargin=20
pp$leftmargin <- 0.05
pp$data1inmargin=0
pp$data2inmargin=0
pp$ideogramheight=0

# Start plot and save as PDF
pdf("plot_gene_annotations.pdf", width=11, height=8)
kp <- plotKaryotype(genome = custom.genome, plot.type = 2, plot.params = pp, labels.plotter=NULL)

# Separate + and - strand annotations
positive <- genes[strand(genes)=="+"]
negative <- genes[strand(genes)=="-"]

# Plot gene annotation intervals
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE, r1=0.05, data.panel=1)
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, r1=0.05, data.panel=2)

# Plot gene name annotations, IMPORTANT: edit label.dist!! 
kpPlotMarkers(kp, data=positive, max.iter=1000, labels=positive$Name, data.panel=1, text.orientation = "vertical", cex=0.4, label.dist=-0.0049, r0=0, r1=0.25, marker.parts = c(0.4,0.4, 0.2), label.margin = 0.1)
kpPlotMarkers(kp, data=negative, max.iter=1000, labels=negative$Name, data.panel=2, text.orientation = "vertical", cex=0.4, label.dist=-0.0049, r0=0, r1=0.25, marker.parts = c(0.4,0.4, 0.2), label.margin = 0.1)
dev.off()

# STEP 5 CREATE TICKMARKS
# Plot parameters
pp <- getDefaultPlotParams(plot.type=2)
pp$topmargin=20
pp$bottomargin=20
pp$leftmargin <- 0.05
pp$data1inmargin=0
pp$data2inmargin=0
pp$ideogramheight=0

# Make sure this genome is the same one from above
custom.genome <- toGRanges(data.frame(chr=c(args$custom_genome_name), start=c(0), end=c(args$custom_genome_end)))

# Start tickmarks plot and save as PDF
pdf("plot_tickmarks.pdf", width=11, height=8)
kp <- plotKaryotype(genome = custom.genome, plot.type = 2, plot.params = pp, labels.plotter=NULL)
kpAddBaseNumbers(kp, tick.dist = 250000, tick.len = 5, cex=0.5, minor.tick.dist = 50000, minor.tick.len = 2, srt=90)
dev.off()
