#!/usr/bin/env Rscript

# Copyright 2023 Sara Llorente-Armijo 

# Rscript to make plot of the enrichment over background for different genomic regions from one or several .tsv file/s (from script `sc_ChIPQC_genomicProfile.R`) and save plot as .pdf

# Arguments of the Rscript are specified with a parser object
# Load "argparse" library
library(argparser, quietly = T)
# Create parser object
parser <- arg_parser("Make plot of the enrichment over background for different genomic regions from one or several .tsv file/s (from script `sc_ChIPQC_genomicProfile.R`) and save plot as .pdf.")
# Specify desired arguments
parser <- add_argument(parser, "--tsv", help = "Input one or several .tsv files separated by spaces. Name with the factor, tissue and replicate. E.g.: --tsv H3K36me3_Endo_rep1_regi.tsv H3K36me3_Endo_rep2_regi.tsv. You can input all files and then select the specific factor with --factor", type = "character", short = "-t", nargs='+')
parser <- add_argument(parser, "--factor", help = "Factor to make plot", type = "character", short = "-f")
parser <- add_argument(parser, "--outdir", help = "Path to folder where to save the output. Default: working directory", default = ".", short = "-d")
parser <- add_argument(parser, "--outfile", help = "Output file name. Default: <factor>_regi_plot.pdf", type = "character", short = "-o")


# E.g.: singularity exec --bind /home/sl6820/mnt/scratch/sl6820 /home/sl6820/mnt/network/DevRegGen/singularity/bioconductor_docker_RELEASE_3_17.sif Rscript --vanilla workflow/scripts/sc_ChIPQC_genomicProfile_plot.R --tsv outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep1_regi.tsv,outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep2_regi.tsv,outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep3_regi.tsv,outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep1_regi.tsv,outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep2_regi.tsv,outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep3_regi.tsv --outfile H3K36me3_regi_plot.pdf

# Parse the command line arguments
argv <- parse_args(parser)

# Output directory
# Check if it exists, or create it
outdir <- argv$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

message("\nOutput directory: ", outdir, "\n")


# Load libraries
library(dplyr)
library(ggplot2)

options(scipen = 99999999)


#################################################
# ChIPQC - Plot enrichment over genomic regions #
#################################################

#==========#
# Get data #
#==========#

message("\n######################")
message("Getting .tsv file/s...")
message("######################\n")

#------------------#
# Load .tsv file/s #
#------------------#

# tsv <- "mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep1_regi.tsv"

# tsv <- "mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep2_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep3_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep2_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep3_regi.tsv"

# tsv <- "mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/Input_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/Input_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me1_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me1_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me2_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me2_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K4me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K9me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K9me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K27ac_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K27ac_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K27me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K27me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K36me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K79me3_Endo_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/H3K79me3_Mesecto_rep1_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/Input_Endo_rep2_regi.tsv,mnt/scratch/sl6820/projects/project_XlaeData/process_ChIP/outputs/chipqc/genomicEnrichment/Input_Mesecto_rep2_regi.tsv"


tsv <- argv$tsv

# Get .tsv file/s
tsv_list <- (tsv %>% strsplit(","))[[1]]
names(tsv_list) <- gsub("_regi.tsv", "", basename(tsv_list))

# Take only those for specified factor
tsv_list_factor <- tsv_list[grep(argv$factor, names(tsv_list))]
tsv_list_factor

# Load .tsv file/s
chipqc_regi_list <- list()
for (i in names(tsv_list_factor)) {
  chipqc_regi_list[[i]] <- read.table(file = tsv_list_factor[[i]], row.names = 1, col.names = c("", i))
}

# Join tables
chipqc_regi_all <- do.call(cbind, chipqc_regi_list) %>%
  tibble::rownames_to_column(var = "GenomicIntervals")

#--------------------------------#
# Create data frame for plotting #
#--------------------------------#

message("\n###################################")
message("Creating data frame for plotting...")
message("###################################\n")

# Join data and organise columns
regiScoresFrame <- reshape2::melt(chipqc_regi_all, id.vars = "GenomicIntervals", variable.name = "Sample", value.name = "log2_Enrichment") %>%
  tidyr::separate(Sample, sep = "_", into = c("Factor", "Tissue", "Replicate"), remove = F)

regiScoresFrame[,"GenomicIntervals"] <- factor(regiScoresFrame[,"GenomicIntervals"],levels=unique(as.vector(regiScoresFrame[,"GenomicIntervals"])))



#=====================================#
# Plot enrichment for genomic regions #
#=====================================#

message("\n###########")
message("Plotting...")
message("###########\n")

# Plot

plot <- ggplot(regiScoresFrame) +
  geom_tile(aes(y=Replicate,x=GenomicIntervals,fill = log2_Enrichment)) +
  facet_grid(Tissue~Factor, switch = "y") +
  scale_fill_gradient2(low="blue",high="yellow",mid="black",midpoint=median(regiScoresFrame$log2_Enrichment)) +
  scale_y_discrete(position = "right") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.ticks.x = element_line(colour = "black"),
        # axis.ticks.length.x = unit(0.1,"cm"),
        axis.text.x = element_text(size = 11, angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 12),
        strip.text.y.left = element_text(size = 11, angle = 0),
        legend.title = element_text(size = 10),
        legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.justification = "top",
        legend.direction = "vertical",
        plot.margin = unit(c(0, 0.1, 0, 0.5), "inches"))

plot

# Save
ggsave(filename = file.path(argv$outdir, paste0(argv$factor, "_regi_plot.pdf")), plot = plot, device = "pdf", width = 6, height = 4)
message("\nOutput saved in ", file.path(argv$outdir, paste0(argv$factor, "_regi_plot.pdf")))

