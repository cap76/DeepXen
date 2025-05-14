#!/usr/bin/env Rscript

# Copyright 2023 Sara Llorente-Armijo 

# Rscript to calculate the enrichment over background for different genomic regions for a BAM file and save output as .tsv

# Arguments of the Rscript are specified with a parser object
# Load "argparse" library
library(argparser, quietly = T)
# Create parser object
parser <- arg_parser("Calculate the enrichment over background for different genomic regions for a BAM file and save output as .tsv.")
# Specify desired arguments
parser <- add_argument(parser, "--bam", help = "Input a BAM file", type = "character", short = "-b")
parser <- add_argument(parser, "--outdir", help = "Path to folder where to save the output. Default: working directory", default = ".", short = "-d")
parser <- add_argument(parser, "--outfile", help = "Output file name. Default: chipqc_regi.tsv", type = "character", short = "-o", default = "chipqc_regi.tsv")

parser <- add_argument(parser, "--gtf", help = "GTF file (.gtf) of genome (for genes, promoters, cds and exons", type = "character", short = "-g")
parser <- add_argument(parser, "--gff", help = "Gzipped GFF file (.gff.gz) of genome (for introns and UTRs)", type = "character", short = "-f")
parser <- add_argument(parser, "--centromeres", help = "BED file with coordinates for centromeres", type = "character", short = "-c")
parser <- add_argument(parser, "--chromsizes", help = "File with a column with chromosomes and other with length of chromosomes (.chrom.sizes file)", type = "character", short = "-s")


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
library(plyranges)
library(GenomicFeatures)
library(ChIPQC)

options(scipen = 99999999)


############################################
# ChIPQC - Enrichment over genomic regions #
############################################

#====================================#
# Get genomic regions for annotation #
#====================================#

message("\n###########################################")
message("Preparing genomic regions for annotation...")
message("###########################################\n")

#-------#
# Genes #
#-------#
gtf_file <- argv$gtf # "../../data/genomes/teleost_amphibia/Xlae/Xlae10.1/Sequence/Transcriptome/GCF_017654675.1_Xenopus_laevis_v10.1_genomic_modified_chrNames.gtf"
gtf_all <- rtracklayer::import.gff(gtf_file, colnames = c("type", "gene_id", "transcript_id", "db_xref", "gene", "gene_biotype", "transcript_biotype", "protein_id", "gene_synonym"))
gtf_gene <- rtracklayer::import.gff(gtf_file, feature.type = "gene", colnames = c("type", "gene_id", "transcript_id", "db_xref", "gene", "gene_biotype", "transcript_biotype", "protein_id", "gene_synonym"))


#-----------#
# Promoters #
#-----------#
promoters200_gr <- promoters(gtf_gene, upstream = 0, downstream = 200)

promoters500_gr <- promoters(gtf_gene, upstream = 0, downstream = 500)

promoters200to500_gr <- resize(promoters500_gr, width = 300, fix = "end")

promotersUp1000_gr <- promoters(gtf_gene, upstream = 1000, downstream = 0)


#-----#
# CDS #
#-----#
gtf_cds <- gtf_all %>% filter(type == "CDS")

#-------#
# Exons #
#-------#
gtf_exons <- gtf_all %>% filter(type == "exon")

#---------#
# Introns #
#---------#
# Get txdb
gff_file <- argv$gff # "../../data/genomes/teleost_amphibia/Xlae/Xlae10.1/Sequence/Transcriptome/GCF_017654675.1_Xenopus_laevis_v10.1_genomic_shortnames.gff.gz"
gr <- BiocIO::import(gff_file)
txdb <- makeTxDbFromGRanges(gr)

# Get introns by transcripts
intronsByTranscript_gr <- intronsByTranscript(txdb)
# Combine all introns in a single GRanges object
intronsAll_gr <- unlist(intronsByTranscript_gr) %>% unique()


#------#
# UTRs #
#------#
# Get 5UTRs
fiveUTRsByTranscript_gr <- fiveUTRsByTranscript(txdb)
# Combine all 5UTRs in a single GRanges object
fiveUTRsAll_gr <- unlist(fiveUTRsByTranscript_gr) %>% unique()

# Get 3UTRs
threeUTRsByTranscript_gr <- threeUTRsByTranscript(txdb)
# Combine all 5UTRs in a single GRanges object
threeUTRsAll_gr <- unlist(threeUTRsByTranscript_gr) %>% unique()


#-------------#
# Centromeres #
#-------------#
centromeres_bed <- argv$centromeres # "../../data/genomes/teleost_amphibia/Xlae/Xlae10.1/Sequence/centromeres_Smith2021/GSE153058_xla_v10.2_cen_mod.bed" # From Smith2021
centromeres_data <- read.table(centromeres_bed, col.names = c("seqnames","start", "end"))
centromeres_gr <- centromeres_data %>% GRanges()

#-----------#
# Telomeres #
#-----------#
# From Bassham1998: telomeres range from 10 to over 50kb 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC121490/

# Load chromosome files
chr_sizes_file <- argv$chromsizes # "../../data/genomes/teleost_amphibia/Xlae/Xlae10.1/Xlae10_chr.chrom.sizes"
chr_sizes_data <- read.table(chr_sizes_file, col.names = c("seqnames", "size"))

# Get fiveTelomeres for 10kb and 50kb length
fiveTelomeres10kb_gr <- GRanges(seqnames = chr_sizes_data$seqnames,
                                ranges = IRanges(1, 10000))

fiveTelomeres50kb_gr <- GRanges(seqnames = chr_sizes_data$seqnames,
                                ranges = IRanges(1, 50000))

# Get threeTelomeres for 10kb and 50kb length
threeTelomeres10kb_gr <- GRanges(seqnames = chr_sizes_data$seqnames,
                                 ranges = IRanges(chr_sizes_data$size-10000, chr_sizes_data$size))

threeTelomeres50kb_gr <- GRanges(seqnames = chr_sizes_data$seqnames,
                                 ranges = IRanges(chr_sizes_data$size-50000, chr_sizes_data$size))


#==========================================#
# Create own annotation of genomic regions #
#==========================================#

message("\n######################")
message("Creating annotation...")
message("######################\n")


xenlae10_annotation <- list(
  version = "SaraCustom",
  promotersUp1000bp = promotersUp1000_gr,
  promoters200to500bp = promoters200to500_gr,
  promoters200bp = promoters200_gr,
  genes = gtf_gene,
  cds = gtf_cds,
  exons = gtf_exons,
  introns = intronsAll_gr,
  fiveUTRs = fiveUTRsAll_gr,
  threeUTRs = threeUTRsAll_gr,
  centromeres = centromeres_gr,
  fiveTelomeres10kb = fiveTelomeres10kb_gr,
  fiveTelomeres50kb = fiveTelomeres50kb_gr,
  threeTelomeres10kb = threeTelomeres10kb_gr,
  threeTelomeres50kb = threeTelomeres50kb_gr
)



#==============================================#
# ChIPQC to get enrichment for genomic regions #
#==============================================#

message("\n###################")
message("Computing ChIPQC...")
message("###################\n")

#--------#
# ChIPQC # 
#--------#
bam_file <- argv$bam

chipqc_out <- ChIPQCsample(bam_file, annotation =  xenlae10_annotation)


#-------------------------------#
# Enrichment of genomic regions # 
#-------------------------------#
# Retrieve genomic profile information in terms of relative enrichment over background genomic distribution
chipqc_regi <- regi(chipqc_out)

# Save output
write.table(chipqc_regi, file = file.path(argv$outdir, argv$outfile), sep = "\t", col.names = F, quote = F)
message("\nOutput saved in ", file.path(argv$outdir, argv$outfile))
