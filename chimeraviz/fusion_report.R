#!/galaxy-central/database/dependencies/_conda/envs/__r-base@3.4.1/bin/Rscript
library(optparse)
library(chimeraviz)
library(Rsamtools)

option_list <- list(
make_option(c("-i", "--input"), dest="input"),
make_option(c("-o", "--output"), dest="output"),
make_option(c("-e", "--edb"), dest="edb", default="/local_tools/chimeraviz/Homo_sapiens.GRCh37.74.sqlite"))

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                                           option_list=option_list)
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
input <- opt$input
output <- opt$output
edb <- opt$edb

fusions <- importFusioncatcher(input, "hg19", 500)

createFusionReport(fusions, output)
