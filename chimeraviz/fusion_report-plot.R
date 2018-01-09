#!/galaxy-central/database/dependencies/_conda/envs/__r-base@3.4.1/bin/Rscript
library(optparse)
library(chimeraviz)
library(Rsamtools)

option_list <- list(
make_option(c("-i", "--input"), dest="input"),
make_option(c("-o", "--output"), dest="output"),
make_option(c("-b", "--bam"), dest="ibam", default=NULL),
make_option(c("-e", "--edb"), dest="edb", default="/local_tools/chimeraviz/Homo_sapiens.GRCh37.74.sqlite"),
make_option(c("-d", "--do"), dest="d"))
# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                                           option_list=option_list)
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
input <- opt$input
output <- opt$output
ibam <- paste0(workDir, '/', opt$ibam)
edb <- opt$edb
id <- opt$d
print(ibam)

# index our bam file
indexBam(ibam)

fusions <- importFusioncatcher(input, "hg19", 500)
edbdb <- ensembldb::EnsDb(edb)

specificFusion <- getFusionById(fusions, id)

png('./plotFusion.png', height=800, width=600)
plotFusion( fusion = specificFusion, bamfile = ibam, edb = edbdb, nonUCSC = FALSE)
dev.off()
png('./transcriptPlot.png', height=800, width=600)
plotTranscripts(fusion = specificFusion, edb = edbdb, bamfile = ibam, nonUCSC = FALSE, ylim=c(0,500))
dev.off()
#png('./fusionTranscriptPlot.png', height=800, width=600)
#plotFusionTranscript(fusion = specificFusion, edb = edbdb, bamfile = ibam)
#dev.off()
png('./FusionTranscriptGraph.png', height=800, width=600)
plotFusionTranscriptsGraph(fusion = specificFusion, edb = edbdb)
dev.off()
