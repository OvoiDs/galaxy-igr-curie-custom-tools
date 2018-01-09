#!/usr/bin/Rscript
# Author: Dayris Thibault This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, 
# commercial or non-commercial, and by any means.
#
# In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public 
# domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt 
# act of relinquishment in perpetuity of all present and future rights to this software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
# OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>
suppressMessages(library("tximport")); suppressMessages(library("argparse")); parser <- ArgumentParser(
  description="Prepares a SArtools process from Salmon files" ); parser$add_argument(
  "-s", "--salmon_paths",
  metavar="PATH",
  type="character",
  nargs="+",
  help="Space separated paths to multiple Salmon quant.sf files (required)",
  default=NULL ); 
parser$add_argument(
  "-c", "--conditions",
  metavar="COND",
  type="character",
  nargs="+",
  help="Corresponding space separated list of conditions (required)",
  default=NULL ); 
parser$add_argument(
  "-n", "--names",
  metavar="NAME",
  type="character",
  nargs="+",
  help="Corresponding space separated list of names (required)",
  default=NULL );
parser$add_argument(
  "-t", "--transcript_2_gene",
  metavar="PATH",
  type="character",
  help="Path to transcript to gene table",
  default=NULL ); 
args <- parser$parse_args(); args$salmon_paths <- sapply(
  args$salmon_paths,
  function(path) as.character(path) ); 
args$transcript_2_gene <- as.character(args$transcript_2_gene); message("Command line parsed"); 
save_columns <- function(dataframe) {
  invisible(
    sapply(
      colnames(dataframe),
      function(idx) write.table(
        dataframe[, idx],
        file=paste0(
          "Genes_Counts_",
          idx,
          ".txt"
        ),
        sep="\t",
        quote=FALSE,
        col.names=FALSE
      )
    )
  );
}
message("banana");
round_df <- function(dataframe, digits) {
    numeric_columns <- sapply(dataframe, mode) == 'numeric'
    dataframe[numeric_columns] <- round(dataframe[numeric_columns], digits)
    dataframe
}
tr2gene <- read.table(args$transcript_2_gene); 
message("banana1.5");
aggregation <- tximport(
  args$salmon_paths,
  type="salmon",
  tx2gene=data.frame(tx=tr2gene$V2, gene=tr2gene$V3) ); 
message("banana2")
agcounts <- aggregation$counts;
colnames(agcounts) <- args$names; 
print(head(agcounts)); 
message("Aggregation performed");
# Never use normalized data with DESeq2 !
save_columns(round_df(agcounts)); counts_files <- sapply(
  args$names,
  function(name) paste0("Genes_Counts_", name, ".txt") ); 
zip(
  zipfile=file.path(getwd(), "Genes_Counts.zip"),
  files=counts_files,
  zip="zip" ); 
message("Counts saved"); design <- data.frame(
  label <- args$names,
  files <- counts_files,
  group <- args$conditions ); 
colnames(design) <- c("label", "files", "group"); 
write.table(
  design,
  file="SARTools_design.txt",
  sep="\t",
  quote=FALSE,
  row.names=FALSE ); 
message("Design file saved");
