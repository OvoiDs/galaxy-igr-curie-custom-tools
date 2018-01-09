#!/usr/bin/Rscript
# Author: Thibault Dayris

suppressMessages(library("sleuth")); suppressMessages(library("argparse")); parser <- ArgumentParser(
  description="Perform your basic sleuth analysis with this script" ); parser$add_argument(
  "design",
  metavar="PATH",
  type="character",
  help="Path to experimental meta design file" ); parser$add_argument(
  "-w", "--workdir",
  metavar="PATH",
  type="character",
  help="Working directory for the R session (default: cwd)",
  default=getwd() ); parser$add_argument(
  "-t", "--tr2gene",
  metavar="PATH",
  type="character",
  help="Use a transcript to gene table (default: None)",
  default=NULL ); parser$add_argument(
  "-o", "--output",
  metavar="BASE",
  type="character",
  help="Name of output file (default: Sleuth_Results)",
  default="Sleuth_Results" ); parser$add_argument(
  "-p", "--threads",
  metavar="N",
  type="integer",
  help="Maximum number of threads to be used in the analysis (default: 1)",
  default=1 ); parser$add_argument(
  "-e", "--gene_analysis",
  action="store_true",
  default=FALSE,
  help="Perform a gene analysis" ); opt <- parser$parse_args(); mc.cores=opt$threads; opt$design <- as.character(opt$design); opt$tr2gene <- as.character(opt$tr2gene); 
opt$output_table <- file.path(opt$workdir, paste0(opt$output, ".tsv")); opt$output_rds <- file.path(opt$workdir, paste0(opt$output, ".rds")); opt$output_genes <- 
file.path(opt$workdir, paste0(opt$output, ".genes.tsv")); opt$output_exp <- file.path(opt$workdir, paste0(opt$output, ".expression.tsv")); opt$args_out <- file.path(opt$workdir, 
paste0(opt$output, ".arguments.rds")); if (is.null(opt$tr2gene)) {
  opt$t2g <- read.table(opt$tr2gene, header = TRUE, sep = "\t");
  colnames(opt$t2g) <- c("ens_gene", "target_id", "ext_gene");
} else {
  opt$t2g <- NULL;
}
# Set tmp directory to avoid /tmp issues
TMP <- file.path(opt$workdir, "tmp"); TEMP <- file.path(opt$workdir, "tmp"); TMPDIR <- file.path(opt$workdir, "tmp");
# print(opt);
saveRDS(opt, file = opt$args_out); meta <- read.table(opt$design, header = TRUE, sep = "\t"); meta$path <- as.character(meta$path); print("Meta data loaded"); if 
(!is.null(opt$genes_analysis)) {
  if (is.null(opt$t2g)) {
    print("Gene Aggregation was not possible without transcript to gene table");
    exit(1);
  }
  so <- sleuth_prep(
      meta,
      ~test,
      target_mapping = opt$t2g,
      aggregation_column = "ens_gene"
  );
} else {
  so <- sleuth_prep(
      meta,
      ~test,
      target_mapping = opt$t2g
  );
}
print("Sleuth object prepared");
# Formula given in sleuth_prep, leaving default arguments
so <- sleuth_fit(so); print("Sleuth model fitted");
# The intercept only model
so <- sleuth_fit(so, ~1, "reduced"); print("Sleuth intercept fitted"); tmp <- sort(meta$test); ref_cond <- paste0("test", tmp[length(meta$test)]); so <- sleuth_wt(so, which_beta = 

ref_cond); print("Wald Test over"); write.table(
    sleuth_results(so, ref_cond),
    file = opt$output_table,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t" ) print("Sleuth Result table saved");

# Geathering transcripts expression
quant <- sleuth::kallisto_table(so);

# Reshaping the dataframe
reduced_quant <- quant[c("target_id", "sample", "tpm")]; re_ordered <- tidyr::spread(reduced_quant, sample, tpm); rownames(re_ordered) <- re_ordered$target_id; re_ordered$target_id 
<- NULL; write.table(re_ordered, file = opt$output_exp,
            quote = FALSE, row.names = TRUE, sep = "\t"); print("Transcripts expression table saved");

# Saving results
saveRDS(so, file = opt$output_rds); print("RDS saved"); write.table(
    sleuth_results(so, ref_cond),
    file = opt$output_table,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t" ); print("Sleuth Result table saved");
