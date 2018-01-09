#!/usr/bin/Rscript
# Author: Thibault Dayris

# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>

# Libraries
suppressMessages(library("sleuth"));
suppressMessages(library("argparse"));
suppressMessages(library("rmarkdown"));
suppressMessages(library("DT"));

sink(NULL, type = "message");
# Additional functions

panel.cor <- function(x, y, digits=2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1));
  r <- abs(cor(x, y));
  txt <- format(c(r, 0.123456789), digits=digits)[1];
  test <- cor.test(x,y);
  Signif <- ifelse(
    round(test$p.value,3)<0.001,
    "p<0.001",
    paste("p=",round(test$p.value,3))
  );
  text(0.5, 0.25, paste("r=",txt));
  text(.5, .75, Signif);
}

# Command line parser
parser <- ArgumentParser(
  description = "Perform your basic sleuth analysis with this script"
);

parser$add_argument(
  "design",
  metavar = "PATH",
  type = "character",
  help = "Path to experimental meta design file"
);

parser$add_argument(
  "report",
  metavar = "PATH",
  type = "character",
  help = "Path to report markdown file"
);

parser$add_argument(
  "counts",
  metavar = "PATH",
  type = "character",
  help = "Path to counts.zip file"
);

parser$add_argument(
  "-w", "--workdir",
  metavar = "PATH",
  type = "character",
  help = "Working directory for the R session (default: cwd)",
  default = getwd()
);

parser$add_argument(
  "-t", "--tr2gene",
  metavar = "PATH",
  type = "character",
  help = "Use a transcript to gene table (default: None)",
  default = NULL
);

parser$add_argument(
  "-o", "--output",
  metavar = "BASE",
  type = "character",
  help = "Name of output file (default: Sleuth_Results)",
  default = "Sleuth_Results"
);

parser$add_argument(
  "-p", "--threads",
  metavar = "N",
  type = "integer",
  help = "Maximum number of threads to be used in the analysis (default: 1)",
  default = 1
);

parser$add_argument(
  "-e", "--gene_analysis",
  action = "store_true",
  default = FALSE,
  help = "Perform a gene analysis (default: false)"
);

# Internal options definition
opt <- parser$parse_args();
mc.cores <- opt$threads;
options(mc.cores = opt$threads);

opt$design <- as.character(opt$design);
opt$tr2gene <- as.character(opt$tr2gene);
opt$counts <- as.character(opt$counts);
opt$gene_analysis <- as.logical(opt$gene_analysis)

opt$output_table <- file.path(
  opt$workdir,
  paste0(opt$output, ".tsv")
);

opt$output_rds <- file.path(
  opt$workdir,
  paste0(opt$output, ".rds")
);

opt$output_exp <- file.path(
  opt$workdir,
  paste0(opt$output, ".expression.tsv")
);

opt$args_out <- file.path(
  opt$workdir,
  paste0(opt$output, ".arguments.rds")
);

opt$output_html <- file.path(
  opt$workdir,
  paste0(opt$output, ".report.html")
);

# Transcript to gene table loading
if (!is.null(opt$tr2gene)) {
  opt$t2g <- read.table(opt$tr2gene, header = FALSE, sep = "\t", quote="");
  colnames(opt$t2g) <- c("ens_gene", "target_id", "ext_gene");
} else {
  opt$t2g <- NULL;
}

# Set tmp directory to avoid /tmp issues
TMP <- file.path(opt$workdir, "tmp");
TEMP <- file.path(opt$workdir, "tmp");
TMPDIR <- file.path(opt$workdir, "tmp");

# Save arguments
saveRDS(opt, file = opt$args_out);

# Load experimental design
meta <- read.table(opt$design, header = TRUE, sep = "\t");
meta$path <- as.character(meta$path);

unzip(
  opt$counts,
  unzip="unzip"
);
message("Metadata loaded");

# Prepare Sleuth run
# target_mapping is NULL by default in Sleuth
if (isTRUE(opt$gene_analysis)) {
  message("Differential Gene Expression Analysis")
  if (is.null(opt$t2g)) {
    message("Gene Aggregation was not possible: no transcript to gene table");
    q(status = 1);
  }
  so <- sleuth_prep(
      meta,
      ~test,
      target_mapping = opt$t2g,
      aggregation_column = "ens_gene",
      num_cores=opt$threads
  );
} else {
  message("Differential Transcripts Expression Analysis")
  so <- sleuth_prep(
      meta,
      ~test,
      target_mapping = opt$t2g,
      num_cores=opt$threads
  );
}
message("Sleuth object prepared");

# Formula given in sleuth_prep, leaving default arguments
so <- sleuth_fit(so);
message("Sleuth model fitted");

# The intercept only model
so <- sleuth_fit(so, ~1, "reduced");
message("Sleuth intercept fitted");

# Perform Wald test to get diff exp
tmp <- sort(meta$test);
ref_cond <- paste0("test", tmp[length(meta$test)]);
so <- sleuth_wt(so, which_beta = ref_cond);
message("Wald Test over");

# Save results
write.table(
    sleuth_results(so, ref_cond),
    file = opt$output_table,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
);

message("Sleuth Result table saved");

# Geathering transcripts expression
quant <- sleuth::kallisto_table(so);

# Reshaping the dataframe
reduced_quant <- quant[c("target_id", "sample", "tpm")];
re_ordered <- tidyr::spread(reduced_quant, sample, tpm);

# Defining target_id as index column
rownames(re_ordered) <- re_ordered$target_id;
re_ordered$target_id <- NULL;

# Save expression table
write.table(
  re_ordered,
  file = opt$output_exp,
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
);
message("Transcripts expression table saved");

# Save results table
write.table(
    sleuth_results(so, ref_cond),
    file = opt$output_table,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
);
message("Sleuth Result table saved");

# Saving RDS, as pairwise scatter plot may be very large for huge
# experimental designs.
sleuth::sleuth_save(so, file = opt$output_rds);
message("RDS saved");

render(
  opt$report,
  output_file = opt$output_html,
  params = list(sleuth = so, ref_cond = ref_cond),
  quiet=TRUE
);
message("Report created");
message("Process is now over");
q(status = 0);
