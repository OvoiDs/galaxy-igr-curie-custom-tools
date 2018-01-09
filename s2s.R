#!/usr/bin/Rscript
# Author: Dayris Thibault This is free and unencumbered software released into 
# the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or distribute 
# this software, either in source code form or as a compiled binary, for any 
# purpose, commercial or non-commercial, and by any means.
#
# In jurisdictions that recognize copyright laws, the author or authors of this 
# software dedicate any and all copyright interest in the software to the public 
# domain. We make this dedication for the benefit of the public at large and to 
# the detriment of our heirs and successors. We intend this dedication to be an 
# overt act of relinquishment in perpetuity of all present and future rights to 
# this software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>
suppressMessages(library("wasabi")); 
suppressMessages(library("argparse")); 
parser <- ArgumentParser(description="Prepares a Sleuth process from Salmon files" ); 
parser$add_argument(
  "-s", "--salmon_paths",
  metavar="PATH",
  type="character",
  nargs="+",
  help="Space separated paths to multiple Salmon zipped files (required)",
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
  "-o", "--output",
  metavar="PATH",
  type="character",
  help="Zip file output path (default: Salmon2Sleuth_archive)",
  default=file.path(getwd(), "Salmon2Sleuth_archive") ) 
args <- parser$parse_args(); 
args$salmon_paths <- as.character(args$salmon_paths); 
message("Argument parsed"); 
design <- data.frame(
  sample=args$names,
  path=sapply(args$salmon_paths, function(path) sub("^([^.]*).*", "\\1", path)),
  test=args$conditions ); 
print(design); 
unzip_and_prepare <- function(path) {
  unzipped_name <- paste0(sub("^([^.]*).*", "\\1", path), "_dir/output");
  suppressMessages(
    unzip(
      zipfile=path,
      unzip="unzip",
      exdir=paste0(sub("^([^.]*).*", "\\1", path), "_dir/")
    )
  );
  suppressMessages(prepare_fish_for_sleuth(unzipped_name));
  unzipped_name
}
tmp <- sapply(args$salmon_paths, function(zipped) unzip_and_prepare(zipped)); 
zip(paste0(args$output, ".zip"),
  files=tmp,
  zip="zip" ); 

message("Files converted"); 
write.table(
  design,
  "sleuth_design.tsv",
  quote=FALSE,
  sep="\t",
  row.names=FALSE ) 
message("Process is now over");
