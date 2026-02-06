#!/usr/bin/env Rscript

library(dplyr)
library(SummarizedExperiment)
library(optparse)
library(utils)


#' Force a vector of values into standard normal distribution
#'
#' @param x numeric vector with arbitrary distribution
#'
#' @return Vector with a standard normal distribution
#' @export
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}

quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

#' Perform inverse normal transformation on the columns of a matrix
#'
#' @param x matrix
#'
#' @return matrix
#' @export
quantileNormaliseCols <- function(matrix,...){
  quantileNormaliseMatrix(matrix, ...)
}

#' Perform inverse normal transformation on the rows of a matrix
#'
#' @param x matrix
#'
#' @return matrix
#' @export
quantileNormaliseRows <- function(matrix,...){
  t(quantileNormaliseMatrix(t(matrix), ...))
}

normaliseSE_quantile <- function(se, assay_name = "usage"){

  #Extract assays
  assay_list = SummarizedExperiment::assays(se)
  assay_matrix = assay_list[[assay_name]]

  #Quantile normalise
  print(dim(assay_matrix))
  qnorm = quantileNormaliseRows(assay_matrix)
  assay_list[["qnorm"]] = qnorm

  #Make ab update se object
  se = SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    colData = SummarizedExperiment::colData(se),
    rowData = SummarizedExperiment::rowData(se))
  return(se)
}

option_list <- list(
  optparse::make_option(c("-p", "--phenotype_matrix"), type="character", default=NULL,
              help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("-m", "--min-exp"), type="integer", default = NULL,
              help="Minimal number of not-NA values per phenotype", metavar = "numeric"),
  optparse::make_option(c("-o", "--output_file_prefix"), type="character", default="phenotype_matrix_norm",
              help="Prefix for output files", metavar = "character")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

expression_matrix_path <- opt$p
output_prefix <- opt$o

message('Reading the phenotype matrix...')
phm <- read.csv(expression_matrix_path, sep='\t')
message(paste0('# phenotypes before min exp filtering: ', nrow(phm)))

if (is.null(opt$m)) {
  sample_count <- ncol(phm) - 1
  minexp <- ceiling(sample_count * 0.01)
  message("min_exp not provided → auto-calculated: ", minexp)
} else {
  minexp <- opt$m
  message("Using provided min_exp: ", minexp)
}

phm <- phm[rowSums(!is.na(phm[ , -1, drop = FALSE])) > minexp, ]
message(paste0('# phenotypes after min exp filtering: ', nrow(phm)))

sample_columns <-  grep('phenotype_id', colnames(phm), value = TRUE, invert = TRUE)

# Create the assay matrix
usage_matrix <- phm[, !colnames(phm) %in% 'phenotype_id', drop = FALSE] |> as.matrix()
rownames(usage_matrix) <- phm$phenotype_id

# Row metadata
row_data <- phm[, c('phenotype_id'), drop = FALSE]
rownames(row_data) <- row_data$phenotype_id

col_data <- DataFrame(sample_id = sample_columns)
rownames(col_data) <- sample_columns

# Create SummarizedExperiment with assay named "usage"
se <- SummarizedExperiment(
  assays = list(usage = usage_matrix),
  rowData = row_data,
  colData = col_data
)
message("Normalising the matrix...")
se_n <- normaliseSE_quantile(se)

qnorm_matrix <- assays(se_n)[["qnorm"]]

qnorm_df <- as.data.frame(qnorm_matrix)

qnorm_df <- tibble::tibble(phenotype_id = rownames(qnorm_df), qnorm_df)
tsv_file <- paste0(output_prefix, ".tsv")
gz_file  <- paste0(output_prefix, ".tsv.gz")
utils::write.table(qnorm_df, file = tsv_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
R.utils::gzip(tsv_file, destname = gz_file, overwrite = TRUE)
