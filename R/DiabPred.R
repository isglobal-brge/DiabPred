#' Prediction of type 2 diabetes based on genomic data with deep learning algorithms
#'
#' The function described here makes use of a deep learning model trained using the h2o package to predict type 2 diabetes based on selected SNPs.
#'
#'
#'
#' @title DiabPred
#' @param genofile An object of class \code{SNPGDSFileClass} loaded into the environment with \code{snpgdsOpen} from the \strong{SNPRelate} package.
#'   
#' @export DiabPred
#' 
#' @return A data frame with individuals in rows and 3 colums: 
#' \itemize{ 
#'   \item \code{prediction}: First column with the actual phenotypic output indicating wether the individual is diabetic or non-diabetic.
#'   \item \code{non.diabetic}: Probability of the individual being non-diabetic
#'   \item \code{diabetic}: Probability of the individual being diabetic.
#'   }

DiabPred <- function(genofile) {
  snps <- as.character(read.delim(system.file("extdata", "snps.txt", package = "DiabPred"), sep = "")$snp_id)
  gds <- snpgdsGetGeno(genofile, snp.id = snps, with.id = TRUE)
  matrix <- gds$genotype
  colnames(matrix) <- gds$snp.id
  rownames(matrix) <- gds$sample.id
  predictions <- h2o.mojo_predict_df(frame = matrix, 
                                     mojo_zip_path=system.file("extdata", "model.zip", package = "DiabPred"), 
                                     genmodel_jar_path=system.file("extdata", "h2o-genmodel.jar", package = "DiabPred"))
  rownames(predictions) <- gds$sample.id
  return(predictions)
}




