#' List of SNPs used to build this package's model
#'
#' @title SNP info
#' @description This dataset consists on the SNPs that were used to build the deep learning model used for prediction.
#'   
#' @format A data frame with 8573 rows and 6 columns
#' 
#' \describe{ 
#'   \item \code{chr}{Chromosome where the SNP is located}
#'   \item \code{snp_id}{SNP reference ID}
#'   \item \code{chr_start}{Start position of the SNP in the chromosome}
#'   \item \code{chr_end}{End position of the SNP in the chromosome}
#'   \item \code{ref}{Reference allele of the SNP}
#'   \item \code{alt}{Alternative allele of the SNP}
#'   }
#'   
#'   
#' @source \url{http://www.type2diabetesgenetics.org/}
"snps"