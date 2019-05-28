#' List of SNPs used to build this package's model
#'
#' @title SNP info
#' @description This dataset consists on the SNPs that were used to build the deep learning model used for prediction. The genome reference used for the annotation of these SNPs was hg19.
#'   
#' @format A data.frame consisting of 8573 rows and 6 columns
#' 
#' \describe{ 
#'   \item{chr}{Chromosome where the SNP is located}
#'   \item{snp_id}{SNP reference ID}
#'   \item{chr_start}{Start position of the SNP in the chromosome}
#'   \item{chr_end}{End position of the SNP in the chromosome}
#'   \item{ref}{Reference allele of the SNP}
#'   \item{alt}{Alternative allele of the SNP}
#'   }
#'   
#'   
#' @source \url{http://www.type2diabetesgenetics.org/}
"snps"