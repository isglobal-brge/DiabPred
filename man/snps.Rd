\name{Model SNPs}
\docType{data}
\title{
sNP data frame
}
\description{
SNP data frame consisting of 8573 rows and 6 columns. The SNPs figuring in this file were used to create the model that's used in this package for predictions.
}
\usage{read.delim(system.file("extdata", "snps.txt", package = DiabPred), sep = "")}
\format{
  A data frame with 8573 rows and 6 columns.
  \describe{
  \item \code{chr}: Chromosome where the SNP is located
  \item \code{snp_id}: SNP reference ID
  \item \code{chr_start}: Start position of the SNP in the chromosome
  \item \code{chr_end}: End position of the SNP in the chromosome
  \item \code{ref}: Reference allele of the SNP
  \item \code{alt}: Alternative allele of the SNP
  }
}
\details{
These SNPs were annotated using the version X
}
\source{
\url{http://www.type2diabetesgenetics.org/}
}
\keyword{datasets}
