#' Detection of individuals with Extreme Downregulation of chromosome Y (EDY)
#' from transcriptomic data
#'
#'
#' Escribir los detalles ...
#'
#'
#'
#' @title DiabPred
#' @param x An ExpressionSet or RangedSummarizedExperiment from microarray or
#'   RNAseq experiments.
#' @param gender.var A string indicating the name of the column in the table of
#'   phenotype (accessed by \code{pData} or \code{colData}) that contains the
#'   gender of the \strong{individuals}.
#'   
#' @export DiabPred
#' 
#' @return A data frame with individuals in rows and 3 colums: 
#' \itemize{ 
#'   \item \code{prediction}: probabi .... .
#'   \item \code{non.diabetic}: probabi .... .
#'   \item \code{diabetic}: probabi ....  . 
#'   }

getEDY <- function(x, gender.var, male.key, gene.key, coef=1.2, 
                   log = TRUE, group.var, control.key, experiment.type, ...){
  
  object.type <- tolower(class(x)[1])
  experiment.type <- tolower(experiment.type)
  
  if (experiment.type!="microarray" && experiment.type!="rnaseq"){
    stop("Invalid experiment.type. Allowed types are 'microarray' or 'RNAseq'")
  }
  
  if (object.type == "expressionset"){
    
      #Filter males in the expression set
      if (!missing(gender.var)) {
        ii <- which(varLabels(x)==gender.var)
        x <- x[,pData(x)[,ii]==male.key]
        x <- x[,!is.na(pData(x)[,ii])]
      } else {warning("male.key not specified. Performing analysis with the whole dataset")}
      
      #Add column with the hgnc symbol information to fData
      if (gene.key != "hgnc_symbol"){
        fData(x)$id.feature <- featureNames(x)
        }
      
      annot.expr <- fData(x)
      
      #Select from genes in gene.expr those that we know the hgnc symbol
      gene.expr <- exprs(x)[rownames(exprs(x))%in%annot.expr$id.feature,]
      #Replace gene ID for hgnc symbol
      rownames(gene.expr) <- annot.expr[, gene.key]
      #Select those genes that belong to chrY
      exprY <- gene.expr[annot.expr[, gene.key]%in%EDY::chrY.genes$hgnc_symbol,]
      exprY <- exprY[complete.cases(exprY),]
      #Select those genes that belong to the rest of the genome
      exprRef <- gene.expr[annot.expr[, gene.key]%in%EDY::autosomal.genes$hgnc_symbol,]
      exprRef <- exprRef[complete.cases(exprRef),]
      
  }
  else if (object.type == "rangedsummarizedexperiment"){
    
    #Filter males in the expression set
    if (!missing(gender.var)){
      ii <- which(names(colData(x)) == gender.var)
      x <- x[, !is.na(colData(x)[, ii])]
      x <- x[, colData(x)[, ii] == male.key]
    } else {warning("male.key not specified. Performing analysis with the whole dataset")}
    
    #Add column with the hgnc symbol information to fData
    if (gene.key != "hgnc_symbol"){
      rowData(x)$id.feature <- rownames(assay(x))
    }
    
    annot.expr <- data.frame(rowData(x))
    
    #Select from genes in gene.expr those that we know the hgnc symbol
    gene.expr <- assay(x)[rownames(assay(x))%in%annot.expr$id.feature,]
    #Replace gene ID for hgnc symbol
    rownames(gene.expr) <- annot.expr[, gene.key]
    #Select those genes that belong to chrY
    exprY <- gene.expr[annot.expr[, gene.key]%in%EDY::chrY.genes$hgnc_symbol,]
    exprY <- exprY[complete.cases(exprY),]
    #Select those genes that belong to the rest of the genome
    exprRef <- gene.expr[annot.expr[, gene.key]%in%EDY::autosomal.genes$hgnc_symbol,]
    exprRef <- exprRef[complete.cases(exprRef),]
    
  }
  
  #Apply EDY formulae: 
  if (experiment.type=="rnaseq"){
    
    if (log) {
      Ry <- sweep((exprY), 2, FUN="-", 
                  apply((exprRef), 2, mean))
      }else{ 
        Ry <- sweep(log2(exprY+1), 2, FUN="-", 
                  apply(log2(exprRef+1), 2, mean))}
    } else if (experiment.type == "microarray"){
        if (log) {
        Ry <- sweep(exprY, 2, FUN="-", 
                  apply(exprRef, 2, mean))
      }else{ 
        Ry <- sweep(log2(exprY), 2, FUN="-", 
                  apply(log2(exprRef), 2, mean))
      }
    }
  
  
  EDYcontinuous <- apply(Ry, 2, mean)
  
  if (!missing(group.var)&&!missing(control.key)){
    if (object.type == "expressionset"){
    controls <- EDYcontinuous[pData(x)[,group.var]==control.key]
    } else if (object.type == "rangedsummarizedexperiment"){
    controls <- EDYcontinuous[colData(x)[,group.var]==control.key]
    }
  } else {
    controls <- EDYcontinuous
    warning("No control group specified")
  }
  
  thresh <- median(controls, na.rm=TRUE) - coef*IQR(controls, na.rm=TRUE)
  EDY <- cut(EDYcontinuous, c(-Inf, thresh, Inf), 
             labels = c("Yes", "No"))
  EDY <- relevel(EDY, 2)
  
  #output
  if (object.type == "expressionset"){
    pData(x)$EDY <- EDY
    names(EDY) <- names(EDYcontinuous)
    ans <- list(EDY=EDY, EDYcontinuous=EDYcontinuous, 
                threshold=thresh, eSet=x)
  } else if (object.type == "rangedsummarizedexperiment"){
    colData(x)$EDY <- EDY
    names(EDY) <- names(EDYcontinuous)
    ans <- list(EDY=EDY, EDYcontinuous=EDYcontinuous, 
                threshold=thresh, RangedSummarizedExperiment=x)
  }
  class(ans) <- "EDY"
  ans
}




