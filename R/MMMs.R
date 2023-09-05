
#' check.missing 
#' 
#' @description 
#' Check missing variables for prediction. The function checks the missing values of
#' the new data with respect to the training data, including the missing background 
#' values and missing model values. Missing background values are not directly involved 
#' in the fitted models, but they may contribute to the calculation of pathway scores, 
#' e.g., gsva. 
#' 
#' @param X new data
#' @param pathways one of the "hall","biocarta","pid","cgp","wiki","tft.legency",
#' "tft.gtrd", "cm", "gobp", "gocc", "gomf", "oncogene", "immune", "vax" and "genes". 
#' If pathways = "genes", the gene based model is selected.
#' @param method one of "ssgsea", "gsva", "zscore", "grp", if models are constructed 
#' based on pathways. 
#' @param allgenes all genes in the training data
#' 
#' @examples
#' NewX <- Xtrain[-c(1:100), 1:10]
#' check.missing(NewX, pathways = "vax", method = "grp", allgenes = rownames(Xtrain)) 
#' 
#' @export check.missing

check.missing <- function(X, pathways, method, allgenes){
  
  #- check missing of background genes
  miss.bg <- setdiff(allgenes, rownames(X))
  
  if (length(miss.bg) > 0) {
    cat("Missed background genes:", paste0(miss.bg, collapse = ","), "\n\n")
  } 
  
  # pathways
  if (pathways == "genes") {
    mod <- modsGenes
    miss.md <- setdiff(names(mod), rownames(X))
  } else if (method == "grp") {
    mod <- modsPath[[pathways]][["grp"]][["beta"]]
    miss.md <- setdiff(names(mod), rownames(X))
  } else {
    mod <- modsPath[[pathways]][[method]]
    sets <- pathSets[[pathways]]
    vars <- names(mod)
    sets <- sets[vars]
    miss.md <- unique(unlist(lapply(sets, setdiff, y = rownames(X))))
  }
  
  if (length(miss.md) > 0) {
    cat("Missed model genes:", paste0(miss.md, collapse = ","), "\n\n")
  } 
  
  return(list(model.genes = miss.md, 
              model.background = miss.bg))
}

#' @export genX
#' @import GSVA
#' 

genX <- function(X, pathways, method){
  
  # method, "ssgsea", "gsva", "zscore", "grp", "gene"
  
  if (pathways == "genes"){
    mod <- modsGenes
    return(X[names(mod),])
  } else if (method == "grp") {
    mod <- modsPath[[pathways]][["grp"]][["beta"]]
    return(X[names(mod),, drop = F])
  } else {
    mod  <- modsPath[[pathways]][[method]]
    vars <- names(mod)
    
    sets <- pathSets[[pathways]][vars]
    newX <- GSVA::gsva(X, sets, min.sz = 1, method = method, ssgsea.norm = F)
    return(newX) 
  }
  
}

#' miss.impute 
#' @description 
#' Imputed missing values using K-nearest neighbors based on the training data (GSE136324).
#'  
#' @param NewX new data
#' @param training training data used for model building, GSE136324 
#' @examples 
#' NewX <- Xtrain[-c(1:100), 1:10]
#' miss.impute(NewX, Xtrain)
#' 
#' @import impute
#' @export miss.impute

miss.impute <- function(NewX, Xtrain) {
  
  #- if individual patients in vector, transform to Matrix 
  if (is.vector(NewX))
    NewX <- matrix(NewX, ncol = 1)
  
  NewX <- as.matrix(NewX)
  Xtrain <- as.matrix(Xtrain)
  #- For valid, select all features including in training features
  sel  <- intersect(rownames(NewX), rownames(Xtrain))
  NewX <- NewX[sel, , drop = F]
  
  #- tranform valid to training data structure
  ind <- match(rownames(NewX), rownames(Xtrain))
  mat <- matrix(NA, nrow = nrow(Xtrain), ncol = ncol(NewX), 
                dimnames = list(rownames(Xtrain),
                                colnames(NewX)))
  mat[ind, ] <- NewX
  
  for(i in seq_len(ncol(NewX))){
    
    cat("Processing:", colnames(NewX)[i], "\n")
    sink("nul")
    out <- impute::impute.knn(cbind(Xtrain, ref = mat[,i]), k = 10)
    sink()
    mat[,i] <- out$data[,"ref"]
  }
  mat
}

#' predictMM 
#' @description 
#' 
#' The function provided a set of prognostic models constructed either based on
#' genes or pathways. For the pathway models, two strategies were implemented to construct predictive 
#' models: pathway score method and group lasso using pathway information. 
#' We considered three methods for pathway score calculation (ssGSEA, GSVA, and z-scores) and 
#' 14 data sources providing pathway information. 
#' We recommended a group lasso model based on the vax pathways. 
#'  
#' @param X variables of new data
#' @param pathways one of 14 sources of pathways: "hall","biocarta","pid","cgp","wiki","tft.legency",
#' "tft.gtrd", "cm", "gobp", "gocc", "gomf", "oncogene", "immune","vax". 
#' If pathways = "genes", a gene based model is selected.
#' @param method one of "ssgsea", "gsva", "zscore", "grp", if models are constructed 
#' based on pathways. 
#' @param miss.impute whether missing data should be imputed.
#' @param type the type of predicted value. Choices are the linear predictor ("lp"), 
#' the risk score exp(lp) ("risk"), the expected number of events given the 
#' covariates and follow-up time ("expected"), and the terms of the linear predictor
#' ("terms"). The survival probability for a subject is equal to exp(-expected).
#' see predict.coxph in survival package. 
#' 
#' @examples 
#' NewX <- Xtrain[-c(1:100), 1:10]
#' predictMM(NewX, pathways = "vax", method = "grp", type = "lp")
#' 
#' @import survival
#' @import impute
#' @import GSVA
#' @export predictMM

predictMM <- function(X, 
                      pathways = "immune",
                      method = "ssgsea",
                      miss.impute = T, 
                      type = c("lp", "risk", "expected", "terms", "survival")){
  
  pathways  <- match.arg(pathways, 
                         c("hall","biocarta","pid","cgp","wiki","tft.legency",
                           "tft.gtrd", "cm", "gobp","gocc","gomf","oncogene",
                           "immune","vax", "genes"))
  
  method <- match.arg(method, c("gsva", "ssgsea", "zscore", "grp"))
  
  # check for missing
  miss <- check.missing(X, pathways, method, allgenes = rownames(Xtrain))
  
  # whether missing check is needed
  miss <- (length(miss$model.genes) > 0) || (length(miss$model.background) > 0 & pathways != "genes")
  
  # impute missing data
  if (miss){
    if (miss.impute) {
      X <- miss.impute(X, Xtrain)
    } else{
      stop("Missing data need to be imputed")
    }
  }
  
  #- generate data
  
  if (pathways == "genes"){
    mod <- modsGenes
  } else if (method == "grp"){
    mod <- modsPath[[pathways]][["grp"]][["beta"]]
  } else {
    mod <- modsPath[[pathways]][[method]]
  }
  
  if (length(mod) == 0) {
    
    print("Model with no variable")
    return("")
    
  } else {
    
    Z <- genX(X, pathways = pathways, method = method)
    Z <- Z[match(names(mod), rownames(Z)), ,drop = F]
    
    #- generate Z train
    Ztrain <- genX(Xtrain, pathways = pathways, method = method)
    Ztrain <- Ztrain[match(names(mod), rownames(Ztrain)), ,drop = F]
    
    m <- survival::coxph(ytrain ~ ., data = as.data.frame(t(Ztrain)), 
                         init = mod, iter.max = 0)
    
    #- prediction 
    out <- predict(m, as.data.frame(t(Z)), type = type)
    return(out)
  }
}





