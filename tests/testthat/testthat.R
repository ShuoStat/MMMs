
# expected Results
# 
# NewX <- Xtrain[-c(1:100), 1:2]
# 
# pathList <- list()
# for(pathways in c(names(pathSets), "genes")){
#   print(pathways)
#   methodList <- list()
#   for(method in c("ssgsea", "gsva", "grp", "zscore")){
#     methodList[[method]]<- predictMM(NewX, pathways = pathways, method = method, type = "lp")
#   }
#   pathList[[pathways]] <- methodList
# }
# 
# expectedPathList <- pathList
# 
# save(list = "expectedPathList", file = "../expectedResults/expectedResults.RData")


#- 


testthat::test_that(
  
  "Test predictMM funciton",
  {
    NewX <- Xtrain[-c(1:100), 1:2]
    pathList <- list()
    for(pathways in c(names(pathSets), "genes")){
      print(pathways)
      methodList <- list()
      for(method in c("ssgsea", "gsva", "grp", "zscore")){
        methodList[[method]]<- predictMM(NewX, pathways = pathways, method = method, type = "lp")
      }
      
      pathList[[pathways]] <- methodList
    }
    
    load( "../expectedResults/expectedResults.RData")
    expect_equal(expectedPathList, pathList)
    
  }
)


