#' Gives summaries of GetSubTite Objects.
#' @param Z List produced by GetSubTite.
#' @export
Print_SubTite = function(Z){

  ###
  DESIGN = Z[[1]]
  Borrow = DESIGN[[length(DESIGN)]]


  cat(DESIGN[[1]])
  cat("
")
  cat(DESIGN[[2]])
  cat("

")
  cat(names(DESIGN)[[3]])
  cat("
")
  cat(DESIGN[[3]])
  cat("

")

  cat(DESIGN[[6]])
  cat("

")

  Target = DESIGN[[4]]

  cat(paste0(names(DESIGN)[[length(DESIGN)-1]], DESIGN[[length(DESIGN)-1]]))
  cat("
")

  cat(names(DESIGN)[[4]])
  cat(DESIGN[[4]])
  cat("

")
  cat(names(DESIGN)[[5]])
  cat(DESIGN[[5]])
  cat("

")
  cat("Hyperparameters:
")
  for(k in 7:8){
    cat(paste0(names(DESIGN)[[k]],DESIGN[[k]],"
"))
  }


  cat(paste0(names(DESIGN)[[9]]))
  cat(paste0(DESIGN[[9]]))
  cat("
")


  cat(paste0(names(DESIGN)[[10]]))
  cat(paste0(DESIGN[[10]]))
  cat("
")



  for(k in 11:12){
    cat(paste0(names(DESIGN)[[k]],DESIGN[[k]],"
"))
  }

  cat("
")
  cat(paste0(names(DESIGN)[[13]],"
",
DESIGN[[13]],"

"))

  cat(paste0(names(DESIGN)[[length(DESIGN)]], DESIGN[[length(DESIGN)]]))

  cat("

Trial Data
")
  print(Z$Data)


  cat("
Summary of data included in model

")

  cat("Number of Treated Patients and DLTs
")

  print(Z$`Number of Treated and DLTs`)

  cat("
Number of completely evaluated patients at each subgroup dose
")
  print(Z$`Number Fully Evaluated At Each Dose`)


  if(Borrow==2){
    cat("Clustering Estimates
")
    print(Z$`Clustering Parameters`)


  }


  ###Lastly do dose assignment

  cat("
Optimal Doses
")
  for(k in 1:length(Target)){
    cat(paste0("Next recommended dose level for subgroup ",k, ": Dose ", Z$`Optimal Dose`[k],"
"))
  }


}
