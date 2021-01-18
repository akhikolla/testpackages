
"checkcont"<-function(cont){
  if(!is.data.frame(cont)){
    cont <- as.data.frame(cont)
  }
  
  if("data.table" %in% class(cont)){
    cont <- as.data.frame(cont)
    setDF(cont)
  }
  
  if(!("age"    %in% colnames(cont))){stop("Column 'age' is missing in data frame cont.\n")}
  if(!is.numeric(cont$age)){          cont$age    <- as.numeric(cont$age)}
  if(any(is.na(cont$age))){           stop("Column 'age' of argument 'cont' contains NA values.\n")}
  if(!(all(cont$age==1:nrow(cont)))){ stop("Column 'age' of argument  'cont' must be identical to the row numbers.\n")}

  sexes <- "male" %in% colnames(cont)
  
  if(sexes){
    if(!("male"   %in% colnames(cont))){stop("Column 'male' is missing in data frame cont.\n")}
    if(!("female" %in% colnames(cont))){stop("Column 'female' is missing in data frame cont.\n")}
    if(!is.numeric(cont$female)){       cont$female <- as.numeric(cont$female)}
    if(!is.numeric(cont$male)){         cont$male   <- as.numeric(cont$male)}
    if(any(is.na(cont$male))){          stop("Column 'male' of argument 'cont' contains NA values.\n")}
    if(any(is.na(cont$female))){        stop("Column 'female' of argument 'cont' contains NA values.\n")}
    if(any(cont$male<0)){               stop("Column 'male' of argument 'cont' contains negative values.\n")}
    if(any(cont$female<0)){             stop("Column 'female' of argument 'cont' contains negative values.\n")}
  
    if(nrow(cont)>1 && abs(sum(cont$male[-1])-sum(cont$female[-1]))>0.00001){
      stop("The values sum(cont$male[-1]) and sum(cont$female[-1]) are not equal.\n")
    }
    if(abs(sum(cont$male)+sum(cont$female)-1)>0.00001){
      stop("The value sum(cont$male)+sum(cont$female) is not equal to 1.\n")
    }
    if(any(cont$male[-nrow(cont)]-cont$male[-1]<0)){
      stop("The contributions of age cohorts are not decreasing for males.\n")
    }
    if(any(cont$female[-nrow(cont)]-cont$female[-1]<0)){
      stop("The contributions of age cohorts are not decreasing for females.\n")
    }
    while(cont$female[nrow(cont)]+cont$male[nrow(cont)]==0){
      cont <- cont[1:(nrow(cont)-1),]
    }
    cont$cohort <- cont$male + cont$female
    
  }else{
    if(!("cohort"   %in% colnames(cont))){stop("Columns 'male' and 'female' or 'cohort' are missing in data frame cont.\n")}
    if(!is.numeric(cont$cohort)){         cont$cohort <- as.numeric(cont$cohort)}
    if(any(is.na(cont$cohort))){          stop("Column 'cohort' of argument 'cont' contains NA values.\n")}
    if(any(cont$cohort<0)){               stop("Column 'cohort' of argument 'cont' contains negative values.\n")}
    if(abs(sum(cont$cohort)-1)>0.00001){  stop("The value sum(cont$cohort) is not equal to 1.\n")}
    if(any(cont$cohort[-nrow(cont)]-cont$cohort[-1]<0)){stop("The contributions of age cohorts are not decreasing.\n")}
    
    while(cont$cohort[nrow(cont)]+cont$cohort[nrow(cont)]==0){
      cont <- cont[1:(nrow(cont)-1),]
    }
  }
  
  
  cont
}