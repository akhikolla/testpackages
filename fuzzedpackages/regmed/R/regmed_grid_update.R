regmed_grid_update <-
function(regmed.fit.obj,regmed.inits){

   regmed.inits$Alpha <- regmed.fit.obj$alpha
   regmed.inits$Beta <- regmed.fit.obj$beta
   regmed.inits$Delta <- regmed.fit.obj$delta
   regmed.inits$varx <- regmed.fit.obj$var.x
   regmed.inits$vary <- regmed.fit.obj$var.y

   return(regmed.inits)

}
