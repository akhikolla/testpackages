 kendall2<-function(x,y) {
   n <- length(x)
   if (n<=0) {
      print('Wektor x jest pusty lub ma jeden element!')
      return(NA)
   } else if(n!=length(y)) {
      print('Wektory sa roznej dlugosci!')
      return(NA)
   } else {
      c<-0
      for (i in 1:(n-1)) {
         for (j in (i+1):n) {
            if(x[i]==x[j]) {
               print('wartosci x nie sa rozne!')
               return(NA)
               break
            }
            if(y[i]==y[j]) {
               print('wartosci y nie sa rozne!')
               return(NA)
               break
            }
            if(x[i]>x[j] & y[i]>y[j]) c <- c+1
            if(x[i]<x[j] & y[i]<y[j]) c <- c+1
         }
      }
     return(2*c/(n*(n-1)/2)-1) 
   }
}
