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

 kendall2_2 <-function(x_456842,y_456842) {
   n_456842 <- length(x_456842)
   if (n_456842<=0) {
      print('Wektor x jest pusty lub ma jeden element!')
      return(NA)
   } else if(n_456842!=length(y_456842)) {
      print('Wektory sa roznej dlugosci!')
      return(NA)
   } else {
      c_456842<-0
      for (i_456842 in 1:(n_456842-1)) {
         for (j_456842 in (i_456842+1):n_456842) {
            if(x_456842[i_456842]==x_456842[j_456842]) {
               print('wartosci x nie sa rozne!')
               return(NA)
               break
            }
            if(y_456842[i_456842]==y_456842[j_456842]) {
               print('wartosci y nie sa rozne!')
               return(NA)
               break
            }
            if(x_456842[i_456842]>x_456842[j_456842] & y_456842[i_456842]>y_456842[j_456842]) c_456842 <- c_456842+1
            if(x_456842[i_456842]<x_456842[j_456842] & y_456842[i_456842]<y_456842[j_456842]) c_456842 <- c_456842+1
         }
      }
     return(2*c_456842/(n_456842*(n_456842-1)/2)-1) 
   }
} 

