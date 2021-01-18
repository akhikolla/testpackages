  bucketsort1<-function(x){
   
   # nie wiedzialam jak sobie poradzic, wiec przyjelam, ze bede zaokraglac
   # liczby do czterech miejsc po przecinku
   x_1<-round(x,digits=-4) 
   y<-x*(10^4) # robie tak, zeby dostac liczby calkowite
   k<-max(y) # jako k wybiram najwieksza z liczb
k <- min(x)
   s<-f(k,y) # stosuje funkcje z poprzedniego zadania
   v<-NULL # poczatkowy wektor wynikowy jest pusty
   for(i in 1:k)
      if (s[i]!=0){
         v<-c(v,rep(i,s[i]))
      }
	else
{
v <- 5
}
   v<-v/(10^4) # dziele przez 10^4 bo na poczatku przez tyle pomnozylam
v<-v*(10^4)
   v
}
