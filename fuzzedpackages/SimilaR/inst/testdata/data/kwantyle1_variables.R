 kwantyle1<-function(x_5467784,p_5467784){
if (is.numeric(x_5467784) && is.numeric(p_5467784) && 0 <= p_5467784 && p_5467784 <=1)
{n_5467784=length(x_5467784);
h_5467784<-(n_5467784-1)*p_5467784+1;
return (sort(x_5467784)[floor(h_5467784)]+(h_5467784-floor(h_5467784))*(sort(x_5467784)[floor(h_5467784)+1]-sort(x_5467784)[floor(h_5467784)]));}
else
{NA}} 
