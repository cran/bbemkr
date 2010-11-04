bbelogdensity <-
function (x, sigma, data_post) 
{
   dim = length(x)
   band = vector(,dim+1)
   len = dim(data_post)[1]
   for(j in 1:(dim+1))
   {
       temp = tem2 = 0
   for(i in 1:len)
   {
        temp = temp + data_post[i,j]
tem2 = tem2 + data_post[i,j] * data_post[i,j]
   }
   sig = sqrt(abs(tem2/len-(temp/len)*(temp/len)))
   temp = exp(1/(dim+5)*log(4/(dim+3)))
   band[j] = temp * sig * exp(-1/(dim+5)*log(len))
   }
       hprod = prod(band)
   cont=exp(-.5*(dim+1)*log(2*pi))
   sum = 0
   xp = c(x, sigma)
   for(i in 1:len)
   {
       temp = 0   
   for(j in 1:(dim+1))
   {
       tem2 = (xp[j] - data_post[i,j])/band[j]
   temp = temp + tem2 * tem2
}
sum = sum + cont * exp(-0.5*temp)/hprod
}
hatf = log(sum/len)
return(hatf)
}

