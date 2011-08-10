loglikelihood <-
function(h2, data_x)
{
dm = dim(data_x)[2]
data_num = dim(data_x)[1]
hn = data_num^(-1/(4+dm)) 
h = vector(,dm)
sigma2 = h2[dm+1]
for(i in 1:dm)
{
h[i] = sqrt(h2[i])*hn
}
hprod = prod(h)
cont = exp(-0.5*dm*log(2.0*pi))
suma = sumb = 0
for(j in 2:data_num)
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[1,k] - data_x[j,k])/h[k])^2
}
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight * data_y[j]
sumb = sumb + weight
}
cv = (data_y[1] - suma/sumb)^2
for(i in 2:(data_num-1))
{
suma = sumb = 0
for(j in 1:(i-1))
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[i,k] - data_x[j,k])/h[k])^2
}
weight = cont * exp(-0.5*temp)/hprod
suma = suma + weight * data_y[j]
sumb = sumb + weight
}
for(j in (i+1):data_num)
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[i,k] - data_x[j,k])/h[k])^2
}
weight = cont * exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
sumb = sumb + weight
}
cv = cv + (data_y[i] - suma/sumb)^2
}
suma = sumb =0
for(j in 1:(data_num-1))
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[data_num, k] - data_x[j,k])/h[k])^2
}
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
sumb = sumb + weight
}
cv = cv + (data_y[data_num] - suma/sumb)^2
logp = -0.5*data_num*log(2.0*pi*sigma2) - cv/(2.0*sigma2)
return(logp)
}

