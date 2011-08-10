cost <-
function(x, data_x, data_y, prior_p = 2, prior_st = 1)
{
data_num = dim(data_x)[1]
dm = dim(data_x)[2]

hn = data_num^(-1/(4+dm)) 
tau2 = h = vector(,dm)
for(k in 1:dm)
{
tau2[k] = exp(x[k])
h[k] = sqrt(tau2[k])*hn
}
hprod = prod(h)
cont = exp(-0.5 * dm * log(2.0 * pi))
suma = sumb = 0 
for(j in 2:data_num)
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[1,k] - data_x[j,k])/h[k])^2
}
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
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
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
sumb = sumb + weight
}
for(j in (i+1):data_num)
{
temp = 0 
for(k in 1:dm)
{
temp = temp + ((data_x[i,k] - data_x[j,k])/h[k])^2
}
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
sumb = sumb + weight
}
cv = cv + (data_y[i] - suma/sumb)^2
}
suma = sumb = 0 
for(j in 1:(data_num-1))
{
temp = 0
for(k in 1:dm)
{
temp = temp + ((data_x[data_num,k] - data_x[j,k])/h[k])^2
}
weight = cont*exp(-0.5*temp)/hprod
suma = suma + weight*data_y[j]
sumb = sumb + weight
}
cv = cv + (data_y[data_num] - suma/sumb)^2
logp = -0.5*(data_num+prior_p)*log(0.5*cv+0.5*prior_st)
for(i in 1:dm)
{
logp = logp + x[i] + logpriorh2(tau2[i]*hn*hn)+log(hn*hn)
}
return(-logp)
}

