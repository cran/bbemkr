GasserMullerkernel <-
function(gridpoint, x, y, h)
{
result=vector(,length(gridpoint))
integrand = function(u)
{    
ker((gridpoint-u)/h)
}
inteval = res = vector(,length(x)-2)
for(i in 3:length(x))
{
j=i-2
s_i=(x[i-1]+x[i])/2
s_iminus1=(x[i-2]+x[i-1])/2
res[j] = integrate(integrand, lower = s_iminus1, upper = s_i)$value*y[i]
}
if(sum(res)/h<=0)
{
da = 0
warning("Density is negative.")
}
else
{
da = sum(res)/h
}
return(da)
}

