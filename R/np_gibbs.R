np_gibbs <-
function(data_x, data_y, xh, inicost, prior_p = 2, sizep)
{
    data_num = dim(data_x)[1]
    dim = dim(data_x)[2]    
    dv = rn = vector(,dim)
    fx = inicost
    sum = 0
    for(i in 1:dim)
    {
        rn[i] = rnorm(1)
     sum = sum + rn[i] * rn[i]
    }
for(i in 1:dim)
{
    dv[i] = rn[i]/sqrt(sum) * rnorm(1) * sizep
xh[i] = xh[i] + dv[i]
}
fy = bbecost(data_x, data_y, xh)
r = -1 * (fy - fx)
if(r > 0)
{
   accept = 1
}
else
{
   un = 0
   while(un <= 0)
   {
      un = runif(1)
   }
   if(un < exp(r))
   {
      accept = 1
   }
else
{ 
   accept = 0
}
}
accept_h = 0
if(accept == 1)
{
   accept_h = accept_h + 1
   inicost = fy
   xh = xh
}
else
{
    xh = xh - dv
    inicost = fx
}
temp = bbecost2(data_x, data_y, xh)
un = rgamma(1, shape = 0.5 * (data_num + prior_p), rate = temp)
sigma = 1/un
return(list(xh = xh, sigma = sigma, inicost = inicost, accept_h = accept_h))
}

