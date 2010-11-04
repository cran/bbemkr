bbelogpriors <-
function(x, sigma, prior_p = 2, prior_st = 0.1)
{
   dim = length(x)
   band = vector(,dim)
   sigma2 = sigma * sigma
   for(i in 1:dim)
   {
       band[i] = x[i]
    }
logf = 0
for(i in 1:dim)
{
    logf = logf - log(1 + band[i] * band[i])
}
logf = logf + 0.5 * prior_p * log(0.5 * prior_st) - lgamma(0.5 * prior_p)
logf = logf - 1.0*(0.5 * prior_p) * log(sigma2) - 0.5 * prior_st/sigma2
return(logf)
}

