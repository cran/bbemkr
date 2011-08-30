LaplaceMetropolis <-
function(theta, data=NULL, prior_p, prior_st, method=c("likelihood","L1center","median"))
{
method = match.arg(method)
theta = as.matrix(theta)
niter = length(theta[,1])
p = length(theta[1,])
if(method == "likelihood")
{
h = NULL
for(t in (1:niter))
{
h = c(h, margin_like(theta[t,], data) + margin_prior(theta[t,], data, prior_p = prior_p, prior_st = prior_st))
hmax = max(h)
}
}
if(method == "L1center")
{
L1sum = NULL
oneniter = as.matrix(rep(1,niter))
onep = as.matrix(rep(1,p))
for(t in (1:niter))
{
thetat = theta[t,]
thetatmat = oneniter %*% thetat
L1sum = c(L1sum, sum(abs((theta-oneniter%*%thetat)%*%onep)))
}
argL1center = min((1:niter)[L1sum == min(L1sum)])
thetaL1center = theta[argL1center,]
hmax = margin_like(thetaL1center, data) + margin_prior(thetaL1center, data, prior_p = prior_p, prior_st = prior_st)
}
if(method == "median")
{
thetamed = apply(theta,2,median)
hmax = margin_like(thetamed, data) + margin_prior(thetamed, data, prior_p = prior_p, prior_st = prior_st)
}
if(p==1)
{
logdetV = 2*log(mad(theta[,1]))
}
else
{
logdetV = sum(log(eigen(cov.mve(theta)$cov)$values))
}
return(hmax + 0.5*p*log(2*pi)+0.5*logdetV)
}

