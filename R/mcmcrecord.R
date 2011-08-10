mcmcrecord <-
function(x, inicost, mutsizp,
  warm = 100, M = 100, prob = 0.234, num_batch = 10, step = 10, data_x, 
  data_y, xm, prior_p = 2, prior_st = 1)
{
dm = dim(data_x)[2]
data_num = dim(data_x)[1]
hn = data_num^(-1/(4+dm)) 
size_batch = M/num_batch
sum_h = total_sd = sif = vector(,(dm+2))
ave_h2 = h2 = vector(,(dm+1))
batch_h = matrix(,num_batch,(dm+2))
for(i in 1:(dm+2))
{
sum_h[i] = 0
for(j in 1:num_batch)
{
batch_h[j,i] = 0
}
total_sd[i] = 0
}
for(i in 1:(dm+1))
{
ave_h2[i] = 0
}
mutsizprecord = sigma2record = costrecord = acceptance = vector(,M)
xrecord = matrix(,M,dm) 
cpost = matrix(,M/step,dm+1)
for(k in 1:M)
{
dummy = np_gibbs(x, inicost, k+warm, mutsizp, prob=prob, data_x = data_x, data_y = data_y)
xrecord[k,] = x = dummy$x
sigma2record[k] = dummy$sigma2
costrecord[k] = inicost = dummy$cost
mutsizprecord[k] = mutsizp = dummy$mutsizp
acceptance[k] = dummy$accept_h
sn = ceiling(k/size_batch)
index = ceiling(k/step)
sum_h[dm+2] = sum_h[dm+2] + costrecord[k]
batch_h[sn,dm+2] = batch_h[sn,dm+2] + costrecord[k]
total_sd[dm+2] = total_sd[dm+2] + (costrecord[k])^2
sum_h[dm+1] = sum_h[dm+1] + sqrt(sigma2record[k])
ave_h2[dm+1] = ave_h2[dm+1] + sigma2record[k]
batch_h[sn,dm+1] = batch_h[sn,dm+1] + sqrt(sigma2record[k])
total_sd[dm+1] = total_sd[dm+1] + sigma2record[k]
cpost[index,dm+1] = sigma2record[k]
for(j in 1:dm)
{
temp = exp(xrecord[k,j])
sum_h[j] = sum_h[j] + sqrt(temp)
ave_h2[j] = ave_h2[j] + temp
batch_h[sn,j] = batch_h[sn,j] + sqrt(temp)
total_sd[j] = total_sd[j] + temp
cpost[index,j] = temp
}
}
std_h = vector(,(dm+2))
for(i in 1:(dm+2))
{
std_h[i] = 0
sum_h[i] = sum_h[i]/M
for(j in 1:num_batch)
{
std_h[i] = std_h[i] + (batch_h[j,i]/size_batch-sum_h[i])^2
}
std_h[i] = sqrt(std_h[i]/(num_batch*num_batch-num_batch))
total_sd[i] = sqrt(M*(total_sd[i]/M-sum_h[i]*sum_h[i])/(M-1))
sif[i] = (std_h[i]/total_sd[i])^2*M
}
h2 = ave_h2/M
logmargin = loglikelihood(h2, data_x)
logmargin = logmargin + logpriors(h2, data_x, prior_p = prior_p, prior_st = prior_st) - 
logdensity(h2, cpost)
if(missing(xm))
{
return(list(sum_h = sum_h, h2 = h2, sif = sif, mutsizp = mutsizprecord, 
cpost = cpost, accept = mean(acceptance), Chib = logmargin))
}
else
{
ker = kern(sum_h, hn = hn, data_x = data_x, xm = xm)
R2 = ker$r2
MSE = ker$mse
return(list(sum_h = sum_h, h2 = h2, sif = sif, mutsizp = mutsizprecord, 
cpost = cpost, accept = mean(acceptance), Chib = logmargin, R2 = R2, MSE = MSE))
}
}

