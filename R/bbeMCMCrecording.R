bbeMCMCrecording <-
function (data_x, data_y, x, costpara, kerntype = c("Gaussian", "Epanechnikov", 
                             "Quartic", "Triweight", "Triangular", "Uniform"), num_batch = 50, 
                             M = 10000, step = 20, sizep = 1.2) 
{
    # Gaussian kernel has infinite support
    # Other types of kernel has [-1, 1] finite support

    kerntype = match.arg(kerntype)
    dim = length(x)
    size_batch = M/num_batch
    length2 = M/step
    sum_h = total_sd = xp = h = sif = std_h = vector(, dim + 
        3)
    batch_h = matrix(, num_batch, dim + 3)
    xhtest = matrix(, M, dim)
    sigmatest = accept_htest = costtest = plike = vector(, M)
    for (i in 1:(dim + 3)) {
        sum_h[i] = 0
        for (j in 1:num_batch) {
            batch_h[j, i] = 0
        }
        total_sd[i] = 0
    }
    data_post = matrix(, length2, dim + 1)
    inicost = costpara
    xh = x
    for (k in 1:M) {
        dummy = np_gibbs(data_x, data_y, xh, inicost, kerntype = kerntype, sizep = sizep)
        xh = dummy$xh
        sigma = dummy$sigma
        inicost = dummy$inicost
        xhtest[k, ] = dummy$xh
        costtest[k] = dummy$inicost
        sigmatest[k] = dummy$sigma
        accept_htest[k] = dummy$accept_h
        sn = ceiling(k/size_batch)
        xp[1:dim] = dummy$xh
        xp[dim + 1] = dummy$sigma
        xp[dim + 2] = plike[k] = bbelike(data_x, data_y, xh, 
            sigma, kerntype = kerntype)
        xp[dim + 3] = dummy$inicost
        for (i in (dim + 2):(dim + 3)) {
            temp = xp[i]
            sum_h[i] = sum_h[i] + temp
            batch_h[sn, i] = batch_h[sn, i] + temp
            total_sd[i] = total_sd[i] + temp * temp
        }
        temp = sqrt(xp[dim + 1])
        sum_h[dim + 1] = sum_h[dim + 1] + temp
        batch_h[sn, dim + 1] = batch_h[sn, dim + 1] + temp
        total_sd[dim + 1] = total_sd[dim + 1] + temp * temp
        for (j in 1:dim) {
            temp = exp(xp[j])
            sum_h[j] = sum_h[j] + temp
            batch_h[sn, j] = batch_h[sn, j] + temp
            total_sd[j] = total_sd[j] + temp * temp
        }
        index = ceiling(k/step)
        data_post[index, dim + 1] = sqrt(xp[dim + 1])
        for (j in 1:dim) {
            data_post[index, j] = exp(xp[j])
        }
    }
    for (i in 1:(dim + 3)) {
        std_h[i] = 0
        sum_h[i] = sum_h[i]/M
        for (j in 1:num_batch) {
            temp = batch_h[j, i]/size_batch - sum_h[i]
            std_h[i] = std_h[i] + temp * temp
        }
        std_h[i] = sqrt(std_h[i]/(num_batch * num_batch - num_batch))
        total_sd[i] = sqrt(abs(M * (total_sd[i]/M - sum_h[i] * 
            sum_h[i])/(M - 1)))
        h[i] = sum_h[i]
    }
    for (i in 1:(dim + 3)) {
        sif[i] = std_h[i] * std_h[i]/(total_sd[i] * total_sd[i]) * 
            M
    }
    temp = 0
    for (k in 1:M) {
        plike[k] = plike[k] - sum_h[dim + 2]
        temp = temp + 1/exp(plike[k])
    }
    temp = temp/M
    logmarginal = log(1/temp) + sum_h[dim + 2]
    logmargin = bbeloglikelihood(data_x, data_y, sum_h[1:dim], 
        sum_h[dim + 1], kerntype = kerntype)
    temp = bbelogpriors(sum_h[1:dim], sum_h[dim + 1])
    temp2 = bbelogdensity(sum_h[1:dim], sum_h[dim + 1], data_post)
    logmarginal2 = logmargin + temp - temp2
    return(list(accept_raterecording = sum(accept_htest)/M, sum_h = sum_h, 
        std_h = std_h, batch_h = batch_h, total_sd = total_sd, 
        sif = sif, R2 = kern(data_x, sum_h[1:dim], kerntype = kerntype), data_post = data_post, 
        logmarginalNR = logmarginal, loglikelihoods = logmargin, 
        logprior = temp, logdensity = temp2, logmarginalChib = logmarginal2))
}

