bbecost <-
function (data_x, data_y, x, kerntype = c("Gaussian", "Epanechnikov", "Quartic", 
    "Triweight", "Triangular", "Uniform"), bandx_priors, prior_p = 2, prior_st = 0.1) 
{
    # Gaussian kernel has infinite support
    # Other types of kernel has [-1, 1] finite support
    
    kerntype = match.arg(kerntype)
    data_num = dim(data_x)[1]
    dim = dim(data_x)[2]
    he = vector(, dim)
    for (k in 1:dim) {
        he[k] = exp(x[k])
    }
    hprod = prod(he)
    suma = sumb = 0
    for (j in 2:data_num) {
         weight = prod(ker(as.numeric((data_x[1,]-data_x[j,])/he), kerntype = kerntype)/he)        
         suma = suma + weight * data_y[j]
         sumb = sumb + weight
    }
    mh = suma/sumb
    cv = (data_y[1] - mh) * (data_y[1] - mh)
    for (i in 2:(data_num - 1)) {
        suma = sumb = 0
        for (j in 1:(i - 1)) {
            weight = prod(ker(as.numeric((data_x[j,]-data_x[i,])/he), kerntype = kerntype)/he)
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        for (j in (i + 1):data_num) {
            weight = prod(ker(as.numeric((data_x[j,]-data_x[i,])/he), kerntype = kerntype)/he)
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        mh = suma/sumb
        cv = cv + (data_y[i] - mh) * (data_y[i] - mh)
    }
    suma = sumb = 0
    for (j in 1:(data_num - 1)) {
        weight = prod(ker(as.numeric((data_x[data_num,]-data_x[j,])/he), kerntype = kerntype)/he)
        suma = suma + weight * data_y[j]
        sumb = sumb + weight
    }
    mh = suma/sumb
    cv = cv + (data_y[data_num] - mh) * (data_y[data_num] - mh)
    logf = -0.5 * (data_num + prior_p) * log(0.5 * cv + 0.5 * 
        prior_st) + sum(x)
    if (missing(bandx_priors)) {
        for (i in 1:dim) {
            logf = logf - 1 * log(1 + he[i] * he[i])
        }
    }
    else {
        logf = logf + bandx_priors
    }
    return(-logf)
}

