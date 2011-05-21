bbeloglikelihood <-
function (data_x, data_y, x, sigma, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight",
          "Triangular", "Uniform")) 
{
    # Gaussian kernel has infinite support
    # Other types of kernel has [-1, 1] finite support

    kerntype = match.arg(kerntype)
    data_num = dim(data_x)[1]
    dim = dim(data_x)[2]
    suma = sumb = 0
    for (j in 2:data_num) {
        weight = prod(ker(as.numeric((data_x[1,]-data_x[j,])/x), kerntype = kerntype)/x)
        suma = suma + weight * data_y[j]
        sumb = sumb + weight
    }
    mh = suma/sumb
    cv = (data_y[1] - mh) * (data_y[1] - mh)
    for (i in 2:(data_num - 1)) {
        suma = sumb = 0
        for (j in 1:(i - 1)) {
            weight = prod(ker(as.numeric((data_x[j,]-data_x[i,])/x), kerntype = kerntype)/x)
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        for (j in (i + 1):data_num) {
            weight = prod(ker(as.numeric((data_x[j,]-data_x[i,])/x), kerntype = kerntype)/x)
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        mh = suma/sumb
        cv = cv + (data_y[i] - mh) * (data_y[i] - mh)
    }
    suma = sumb = 0
    for (j in 1:(data_num - 1)) {
        weight = prod(ker(as.numeric((data_x[data_num,]-data_x[j,])/x), kerntype = kerntype)/x)
        suma = suma + weight * data_y[j]
        sumb = sumb + weight
    }
    mh = suma/sumb
    cv = cv + (data_y[data_num] - mh) * (data_y[data_num] - mh)
    logf = -0.5 * data_num * log(2 * pi * sigma * sigma) - cv/(2 * 
        sigma * sigma)
    return(logf)
}

