kern <-
function (data_x, x, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight",
                 "Triangular", "Uniform")) 
{
    # Gaussian kernel has infinite support
    # Other types of kernel has [-1, 1] finite support

    kerntype = match.arg(kerntype)
    data_num = dim(data_x)[1]
    dim = dim(data_x)[2]
    cont = exp(-0.5 * dim * log(2 * pi))
    suma = sumb = 0
    for (i in 1:data_num) {
        suma = suma + data_y[i]
        sumb = sumb + data_y[i] * data_y[i]
    }
    suma = suma/data_num
    sst = sumb - suma * suma * data_num
    sse = 0
    for (i in 1:data_num) {
        suma = sumb = 0
        for (j in 1:data_num) {
            weight = prod(ker(as.numeric((data_x[i,]-data_x[j,])/x), kerntype = kerntype)/x)
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        mh = suma/sumb
        temp = data_y[i] - mh
        sse = sse + temp * temp
    }
    return(1 - sse/sst)
}

