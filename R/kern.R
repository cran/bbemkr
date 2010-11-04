kern <-
function (data_x, x) 
{
    data_num = dim(data_x)[1]
    dim = dim(data_x)[2]
    hprod = prod(x)
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
            temp = 0
            for (k in 1:dim) {
                xa = (data_x[i, k] - data_x[j, k])/x[k]
                temp = temp + xa * xa
            }
            weight = cont * exp(-0.5 * temp)/hprod
            suma = suma + weight * data_y[j]
            sumb = sumb + weight
        }
        mh = suma/sumb
        temp = data_y[i] - mh
        sse = sse + temp * temp
    }
    return(1 - sse/sst)
}

