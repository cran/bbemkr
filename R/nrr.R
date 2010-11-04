nrr <-
function (data_x, logband = TRUE) 
{
    dim = dim(data_x)[2]
    data_num = dim(data_x)[1]
    h = vector(, dim)
    for (j in 1:dim) {
        temp = sum(data_x[, j])
        temp2 = sum(data_x[, j] * data_x[, j])
        sigma = sqrt(temp2/data_num - (temp/data_num) * (temp/data_num))
        temp = exp(1/(dim + 4) * log(4/(dim + 2)))
        h[j] = temp * sigma * exp(-1/(dim + 4) * log(data_num))
    }
    ifelse(logband == TRUE, return(log(h)), return(h))
}

