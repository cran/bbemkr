bbewarmup <-
function (x, costpara, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight",
                      "Triangular", "Uniform"), warm = 4, sizep = 2) 
{
    # Gaussian kernel has infinite support
    # Other types of kernel has [-1, 1] finite support

    kerntype = match.arg(kerntype)
    dim = length(x)
    inicost = costpara
    xh = x
    xhtest = matrix(, warm, dim)
    sigmatest = accept_htest = costtest = vector(, warm)
    for (i in 1:warm) {
        dummy = np_gibbs(data_x, data_y, xh, inicost, kerntype = kerntype, sizep = sizep)
        xh = dummy$xh
        inicost = dummy$inicost
        xhtest[i, ] = dummy$xh
        costtest[i] = dummy$inicost
        sigmatest[i] = dummy$sigma
        accept_htest[i] = dummy$accept_h
    }
    return(list(xh = xhtest[warm, ], xhtest = xhtest, sigma = sigmatest[warm], 
        sigmatest = sigmatest, cost = costtest[warm], costtest = costtest, 
        accept_rate = sum(accept_htest)/warm))
}

