PriestleyChaokernel <-
function(x, y, h, gridpoint)
{
  n = length(y)
  mh = vector(,length(gridpoint))
  for(j in 1:length(gridpoint))
  {
      suma = vector(, n)
      suma[1] = x[1]*ker((gridpoint[j] - x[1])/h)*y[1]
      for(i in 2:n)
      {
          suma[i] = (x[i]-x[i-1])*ker((gridpoint[j] - x[i])/h)*y[i]
      }
      mh[j] = sum(suma)/h
  }
  return(list(gridpoint=gridpoint, mh=mh))
}

