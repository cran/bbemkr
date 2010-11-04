int kern(fit,h)
double *fit,*h; 
{ 
	int i,j,k; 
	double hprod,suma,sumb,cont,mh,temp; 
	double weight,xa,r2,sst,sse,ase; 

	/*11/05/2008
	compute R^2 for a multivariate kernel regression model*/

	hprod=1.0;
	for(k=1;k<=dim;k++) hprod *= h[k];

	cont=exp(-0.5*dim*log(2.0*pi)); /*Gaussian kernel function*/
        suma=0.0;
        sumb=0.0;
        for(i=1;i<=data_num;i++)
        {
                suma += data_y[i];
                sumb += data_y[i]*data_y[i];
        }
        suma=suma/data_num;
        sst=sumb-suma*suma*data_num; /*end SST*/

	sse=0.0; 
	ase=0.0; 
	for(i=1;i<=data_num;i++) 
	{ 
		suma=0.0; 
		sumb=0.0; 
		for(j=1;j<=data_num;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/h[k]; 
				temp += xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma +=weight*data_y[j]; 
			sumb +=weight; 
		} 
		mh=suma/sumb; 
		temp=data_y[i]-mh; 
		sse += temp*temp; 
		temp=xm[i]-mh;
		ase += temp*temp;
	} 

	r2=1.0-sse/sst; 
	temp=ase/data_num;
	fit[1]=r2;
	fit[2]=temp;

	return 0; 
} 
