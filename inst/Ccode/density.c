double logdensity(xp,data_post)
double *xp,**data_post; 
{ 
	int i,j; 
	double hprod,cont,sum,temp,tem2,sigma,hatf,*band; 

	//08/09/2010: log multivariate kernel density 

	band=dvector(0,dim); 
	//normal reference rule 
	for(temp=0.0,j=0;j<=dim;j++) 
	{ 
		temp=0.0; 
		tem2=0.0; 
		for(i=1;i<=length;i++) 
		{ 
			temp += data_post[i][j]; 
			tem2 += data_post[i][j]*data_post[i][j]; 
		} 
		sigma=sqrt(tem2/length-(temp/length)*(temp/length)); 
		temp=exp(1.0/(1.0*(dim+1)+4.0)*log(4.0/(1.0*(dim+1)+2.0))); 
		band[j]=temp*sigma*exp(-1.0/(dim+1+4)*log(1.0*length)); 
	} 
	//compute log density
	for(hprod=1.0,j=0;j<=dim;j++) hprod *= band[j]; 
	cont=exp(-0.5*(dim+1)*log(2.0*pi)); 
	sum=0.0; 	
	for(i=1;i<=length;i++)
	{ 
		for(temp=0.0,j=0;j<=dim;j++) 
		{ 
				tem2=(xp[j]-data_post[i][j])/band[j]; 
				temp += tem2*tem2; 
		} 
		sum += cont*exp(-0.5*temp)/hprod; 
	}
	hatf=log(sum/length); 

	free_dvector(band,0,dim);

	return hatf; 
} 

double logpriors(xp) 
double *xp; 
{ 
	int i; 
	double sigma2,logf,*band; 

	//08/09/2010: log priors computed at posterior estimate
	band=dvector(1,dim); 
	sigma2=xp[0]*xp[0]; 

	for(i=1;i<=dim;i++) band[i]=xp[i]; 
	for(logf=0.0,i=1;i<=dim;i++) logf += -1.0*log(1.0+band[i]*band[i]); 
	logf += 0.5*prior_p*log(0.5*prior_st)-gammln(0.5*prior_p); 
	logf += -1.0*(0.5*prior_p)*log(sigma2)-0.5*prior_st/sigma2; 
	
	free_dvector(band,1,dim); 
	return logf; 
} 
//loglikelihood computed at the posterior estimate
double loglikelihood(he)
double *he;
{ 
	int i,j,k; 
	double hprod,cv,suma,sumb,cont,mh,temp,logf; 
	double weight,xa; 

	for(hprod=1.0,i=1;i<=dim;i++) hprod*=he[i]; 
	cont=exp(-0.5*dim*log(2.0*pi)); 
	cv=0.0; 
	suma=0.0; 
	sumb=0.0;     
	for(j=2;j<=data_num;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[1][k]-data_x[j][k])/he[k]; 
			temp+=xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma+=weight*data_y[j]; 
		sumb+=weight; 
	} 
	mh=suma/sumb; 
	cv+=(data_y[1]-mh)*(data_y[1]-mh); 
	
	for(i=2;i<=data_num-1;i++) 
	{ 
		suma=0.0; 
		sumb=0.0; 
		for(j=1;j<=i-1;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp+=xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma+=weight*data_y[j]; 
			sumb+=weight; 
		} 
		for(j=i+1;j<=data_num;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp+=xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma+=weight*data_y[j]; 
			sumb+=weight; 
		} 
		mh=suma/sumb; 
		cv += (data_y[i]-mh)*(data_y[i]-mh); 
	} 	
	suma=0.0; 
	sumb=0.0; 
	for(j=1;j<=data_num-1;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[data_num][k]-data_x[j][k])/he[k]; 
			temp += xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma+=weight*data_y[j]; 
		sumb+=weight; 
	} 
	mh=suma/sumb; 
	cv += (data_y[data_num]-mh)*(data_y[data_num]-mh); 
	//sigma2=he[0]*he[0] 
	logf = -0.5*data_num*log(2.0*pi*he[0]*he[0])-cv/(2.0*he[0]*he[0]); 

	return logf; 
} 

