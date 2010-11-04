double cost(x) 
double *x; 
{ 
	int i,j,k; 
	double hprod,cv,suma,sumb,cont,mh,temp,logf; 
	double weight,xa,*he,lambda1,lambda2; 

	lambda1=0.25; 
	lambda2=0.25; 

	/*11/05/2008: likelihood cv: sample the bandwidth 
	for a nonparametric regression model*/ 

	he=dvector(1,dim); 

	hprod=1.0; 
	for(k=1;k<=dim;k++) 
	{ 
		he[k]=exp(x[k]); 
		hprod *= he[k]; 
	} 
	/*	
	if(he[1]<=0.001) return 1.0*exp(20.0); 
	if(he[2]<=0.001) return 1.0*exp(20.0); 
	*/ 
	cont=exp(-0.5*dim*log(2.0*pi)); /*Gaussian kernel function*/ 
	cv=0.0; 
	suma=0.0; 
	sumb=0.0; 
	/*sumb=cont/hprod;*/  
	for(j=2;j<=data_num;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[1][k]-data_x[j][k])/he[k]; 
			temp += xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma += weight*data_y[j]; 
		sumb += weight; 
	} 
	/*mh=(suma/(data_num-1))/(sumb/data_num);*/ 
	mh=suma/sumb; 
	cv += (data_y[1]-mh)*(data_y[1]-mh); 

	for(i=2;i<=data_num-1;i++) 
	{ 
		suma=0.0; 
		sumb=0.0; 
		/*sumb=cont/hprod;*/  
		for(j=1;j<=i-1;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp += xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma +=weight*data_y[j]; 
			sumb +=weight; 
		} 
		for(j=i+1;j<=data_num;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp += xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma +=weight*data_y[j]; 
			sumb +=weight; 
		} 
		/*mh=(suma/(data_num-1))/(sumb/data_num);*/ 
		mh=suma/sumb; 
		cv += (data_y[i]-mh)*(data_y[i]-mh); 
	} 
	suma=0.0; 
	sumb=0.0; 
	/*sumb=cont/hprod;*/ 
	for(j=1;j<=data_num-1;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[data_num][k]-data_x[j][k])/he[k]; 
			temp += xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma += weight*data_y[j]; 
		sumb += weight; 
	} 
	/*mh=(suma/(data_num-1))/(sumb/data_num);*/ 
	mh=suma/sumb; 
	cv += (data_y[data_num]-mh)*(data_y[data_num]-mh); 

	logf=-0.5*(1.0*data_num+prior_p)*log(0.5*cv+0.5*prior_st); 


	for(i=1;i<=dim;i++) logf += x[i]; 
	for(i=1;i<=dim;i++) logf += -1.0*log(1.0+he[i]*he[i]); 

	free_dvector(he,1,dim); 
	return -1.0*logf; 
} 

double cost2(x) 
double *x; 
{ 
	int i,j,k; 
	double hprod,cv,suma,sumb,cont,cont2,mh,temp; 
	double weight,xa,*he; 

	/*11/05/2008: cost fnction for sigma2*/ 

	he=dvector(1,dim); 

	hprod=1.0; 
	for(k=1;k<=dim;k++) 
	{ 
		he[k]=exp(x[k]); 
		hprod *= he[k]; 
	} 

	cont=exp(-0.5*dim*log(2.0*pi)); /*Gaussian kernel function*/ 
	cv=0.0; 
	suma=0.0; 
	sumb=0.0; 
	/*sumb=cont/hprod;*/ 
	for(j=2;j<=data_num;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[1][k]-data_x[j][k])/he[k]; 
			temp += xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma += weight*data_y[j]; 
		sumb += weight; 
	} 
	/*mh=(suma/(data_num-1))/(sumb/data_num);*/
	mh=suma/sumb; 
	cv += (data_y[1]-mh)*(data_y[1]-mh); 

	for(i=2;i<=data_num-1;i++) 
	{ 
		suma=0.0; 
		/*sumb=cont/hprod;*/ 
		sumb=0.0; 
		for(j=1;j<=i-1;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp += xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma +=weight*data_y[j]; 
			sumb +=weight; 
		} 
		for(j=i+1;j<=data_num;j++) 
		{ 
			for(temp=0.0,k=1;k<=dim;k++) 
			{ 
				xa=(data_x[i][k]-data_x[j][k])/he[k]; 
				temp += xa*xa; 
			} 
			weight=cont*exp(-0.5*temp)/hprod; 
			suma +=weight*data_y[j]; 
			sumb +=weight; 
		} 
		/*mh=(suma/(data_num-1))/(sumb/data_num);*/ 
		mh=suma/sumb; 
		cv += (data_y[i]-mh)*(data_y[i]-mh); 
	} 
	suma=0.0; 
	sumb=0.0; 
	/*sumb=cont/hprod;*/ 
	for(j=1;j<=data_num-1;j++) 
	{ 
		for(temp=0.0,k=1;k<=dim;k++) 
		{ 
			xa=(data_x[data_num][k]-data_x[j][k])/he[k]; 
			temp += xa*xa; 
		} 
		weight=cont*exp(-0.5*temp)/hprod; 
		suma += weight*data_y[j]; 
		sumb += weight; 
	} 
	/*mh=(suma/(data_num-1))/(sumb/data_num);*/ 
	mh=suma/sumb; 
	cv += (data_y[data_num]-mh)*(data_y[data_num]-mh); 
	cv += prior_st; 

	free_dvector(he,1,dim); 

	return 0.5*cv; 
} 

