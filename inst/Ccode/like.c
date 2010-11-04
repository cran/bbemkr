double like(xp)
double *xp;
{ 
	int i,j,k; 
	double hprod,cv,suma,sumb,cont,mh,temp,logf; 
	double weight,xa,*he; 

	he=dvector(1,dim); 
	hprod=1.0; 
	for(i=1;i<=dim;i++) 
	{ 
		he[i]=exp(xp[i]); 
		hprod*=he[i]; 
	} 
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
	//sigma2=xp[0] 
	logf = -0.5*data_num*log(2.0*pi*xp[0])-cv/(2.0*xp[0]); 

	free_dvector(he,1,dim); 

	return logf; 
} 

