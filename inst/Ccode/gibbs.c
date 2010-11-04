int np_gibbs(xh) 
double *xh; 
{
	int i,accept; 
	double fx,sum,fy,r,un,temp; 
	double *dv,*rn; 

	/*14/03/2008: block-move for sampling bandwidths*/

	dv=dvector(1,dim);
	rn=dvector(1,dim);

	fx=xh[-1]; 
	for(sum=0.0,i=1;i<=dim;i++) 
	{ 
		rn[i]=gasdev(); 
		sum+=rn[i]*rn[i]; 
	} 
	for(i=1;i<=dim;i++) 
	{ 
		dv[i]=rn[i]/sqrt(sum)*gasdev()*mutsizp;  
		xh[i] += dv[i]; 
	} 
	fy=cost(xh); 

	r=-1.0*(fy-fx); 
	if(r>0.0) accept=1; 
	else 
	{ 
		un=0.0; 
		while(un<=0.0) un=rand()*1.0/RAND_MAX; 
		if(un<exp(r)) accept=1; 
		else accept=0; 
	} 
	if(accept==1) 
	{ 
		accept_h++; 
		xh[-1]=fy; 
	} 
	else 
	{ 
		for(i=1;i<=dim;i++) xh[i] -= dv[i]; 
	} 
	total_h++; 

	temp=cost2(xh); 
	un=0.0; 
	while(un<=0.0) un=Rgamma(0.5*(1.0*data_num+prior_p),temp); 
	xh[0]=1.0/un; 

	free_dvector(dv,1,dim); 
	free_dvector(rn,1,dim); 

	return 0; 
} 

