/*generate random numbers form Gamma(a,b)*/
float Rgamma(a,b) 
float a, b; 
{ 
	int ok; 
	double d,q,un,u1,y,z; 

	if(a<=0.0 || b<=0.0) 
	{ 
		printf("Gamma parameter error (<0.0)\n"); 
		return 0.0;  
	} 

	if(a<1.0) 
	{ 
		/* Ahrens, P.213 */ 
		ok=0; 
		while(ok==0) 
		{ 
			un=0.0; 
			while(un<=0.0 || un>=1.0) un=rand()*1.0/RAND_MAX; 
			d=(2.718282+a)/2.718282; 
			q=d*un; 

			if(q<=1.0) 
			{ 
				z=exp(1.0/a*log(q)); 
				u1=rand()*1.0/RAND_MAX; 
				if(u1<exp(-z)) ok=1; 
			} 
			else 
			{ 
				z=-1.0*log((d-q)/a); 
				u1=rand()*1.0/RAND_MAX; 
				if(u1<exp((a-1)*log(z))) ok=1; 
			} 
		} /* end ok */ 
	} 
	else 
	{ 
		/* a>=1.0 Fishman, P.214 */ 
		ok=0; 
		while(ok==0) 
		{ 
			un=0.0; 
			while(un<=0.0 || un>=1.0) un=rand()*1.0/RAND_MAX; 
			y=-1.0*log(un); 

			u1=rand()*1.0/RAND_MAX; 
			if(u1<exp((a-1)*(log(y)-(y-1)))) { z=a*y; ok=1; } 
		} 
	} 
	return z/b; 
} 

double dloggauss(z,mu,sd) 
double z,mu,sd; 
{ 
	double sum; 
	
	sum=-0.5*log(2.0*pi*sd*sd); 
	sum+=-0.5*(z-mu)*(z-mu)/sd/sd; 

	return sum; 
} 


double gasdev() 
{ 
	static int iset=0; 
	static double gset; 
	double fac,r,v1,v2; 

	if(iset==0) 
	{ 
		do 
		{ 
			v1=rand()*2.0/RAND_MAX-1.0; 
			v2=rand()*2.0/RAND_MAX-1.0; 
			r=v1*v1+v2*v2; 
		} 
		while(r>=1.0); 
		fac=sqrt(-2.0*log(r)/r); 
		gset=v1*fac; 
		iset=1; 
		return v2*fac; 
	} 
	else 
	{ 
		iset=0; 
		return gset; 
	} 
} 

float gammln(xx)
double xx;
{
    double x,y,tmp,ser;
    double cof[6];
        int j;
 
        cof[0]=76.18009172947146;
        cof[1]=-86.5053203294167;
        cof[2]=24.01409824083091;
        cof[3]=-1.231739572450155;
        cof[4]=0.1208650973866179e-2;
        cof[5]=-0.5395239384953e-5;

    y=xx;
    x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
}
