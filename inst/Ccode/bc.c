#include "comm.h" 
#include "nrutil.h" 
#include "nrutil.c" 
#include "lib.c" 
#include "cost.c" 
#include "gibbs.c" 
#include "like.c" 
#include "kern.c" 
#include "density.c" 

int main()
{
	long ltime; 
	int warm=4,M=10,step=2; 
       int ik_start, ik_end;
	int stime,i,j,k,ik,num_batch,size_batch,sn,index,length; 
	double *x,*sum_h,**batch_h,fvalue,*std_h,*h,*total_sd,*sif; 
	double *ave,**cov,**res,**data_yx,*plike,**data_post,**simu_res,*fit; 
	double un,temp,temp2,sigma,suma,sumb,sig2,r2,sse,logmargin; 
	FILE *ins; 
	char ap_filename[100];
	char ar_filename[100];

	/*01/09/2010: Monte Carlo simulation
	Nonparametric regression with unknown error density
	x[-2]: log likelihood
	x[-1]: negative log posterior
	x[0] : sigma
	x[1:dim]: log bandwidths 
       x[dim+1]: Newton-Raftery (1994) log marginal likelihood
       x[dim+2]: Chib (1995) log marginal likelihood */      

	num_batch=5; 
	size_batch=M/num_batch; 
	length=M/step; 
	
	ave=dvector(0,dim); 
	cov=dmatrix(0,dim,0,dim); 
	x=dvector(-2,dim); 
	h=dvector(-2,dim); 
	batch_h=dmatrix(1,num_batch,-2,dim); 

	total_sd=dvector(-2,dim); 
	sif=dvector(-2,dim); 
	sum_h=dvector(-2,dim); 
	std_h=dvector(-2,dim); 
	data_yx=dmatrix(1,data_num,1,dim+1); 
	plike=dvector(1,M); 
	data_post=dmatrix(1,length,0,dim); 
	simu_res=dmatrix(1,mc,-2,dim+5); 
	fit=dvector(1,2); 

	ltime=time(NULL); 
	stime=(unsigned int)ltime/2; 
	//stime=585933673;
	srand(stime); 
		
	sprintf(ap_filename, "ap_%04d-%04d.txt", ik_start, ik_end);
	sprintf(ar_filename, "ar_%04d-%04d.txt", ik_start, ik_end);
	
	ins=fopen("/srv/home/hanshang/data_kernel/data2Dfirst50.txt","r"); 
	if(ins==NULL) { printf("cannot open data file\n"); return; }  
	for(i=1;i<=data_num;i++) 
	{ 
		for(j=1;j<=dim+1;j++) fscanf(ins," %lf",&data_yx[i][j]); 
	} 
	fclose(ins); 
       
		for(i=1;i<=data_num;i++) 
		{ 
			data_y[i]=data_yx[i][1]; 
			for(j=1;j<=dim;j++) data_x[i][j]=data_yx[i][j+1]; 
		} 
		//normal reference rule
		for(temp=0.0,j=1;j<=dim;j++) 
		{ 
			temp=0.0; 
			temp2=0.0; 
			for(i=1;i<=data_num;i++) 
			{ 
				temp += data_x[i][j]; 
				temp2 += data_x[i][j]*data_x[i][j]; 
			} 
			sigma=sqrt(temp2/data_num-(temp/data_num)*(temp/data_num)); 
			temp=exp(1.0/(1.0*dim+4.0)*log(4.0/(1.0*dim+2.0))); 
			h[j]=temp*sigma*exp(-1.0/(dim+4)*log(1.0*data_num)); 
		} 
		ins=fopen(ap_filename,"a"); 
		fprintf(ins,"random seed=%d\n",stime); 
		fprintf(ins,"Normal Reference Rule:\n"); 
		for(j=1;j<=dim;j++) fprintf(ins,"%16.8f\n",h[j]); 
		fclose(ins); 

		//initial values using NRR
		x[0]=log(0.5);  
		for(j=1;j<=dim;j++) x[j]=log(h[j]); 
		x[-1]=cost(x); 
		printf("cost = %f\n",x[-1]); 
              printf("cost2 = %f\n", cost2(x));
		printf("Initialization finished\n"); 

		sizep[1]=0.75; 
		sizep[2]=0.80; 

		for(k=1;k<=warm;k++) np_gibbs(x);
		//printf("warm finished\n"); 

              for(k=0;k<=dim;k++) printf("para=%f\n", x[k]);

		for(i=-2;i<=dim;i++) 
		{ 
			sum_h[i]=0.0; 
			for(j=1;j<=num_batch;j++) batch_h[j][i]=0.0; 
			total_sd[i]=0.0; 
		} 
		for(k=1;k<=M;k++) 
		{ 
			np_gibbs(x);  

			sn=ceil(1.0*k/size_batch); 
			x[-2]=like(x); 
			plike[k]=x[-2]; 
			for(i=-2;i<=-1;i++) 
			{ 
				temp=x[i]; 
				sum_h[i] += temp; 
				batch_h[sn][i] += temp; 
				total_sd[i] += temp*temp; 
			} 
			temp=sqrt(x[0]); 
			sum_h[0] += temp; 
			batch_h[sn][0] += temp; 
			total_sd[0] += temp*temp; 
			for(j=1;j<=dim;j++) 
			{ 
				temp=exp(x[j]); 
				sum_h[j] += temp; 
				batch_h[sn][j] += temp; 
				total_sd[j] += temp*temp; 
			} 
			index=ceil(1.0*k/step);  
			data_post[index][0]=sqrt(x[0]); 
			for(j=1;j<=dim;j++) data_post[index][j]=exp(x[j]);
		} 
		//compute marginal likelihood 
              ik=1;
		for(i=-2;i<=dim;i++) 
		{ 
			std_h[i]=0.0; 
			sum_h[i]=sum_h[i]/M; 
			for(j=1;j<=num_batch;j++) 
			{ 
				temp=batch_h[j][i]/size_batch-sum_h[i]; 
				std_h[i] += temp*temp; 
			} 
			std_h[i]=sqrt(std_h[i]/(num_batch*num_batch-num_batch)); 
			total_sd[i]=sqrt(1.0*M*(total_sd[i]/M-sum_h[i]*sum_h[i])/(M-1)); 
			h[i]=sum_h[i]; 
			simu_res[ik][i]=h[i]; 
		}  
		for(i=-2;i<=dim;i++) 
		{ 
			sif[i]=std_h[i]*std_h[i]/(total_sd[i]*total_sd[i])*M; 
		} 
		simu_res[ik][dim+1]=1.0*accept_h/total_h; 
		ins=fopen(ap_filename,"a"); 
		fprintf(ins,"Block-wise Metropolis-Hastings Algorithm\n"); 
		fprintf(ins,"para: accept_rate=%g\n",1.0*accept_h/total_h); 
		fprintf(ins,"      average      batch_mean_sd  total_sd      SIF\n"); 
		for(i=-2;i<=-1;i++) 
		{ 
			fprintf(ins,"%d %12.2f %13.6f %13.6f %9.2f",-i,sum_h[i],std_h[i],total_sd[i],sif[i]); 
			fprintf(ins,"\n"); 
		} 
		for(i=0;i<=dim;i++) 
		{ 
			fprintf(ins,"%d %12.6f  %12.6f  %12.6f  %8.2f",i,sum_h[i],std_h[i],total_sd[i],sif[i]); 
			fprintf(ins,"\n"); 
		} 
		for(temp=0.0,k=1;k<=M;k++) 
		{ 
			plike[k] -= sum_h[-2]; 
			temp += 1.0/exp(plike[k]); 
		} 
		fclose(ins);  
		temp=temp/M; 
		temp = log(1.0/temp)+sum_h[-2]; 
		//harmonic mean of likelihood values 
		simu_res[ik][dim+2]=temp; 
		ins=fopen(ap_filename,"a"); 
		fprintf(ins,"log marginal likelihood=%g\n",temp); 
		fclose(ins); 
		//compute marginal likelihood 
		logmargin=loglikelihood(h);
		printf("logmargin=%f\n",logmargin);
		
		temp=logpriors(h);
		printf("logprior=%f\n",temp);
		
		temp2=logdensity(h,data_post); 
		printf("logdensity=%f\n",temp2);
		
		logmargin += temp-temp2; 
		simu_res[ik][dim+3]=logmargin; 
		ins=fopen(ap_filename,"a"); 
		fprintf(ins,"Chib's (1995, JASA) marginal likelihood=%g\n",logmargin); 
		fclose(ins); 
		
		//compute R2 
		kern(fit,h); 
		ins=fopen(ap_filename,"a"); 
		fprintf(ins,"R2=%6.4f\n",fit[1]); 
		fprintf(ins,"Mean Squared Error=%6.4f\n",fit[2]); 
		fprintf(ins,"\n"); 
		fclose(ins); 

		simu_res[ik][dim+4]=fit[1]; 
		simu_res[ik][dim+5]=fit[2]; 
		ins=fopen(ar_filename,"a"); 
		for(i=-2;i<=dim+5;i++) fprintf(ins,"%12.6f ",simu_res[ik][i]); 
		fprintf(ins,"\n"); 
		fclose(ins); 
	
	free_dvector(ave,0,dim); 
	free_dmatrix(cov,0,dim,0,dim); 
	free_dvector(x,-2,dim); 
	free_dvector(h,-2,dim); 
	free_dmatrix(batch_h,1,num_batch,-2,dim); 
	free_dvector(total_sd,-2,dim); 
	free_dvector(sif,-2,dim); 
	free_dvector(sum_h,-2,dim); 
	free_dvector(std_h,-2,dim); 
	free_dmatrix(data_yx,1,data_num,1,dim+1); 
	free_dvector(plike,1,M); 
	free_dmatrix(data_post,1,length,0,dim); 
	free_dmatrix(simu_res,1,mc,-2,dim+5); 
	free_dvector(fit,1,2); 

	//printf("end of main program\n"); 
	return 0; 
} 

