#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
#define pi 3.14159265
 
#define data_num 50    
#define dim 2

static int length=5; 
static int mc=1000; 
static double mutsizp=2; 
static double prior_p=2.0; 
static double prior_st=0.1; 

static double data_y[1+data_num]; 
static double ystar[1+data_num]; 
static double mhat[1+data_num]; 
static double m2[1+data_num]; 
static double data_x[1+data_num][dim+1]; 

static double xm[1+data_num]; 

static double sizep[dim+1]; 
static double accept_p[dim+1]; 
static double total_p=1.0; 
static double Gamma_table[1+10000]; 

static double total_h=1.0;
static double accept_h=0.0;


