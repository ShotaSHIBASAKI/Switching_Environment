//
//  parm_generator.c
//  Sample parameters for simulating communities
//  One can check whether this code has bugs or not  by
//  plotting  parameter values' distributions
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "MT.h"
#include <unistd.h>

void Parm_Gen(int num, double delta);
double Uniform( void );
double rand_normal( double mu, double sigma );
double rand_gamma( double theta, double kappa );
double rand_beta(double a, double b);
double *malloc_vector(int n);
void free_vector(double *a);
int *malloc_vector_int(int n);
void free_vector_int(int *a);
double **malloc_matrix(int nrow, int ncol);
void free_matrix(double **a, int nrow);

void Parm_Gen(int num, double delta)
{
    if (delta>=1.00)
    {
        delta=0.99;
    }
    int type = 100;  // sampleing number for parameter values
    int n, i, j, k, l;
    double **beta; //matrix for species per capita effects on compounds
    beta=malloc_matrix(num, num); // beta[i,j] is effect of species i on compund j
    double **K; // amount of compound that gives harf-max effect on species
    K=malloc_matrix(num, num); //
    char path[256];
    char fname1[256];
    char fname2[256];
    printf("Death rate %.2f Start\n", delta);
    for (i=1;i<=type;i++)
    {
        //printf("paramerter set %d \n", i);
        // generate parameter values
        for (k=1;k<=num; k++)
        {
            for(l=1;l<=num;l++)
            {
                K[k][l]=rand_normal(100,10); // rGaussin  whose mean is 100 and std=10
                if (l%2==1)
                {
                    // resource follows Gaussin
                    beta[k][l]=rand_normal(1.0, 0.1); // Gaussin  whose mean is 1.0 and std=0.1
                }
                else
                {
                    //toxin follows Beta dist
                    beta[k][l]=-rand_beta(100*delta, 100*(1-delta));
                }
            }
        }
        sprintf(path, "./parameter%d", i);
        mkdir(path, 0777);
        chdir(path);

        // write cev file for parameter list beta and K
        FILE *fp1;
        sprintf(fname1,"Prameter_beta_set%d.csv",i);
        if ((fp1 = fopen(fname1,"w"))== NULL)
        {
            printf("unable to make result file. \n");
            exit(1);
        }
        fp1 = fopen(fname1,"w");
        for (k=1; k<=num; k++)
        {
            for(l=1;l<=num-1;l++)
            {
                fprintf(fp1,"%.3f,",beta[k][l]);
            }
            if(k<num)
            {
                fprintf(fp1,"%.3f\n",beta[k][num]);
            }
            else
            {
                fprintf(fp1,"%.3f",beta[k][num]);
            }
        }
        fclose(fp1);
        
        FILE *fp2;
        sprintf(fname2,"Prameter_K_set%d.csv",i);
        if ((fp2 = fopen(fname2,"w"))== NULL)
        {
            printf("unable to make result file. \n");
            exit(1);
        }
        fp2 = fopen(fname2,"w");
        for (k=1; k<=num; k++)
        {
            for(l=1;l<=num-1;l++)
            {
                fprintf(fp2,"%.3f,",K[k][l]);
            }
           if(k<num)
           {
               fprintf(fp2,"%.3f\n",K[k][num]);
           }
           else
           {
               fprintf(fp2,"%.3f",K[k][num]);
           }
        }
        fclose(fp2);
        chdir("../");
    }
    free_matrix(beta, num);
    free_matrix(K, num);
}

//probabvility distributions
double Uniform( void )
{
    return genrand_real3();
}

double rand_normal( double mu, double sigma )
{
   double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
   return mu + sigma*z;
}

double rand_gamma( double theta, double kappa )
{

   int int_kappa;
   double frac_kappa;

   int_kappa  = (int)kappa;
   frac_kappa = kappa - (double)int_kappa;

   double u,uu;
   double b,p,x_frac,x_int;
   int i;

   /*integer part*/
   x_int=0;
   for(i=0;i<int_kappa;i++){
       x_int+=-log(Uniform()); // add expnential random number with mean 1
   }

   /*fractional part*/
   if( fabs(frac_kappa) < 0.01 ) x_frac=0;

   else{
       b=(exp(1.0)+frac_kappa)/exp(1.0);
       while(1){

           u=Uniform();
           p=b*u;

           uu=Uniform();

           if(p<=1.0){
               x_frac=pow(p,1.0/frac_kappa);
               if(uu<=exp(-x_frac)) break;
           }

           else{
               x_frac=-log((b-p)/frac_kappa);
               if(uu<=pow(x_frac,frac_kappa-1.0)) break;
           }

       }
   }

   return (x_int+x_frac)*theta;
}
double rand_beta(double a, double b)
{
    /*
     Use the transiformation below:
     Suppose X and Y follow above Gamma(theta, a) and  Gamma(theta, b) respectively,
     then the fraction Z=X/(X+Y) follows Beta distribution (a, b).
     */
    double theta=50.0; // value of theta does not matter
    double x,y;
    x=rand_gamma(theta, a);
    y=rand_gamma(theta, b);
    return x/(x+y);
}
/*vecotr ver*/
double *malloc_vector(int n)
{
    double *a;
    if ((a=malloc(sizeof(double)*n))==NULL)
    {
        printf("unable to keep memory for vector, sorry\n");
        exit(1);
    }
    return (a-1);/*row number strats　with 1*/
}
void free_vector(double *a)
/*free the memory*/
{
    free(a+1);
}

/*vecotr ver*/
int *malloc_vector_int(int n)
{
    int *a;
    if ((a=malloc(sizeof(double)*n))==NULL)
    {
        printf("unable to keep memory for vector, sorry\n");
        exit(1);
    }
    return (a-1);/*row number stratswith number 1*/
}
void free_vector_int(int *a)
/*free the memory*/
{
    free(a+1);
}

/*matrix ver*/
double **malloc_matrix(int nrow, int ncol)
{
    double **a;
    int i;
    if ((a=malloc(nrow*sizeof(double *)))==NULL)
    {
        printf("unable to keep memory for matrix sorry");
        exit(1);
    }
    a = a-1;
    for (i=1;i<=nrow;i++) a[i]=malloc(ncol*sizeof(double));
    for (i=1;i<=nrow;i++) a[i]=a[i]-1;/* move row direction*/
    
    return a;
}
void free_matrix(double **a, int nrow)
{
    int i;
    for ( i=1; i<=nrow; i++) free((void *)(a[i]+1));
    free((void *)(a+1));}

int  main(void){
    int *NUM;
    int count=5;
    int i,j;
    double *DEATH;
    char path[256];
    NUM=malloc_vector_int(count);//species number vector
    NUM[1]=2;
    NUM[2]=4;
    NUM[3]=6;
    NUM[4]=8;
    NUM[5]=10;
    DEATH=malloc_vector(count);// mean maximum death rate
    DEATH[1]=0.1;
    DEATH[2]=0.2;
    DEATH[3]=0.4;
    DEATH[4]=0.6;
    DEATH[5]=1.0;
    for(i=1;i<=count;i++)
    {
        //printf("Death rate %.1f Start\n", DEATH[i]);
        sprintf(path,"./death%.1f", DEATH[i]);
        mkdir(path, 0777);
        chdir(path);
        for (j=1;j<=count;j++)
        {
            //move directory
            printf("=> Number of species is %d\n", NUM[j]);
            sprintf(path,"./Species_number%d",NUM[j]);
            mkdir(path, 0777);
            chdir(path);
            Parm_Gen(NUM[j],DEATH[i]);// parm generation
            chdir("../");
        }
        chdir("../");
    }
    free_vector_int(NUM);
    free_vector(DEATH);
    return 0;
}
