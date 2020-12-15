//
//  main.c
//  resource, toxin, and consumer model
//  only resource supply switching (model 1)
//  check the convergence to quasi-stationary distribution
//  Created by Shota on 14. 11. 19
//  Copyright © 2019 ShotaS. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include"MT.h"//Mersenne twister


//  Defining constant parameters
#define replicate 1000// number of replicate in each case
#define dilution 0.1// dilution rate of  ystem
#define Yx1 1.0// biomass yield of 1 for resource
#define Yz1 1.0// biomass yield of 1for toxin
#define Yx2 1.0// biomass yield of 2 for resource
#define Yz2 1.0// biomass yield of 2 fot toxin
#define x_init 150// initial resoruce amount
#define y_init 10// initial species number (for each species)
#define z_init 125//initial toxin amount
#define change 1// number of scenario for changing consumer 2's parameter value
#define switch_num 7// number of swtiching rate nu
#define xmax 200// resource supply at rich
#define xmin 50// resource supply at poor
#define zmean 125// mean toxin supply

double Uniform( void );
double rand_exp( double lambda );
int *Gillespie1(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu, int tri);
int *Gillespie2(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu, int model, int tri);
double SumTransition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu);
int *Transition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu);
char fname[256];


/*main function*/
int main(void) {
    init_genrand((unsigned)time(NULL));// generate seed for random var.
    /**Define parameter of consumer 1*/
    double mu1=1.0; // maximum growth rate of consumer 1
    double Kx1=100; // median-effect ressource amount for consumer 1
    double d1=0.2;// max death rate of consumer 1
    double Kz1=100;//median-effect toxin amount
    /*Define parameter of consuemr 2*/
    double mu2=mu1; // maximum growth rate of consumer 2
    double Kx2=Kx1; // median-effect resource amount for consumer 2
    double d2=d1;// death rate of consumer 2
    double Kz2=Kz1;//median-effect toxin amount
    
    /*Parameters used in simulations*/
    double nu;// switching rate
    int i;// index for switching rate
    int j;// index for replicate
    int k;// index for parameter value of consumer 2
    double nu_array[switch_num];//nu array
    nu_array[0]=pow(10, -5);
    nu_array[1]=pow(10, -4);
    nu_array[2]=pow(10, -3);
    nu_array[3]=pow(10, -2);
    nu_array[4]=pow(10, -1);
    nu_array[5]=pow(10, 0);
    nu_array[6]=pow(10, 1);
    int model;/*model type. 0: no consumer 2, 1: cnosumer 2 with smaller mu
                    2: cnosumer 2 with larger K, 3: cnosumer 2 with larger d */
    printf("Please choose a model; model type 0: no consumer 2, 1: cnosumer 2 with smaller mu 2: cnosumer 2 with larger Kx, 3: cnosumer 2 with larger d, 4: consumer 2 with smaller Kz\n --->");
    scanf("%d", &model);
    if(model>4 )
    {
        printf("model should be 0, 1, 2, 3, or 4.\n");
        return 1;
    }
    
    /*Simulation starts*/
    if(model==0)
    {
        // consumer 1 alone
        for (i=0;i<switch_num;i++)
        {
            nu=nu_array[i];
            printf("switching rate nu= 10^ %.0f \n", log10(nu));
            for (j=0;j<replicate;j++)
            {
                printf("Troal %d \n", j);
                FILE *fp;
                sprintf(fname,"Convergence_RTC_switching10^%.0f_trial%d.csv",log10(nu), j);
                fp=fopen(fname, "w");
                fprintf(fp, "time, resource, consumer1, toxin, environment\n");
                fclose(fp);
                /*run Gillispie algorithm and save extinction time*/
                Gillespie1(x_init, y_init, 0, z_init, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu,j);
            }
        }
    }
    else
    {
        //vs consumer2
        double change_array[change];
        // to change the parameter value of consumer 2
        change_array[0]=1.1;

        /*Co-culture with consumer 2*/
        for(k=0;k<change;k++)
        {
            //Change parameter value of consumer 2
            if (model==1)
            {
                mu2=mu1/change_array[k];// smaller growth rate
            }
            else if(model==2)
            {
                Kx2=Kx1*change_array[k];// larger median-effect resoruce amount
            }
            else if(model==3)
            {
                d2=d1*change_array[k];// larger death rate
            }
            else
            {
                Kz2=Kz1/change_array[k];//smaller median-effect toxi amount
            }
            /* consumer 1 vs consuemr 2*/
	    //for (i=0;i<num_sim;i++)
            for (i=4;i<5;i++)
            {
                nu=nu_array[i];
                for(j=0;j<replicate;j++)
                {
                    FILE *fp;
                    sprintf(fname,"Convergence_RTC_model%d_switching10^%.0f_trial%d.csv",model, log10(nu), j);
                    fp=fopen(fname, "w");
                    fprintf(fp, "time, resource, consumer1, consumer2, toxin, environment\n");
                    fclose(fp);
                    
                    /*run Gillispie algorithm and save extinction time*/
                    Gillespie2(x_init, y_init, y_init, z_init, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu, model, j);
                }
            }
        }
    }
    return 0;
}

/*Unifrom distribution*/
double Uniform( void )
{
    return genrand_real3();
}
/*Exponential distribution*/
double rand_exp( double lambda )
{
   return -log(Uniform())/lambda;
}

/*Gillespie alg*/
int *Gillespie1(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu, int tri)
{
    int x=x0;// amount of resource
    int y1=y10;//amount of consumer 1
    int y2=y20;//amount of consumer 2
    int z=z0;// amount of toxin
    int xi;//environmental state (xi=1; rich xi=-1 poor resource supply)
    int *result;// result of transition
    if (genrand_real3()<0.5)
    {
        xi=1;// aboundant resource
    }
    else
    {
        xi=-1;//scarece resource supply
    }
    double t=0;// index for time
    double T_end = 200.0;//end of the simulation in coutinuous time scale
    
    double waiting_time;//waiting time until next event happens
    double lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1,Kz1,  mu2, Kx2, d2, Kz2, nu);// sum of transtion rates over all events
    waiting_time=rand_exp(lambda);
    FILE *fp;
    sprintf(fname,"Convergence_RTC_switching10^%.0f_trial%d.csv",log10(nu), tri);
    fp=fopen(fname, "a");
    fprintf(fp, "%.3f, %d, %d,  %d, %d\n", 0.0, x_init, y_init, z_init, xi);
    while(t+waiting_time<T_end)
    {
        result=Transition(x, y1, y2, z, xi, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu );
        // Update process
        x=*(result+0);
        y1=*(result+1);
        y2=*(result+2);
        z=*(result+3);
        xi=*(result+4);
        t+=waiting_time;
        //write csv file
        fprintf(fp, "%.3f, %d, %d, %d, %d\n", t, x, y1, z, xi);
        lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu);
        waiting_time=rand_exp(lambda);
    }
    fclose(fp);
    return result;
}
int *Gillespie2(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu, int model, int tri)
{
    int x=x0;// amount of resource
    int y1=y10;//amount of consumer 1
    int y2=y20;//amount of consumer 2
    int z=z0;// amount of toxin
    int xi;//environmental state (xi=1; rich xi=-1 poor resource supply)
    int *result;// result of transition
    if (genrand_real3()<0.5)
    {
        xi=1;// aboundant resource
    }
    else
    {
        xi=-1;//scarece resource supply
    }
    double t=0;// index for time
    double T_end = 200.0;//end of the simulation in coutinuous time scale
    int flag2=0;
    double waiting_time;//waiting time until next event happens
    double lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1,Kz1,  mu2, Kx2, d2, Kz2, nu);// sum of transition rates over all events
    waiting_time=rand_exp(lambda);
    FILE *fp;
    sprintf(fname,"Convergence_RTC_model%d_switching10^%.0f_trial%d.csv",model, log10(nu), tri);
    fp=fopen(fname, "a");
    fprintf(fp, "%.3f, %d, %d, %d, %d, %d\n", 0.0, x_init, y_init, y_init, z_init, xi);
    while(t+waiting_time<T_end)
    {
        result=Transition(x, y1, y2, z, xi, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu );
        // Update process
        x=*(result+0);
        y1=*(result+1);
        y2=*(result+2);
        z=*(result+3);
        xi=*(result+4);
        t+=waiting_time;
        //write csv file
        fprintf(fp, "%.3f, %d, %d, %d, %d, %d\n", t, x, y1, y2, z, xi);
        lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu);
        waiting_time=rand_exp(lambda);
    }
    fclose(fp);
    return result;
}
double SumTransition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu)
{
    /*Sum of transition rate for calculating waiting time*/
    double prob[9];
    double sum;
    int i;
    if(xi==1)
    {
        prob[0]=dilution*xmax;// increase of resoruce at rich
    }
    else
    {
        prob[0]=dilution*xmin;// increase of resource at poor
    }
    prob[1]=dilution*x+mu1/Yx1*x/(x+Kx1)*y1+mu2/Yx2*x/(x+Kx2)*y2;// decrease of  resource
    prob[2]=mu1*x/(x+Kx1)*y1;// increase of consumer 1
    prob[3]=(dilution+d1*z/(z+Kz1))*y1;// decrease of consumer 1
    prob[4]=mu2*x/(x+Kx2)*y2;// increase of consumer 2
    prob[5]=(dilution+d2*z/(z+Kz2))*y2;// decrease of consumer 2
    prob[6]=dilution*zmean;// increse of toxin
    prob[7]=dilution*z+d1/Yz1*z/(z+Kz1)*y1+d2/Yz2*z/(z+Kz2)*y2;//decreasing toxin
    prob[8]=nu;//switching rate
    sum=0.0;
    for(i=0;i<9;i++)
    {
        sum+=prob[i];
    }
    return sum;
}

int *Transition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu)
{
    /*Calculate which event happens*/
    double prob[9];
    double sum;
    static int result[5];// array of x, y1, y2, z, xi at certain time point
    result[0]=x;
    result[1]=y1;
    result[2]=y2;
    result[3]=z;
    result[4]=xi;
    //Calculate probability that each event happens
    int i;
    if(xi==1)
    {
        prob[0]=dilution*xmax;// increase of resoruce at rich
    }
    else
    {
        prob[0]=dilution*xmin;// increase of resource at poor
    }
    prob[1]=dilution*x+mu1/Yx1*x/(x+Kx1)*y1+mu2/Yx2*x/(x+Kx2)*y2;// decrease of  resource
    prob[2]=mu1*x/(x+Kx1)*y1;// increase of consumer 1
    prob[3]=(dilution+d1*z/(z+Kz1))*y1;// decrease of consumer 1
    prob[4]=mu2*x/(x+Kx2)*y2;// increase of consumer 2
    prob[5]=(dilution+d2*z/(z+Kz2))*y2;// decrease of consumer 2
    prob[6]=dilution*zmean;// increse of toxin
    prob[7]=dilution*z+d1/Yz1*z/(z+Kz1)*y1+d2/Yz2*z/(z+Kz2)*y2;//decreasing toxin
    prob[8]=nu;//switching rateecrease of consumer 2
    sum=0.0;
    for(i=0;i<9;i++)
    {
        sum+=prob[i];
    }
    if (sum<=0)
    {
        printf("Error in sum of rates over events");
    }
    for(i=0;i<9;i++)
    {
        prob[i]=prob[i]/sum;// normalize
    }
    // Determine which event happens
    int event=0;
    double p=prob[event];
    double q=genrand_real1();
    while(p<q)
    {
        event+=1;
        p+=prob[event];
    }
    if (event==0)
    {
        result[0]+=1;// increasing resource
    }
    else if (event==1)
    {
        result[0]-=1;//decreasing resource
    }
    else if (event==2)
    {
        result[1]+=1;//increasing consumer 1
    }
    else if (event==3)
    {
        result[1]-=1;//decreasing consumer 1
    }
    else if (event==4)
    {
        result[2]+=1;//increasing consumer 2
    }
    else if (event==5)
    {
        result[2]-=1;//decreasing consumer 2
    }
    else if(event==6)
    {
        result[3]+=1;//increase toxin
    }
    else if(event==7)
    {
        result[3]-=1;//decreasing toxin
    }
    else
    {
        result[4]=-result[4];//switching environment
    }
    return result;
}
