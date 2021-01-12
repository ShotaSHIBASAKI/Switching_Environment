//
//  main.c
//  resource, toxin, and species model
//  resource  switching model 1
//  Four (NOT Two) environmental states and seasonal change (x1->x2->x3->x4->x1) with stochastic durations
//  Createdon 2. 11. 19
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include"MT.h"//Mersenne twister


//  Defining constant parameters
#define replicate 100000// number of replicate in each case
#define dilution 0.1// dilution rate of  ystem
#define Yx1 1.0// biomass yield of 1 for resource
#define Yz1 1.0// biomass yield of 1for toxin
#define Yx2 1.0// biomass yield of 2 for resource
#define Yz2 1.0// biomass yield of 2 fot toxin
#define x_init 150// initial resoruce amount
#define y_init 10// initial species number (for each species)
#define z_init 100//initial toxin amount
#define change 1// number of scenario for changing consumer 2's parameter value
#define switch_num 9// number of swtiching rate nu
#define x1 200//resource supply season 1
#define x2 150//resource supply season 2
#define x3 100//resource supply season 3
#define x4 50// resource supply season 4
#define zmax 200 // max toxin supply
#define zmin 50 // min toxin supply
#define zmean (zmax+zmin)/2// mean toxin supply
#define d_num 10// variation of toxicity

double Uniform( void );
double rand_exp( double lambda );
int *Gillespie(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu);
double SumTransition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu);
int *Transition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu);
char fname[256];
char newdir[256];
char dirname[256];


/*main function*/
int main(void) {
    init_genrand((unsigned)time(NULL));
    /**Define parameter of consumer 1*/
    double mu1=1.0; // maximum growth rate of consumer 1
    double Kx1=100; // median-effect ressource amount for consumer 1
    double d1;//death rate of consumer 1
    double Kz1=100;//median-effect toxin amount
    /*Define parameter of consuemr 2*/
    double mu2=mu1; // maximum growth rate of consumer 2
    double Kx2=Kx1; // median-effect resource amount for consumer 2
    double d2;// death rate of consumer 2
    double Kz2=Kz1;//median-effect toxin amount
    
    /*Parameters used in simulations*/
    double nu;// switching rate
    int i;// index for switching rate
    int j;// index for replicate
    int k;// index for parameter value of consumer 2
    int l;// index of death rate array
    double nu_array[switch_num];//nu array
    nu_array[0]=pow(10, -5);
    nu_array[1]=pow(10, -4);
    nu_array[2]=pow(10, -3);
    nu_array[3]=pow(10, -2);
    nu_array[4]=pow(10, -1);
    nu_array[5]=pow(10, 0);
    nu_array[6]=pow(10, 1);
    nu_array[7]=pow(10, 2);
    nu_array[8]=pow(10, 3);
    int model=1; // specify the model of environmental switching
    int *data;// result of G algotithm
    int type;/*model type. 0: no consumer 2, 1: cnosumer 2 with smaller mu
                    2: cnosumer 2 with larger K, 3: cnosumer 2 with larger d */
    int result[5];// extinction time
    double one_extinct1[switch_num]; // extinction prob of species 1 in one consumer model
    double two_extinct1[change][switch_num]; // extinction prob of species 1 in two consumer model
    double both_extinct[change][switch_num]; // both species go extinct
    double comp_exclude1[change][switch_num]; // only species 1 goes extinct in two consumer model
    double coexist[change][switch_num];// conexistence of 2 species
    /*Simulation starts*/
    for (l=0;l<d_num;l++)
    {
        d1=0.1*(l+1);// define death rate of species 1
        d2=d1; // death rate of species 2 should be the same except for case where speceis 2 has larger death rate tna species 1
        // move directory
        sprintf(dirname,"./death%.1f",d1);
        mkdir(dirname, 0700); // create a folder
        chdir(dirname);
        sprintf(dirname,"./Model%d",model);
        mkdir(dirname, 0700); // create a folder
        chdir(dirname);
        printf("Death rate%.1f of model %d \n", d1, model);
        mkdir("./OneConsumer", 0700); // create a folder
        mkdir("./TwoConsumer", 0700); // create a folder
        for (type=0; type<2;type++)
        {
            if(type==0)
            {
                printf("Start simulation in absence of species2 \n");
                // consumer 1 alone
                for (i=0;i<switch_num;i++)
                {
                    nu=nu_array[i];
                    one_extinct1[i]=0.0; // extinction prob of species 1 in one consumer model
                    
                    for (j=0;j<replicate;j++)
                    {
                        //run Gillispie algorithm and save extinction time
                        data=Gillespie(x_init, y_init, 0, z_init, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu);
                        result[0]=*(data+0);//resource
                        result[1]=*(data+1);//consumer 1
                        result[2]=*(data+2);//consumer 2
                        result[3]=*(data+3);//toxin
                        result[4]=*(data+4);//environment
                        if (result[1]<=0)
                        {
                            one_extinct1[i]+=1;
                        }
                    }
                    one_extinct1[i]=one_extinct1[i]/replicate;
                }
               
                
                
                FILE *fp;
                sprintf(fname,"OneConsumer_Extinction_model%d.csv",model);
                fp=fopen(fname, "w");
                fprintf(fp, "nu=10^-5, nu=10^-4, nu=10^-3, nu=10^-2, nu=10^-1m nu=1, nu=10, nu=10^2, nu=10^3\n"); // header
                for(i=0;i<switch_num-1;i++)
                {
                    fprintf(fp, "%.5f,", one_extinct1[i]);
                }
                fprintf(fp, "%.5f\n", one_extinct1[switch_num-1]);
                fclose(fp);
                sprintf(newdir, "./OneConsumer/OneConsumer_Extinction_model%d.csv",model);
                rename(fname, newdir);// move file
            }
            else
            {
                printf("Start simulation competing with species2. Type %d \n", type);
                sprintf(newdir, "./TwoConsumer/type%d", type);
                mkdir(newdir, 0700); // create a folder
                //vs consumer2
                double change_array[change];
                // to change the parameter value of consumer 2
                change_array[0]=1.1;
                //Co-culture with species 2/
                for(k=0;k<change;k++)
                {
                    //Change parameter value of consumer 2
                    if (type==1)
                    {
                        mu2=mu1/change_array[k];// smaller growth rate
                    }
                    else if(type==2)
                    {
                        Kx2=Kx1*change_array[k];// larger median-effect resoruce amount
                    }
                    else if(type==3)
                    {
                        d2=d1*change_array[k];// larger death rate
                    }
                    else
                    {
                        Kz2=Kz1/change_array[k];//smaller median-effect toxi amount
                    }
                    
                    for (i=0;i<switch_num;i++)
                    {
                        nu=nu_array[i];
                        two_extinct1[k][i]=0.0; // extinction prob of species 1 in two consumer model
                        both_extinct[k][i]=0.0; // both species go extinct
                        comp_exclude1[k][i]=0.0; // only species 1 goes extinct in two consumer model
                        coexist[k][i]=0.0; // coexistence of 2 species
                        for(j=0;j<replicate;j++)
                        {
                            //run Gillispie algorithm and save extinction time
                            data=Gillespie(x_init, y_init, y_init, z_init, mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu);
                            result[0]=*(data+0);//resource
                            result[1]=*(data+1);//consumer 1
                            result[2]=*(data+2);//consumer 2
                            result[3]=*(data+3);//toxin
                            result[4]=*(data+4);//environment
                            if (result[1]<=0)
                            {
                                two_extinct1[k][i]+=1;
                                if (result[2]<=0)
                                {
                                    both_extinct[k][i]+=1;
                                }
                                else
                                {
                                    comp_exclude1[k][i]+=1;
                                }
                            }
                            else
                            {
                                if (result[2]>0)
                                {
                                    coexist[k][i]+=1;
                                }
                            }
                        }
                        two_extinct1[k][i]=two_extinct1[k][i]/replicate;
                        both_extinct[k][i]=both_extinct[k][i]/replicate;
                        comp_exclude1[k][i]=comp_exclude1[k][i]/replicate;
                        coexist[k][i]=coexist[k][i]/replicate;
                    }
                }
                //save csv files
                //extinction prob of species 1
                FILE *fp1;
                sprintf(fname,"TwoConsumer_Extinction_model%d_competitor%d.csv",model, type);
                fp1=fopen(fname, "w");
                fprintf(fp1, "nu=10^-5, nu=10^-4, nu=10^-3, nu=10^-2, nu=10^-1m nu=1, nu=10, nu=10^2, nu=10^3\n"); // header
                for(k=0;k<change;k++)
                {
                    for (i=0;i<switch_num-1;i++)
                    {
                        fprintf(fp1, "%.5f,", two_extinct1[k][i]);
                    }
                    fprintf(fp1, "%.5f\n", two_extinct1[k][switch_num-1]);
                }
                fclose(fp1);
                sprintf(newdir, "./TwoConsumer/type%d/TwoConsumer_Extinction_model%d_competitor%d.csv", type,model,type);
                rename(fname, newdir);// move file
                
                //prob of both extinction
                FILE *fp2;
                sprintf(fname,"TwoConsumer_BothExtinction_model%d_competitor%d.csv",model, type);
                fp2=fopen(fname, "w");
                fprintf(fp2, "nu=10^-5, nu=10^-4, nu=10^-3, nu=10^-2, nu=10^-1m nu=1, nu=10, nu=10^2, nu=10^3\n"); // header
                for(k=0;k<change;k++)
                {
                    for (i=0;i<switch_num-1;i++)
                    {
                        fprintf(fp2, "%.5f,", both_extinct[k][i]);
                    }
                    fprintf(fp2, "%.5f\n", both_extinct[k][switch_num-1]);
                }
                fclose(fp2);
                sprintf(newdir, "./TwoConsumer/type%d/TwoConsumer_BothExtinction_model%d_competitor%d.csv", type,model,type);
                rename(fname, newdir);// move file
                
                // competitive exclusion
                FILE *fp3;
                sprintf(fname,"TwoConsumer_CompetitiveExclusion_model%d_competitor%d.csv",model, type);
                fp3=fopen(fname, "w");
                fprintf(fp3, "nu=10^-5, nu=10^-4, nu=10^-3, nu=10^-2, nu=10^-1m nu=1, nu=10, nu=10^2, nu=10^3\n"); // header
                for(k=0;k<change;k++)
                {
                    for (i=0;i<switch_num-1;i++)
                    {
                        fprintf(fp3, "%.5f,", comp_exclude1[k][i]);
                    }
                    fprintf(fp3, "%.5f\n", comp_exclude1[k][switch_num-1]);
                }
                fclose(fp3);
                sprintf(newdir, "./TwoConsumer/type%d/TwoConsumer_CompetitiveExclusion_model%d_competitor%d.csv", type,model,type);
                rename(fname, newdir);// move file
                
                //coexistence
                FILE *fp4;
                sprintf(fname,"TwoConsumer_Coexistence_model%d_competitor%d.csv",model, type);
                fp4=fopen(fname, "w");
                fprintf(fp4, "nu=10^-5, nu=10^-4, nu=10^-3, nu=10^-2, nu=10^-1m nu=1, nu=10, nu=10^2, nu=10^3\n"); // header
                for(k=0;k<change;k++)
                {
                    for (i=0;i<switch_num-1;i++)
                    {
                        fprintf(fp4, "%.5f,", coexist[k][i]);
                    }
                    fprintf(fp4, "%.5f\n", coexist[k][switch_num-1]);
                }
                fclose(fp4);
                sprintf(newdir, "./TwoConsumer/type%d/TwoConsumer_Coexistence_model%d_competitor%d.csv", type,model,type);
                rename(fname, newdir);// move file
                
            }
        }
        chdir("../../"); // return to the upper dirctory
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
int *Gillespie(int x0, int y10, int y20, int z0, double mu1, double Kx1, double d1, double Kz1, double mu2, double Kx2, double d2, double Kz2, double nu)
{
    int x=x0;// amount of resource
    int y1=y10;//amount of consumer 1
    int y2=y20;//amount of consumer 2
    int z=z0;// amount of toxin
    int xi;//environmental state (xi=1; rich xi=-1 poor resource supply)
    int *result;// result of transition
    double p=genrand_real3();
    if(p<0.2)
    {
        xi=1;// supply x1
    }
    
     else if(p<0.5)
    {
        xi=2;//supply s2
    }
    else if(p<0.75)
    {
        xi=3;//supply x3
    }
    else
    {
        xi=-1;
    }
    double t=0;// index for time
    double T_end = 200.0;//end of the simulation in coutinuous time scale
    int flag2=0;
    double waiting_time;//waiting time until next event happens
    double lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1,Kz1,  mu2, Kx2, d2, Kz2, nu);// sum of transtion rates over all events
    waiting_time=rand_exp(lambda);
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
        lambda=SumTransition(x, y1, y2, z, xi,mu1, Kx1, d1, Kz1, mu2, Kx2, d2, Kz2, nu);
        waiting_time=rand_exp(lambda);
    }
    // write csv files of final states: for debugging
    return result;
}

double SumTransition(int x, int y1, int y2, int z, int xi,  double mu1, double Kx1, double d1, double Kz1, double mu2,double Kx2 ,double d2, double Kz2, double nu)
{
    /*Sum of transition rate for calculating waiting time*/
    double prob[9];
    double sum;
    int i;
   if (xi==1)
   {
       prob[0]=dilution*x1;// increase of resoruce with x1 supply
   }
   else if(xi==2)
   {
       prob[0]=dilution*x2;// increase of resoruce with x2 supply
   }
   else if(xi==3)
   {
       prob[0]=dilution*x3;// increase of resoruce with x3 supply
   }
   else
   {
        prob[0]=dilution*x4;// increase of resoruce with x4 supply
   }
   prob[1]=dilution*x+mu1/Yx1*x/(x+Kx1)*y1+mu2/Yx2*x/(x+Kx2)*y2;// decrease of  resource
   prob[2]=mu1*x/(x+Kx1)*y1;// increase of consumer 1
   prob[3]=(dilution+d1*z/(z+Kz1))*y1;// decrease of consumer 1
   prob[4]=mu2*x/(x+Kx2)*y2;// increase of consumer 2
   prob[5]=(dilution+d2*z/(z+Kz2))*y2;// decrease of consumer 2
   prob[6]=dilution*zmean;// increse of toxin with  mean supply
   prob[7]=dilution*z+d1/Yz1*z/(z+Kz1)*y1+d2/Yz2*z/(z+Kz2)*y2;//decreasing toxin
   prob[8]=nu;//switching rateecrease of consumer 2
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
    if (xi==1)
    {
        prob[0]=dilution*x1;// increase of resoruce with x1 supply
    }
    else if(xi==2)
    {
        prob[0]=dilution*x2;// increase of resoruce with x2 supply
    }
    else if(xi==3)
    {
        prob[0]=dilution*x3;// increase of resoruce with x2 supply
    }
     else
     {
         prob[0]=dilution*x4;// increase of resoruce with x4 supply
     }
   prob[1]=dilution*x+mu1/Yx1*x/(x+Kx1)*y1+mu2/Yx2*x/(x+Kx2)*y2;// decrease of  resource
   prob[2]=mu1*x/(x+Kx1)*y1;// increase of consumer 1
   prob[3]=(dilution+d1*z/(z+Kz1))*y1;// decrease of consumer 1
   prob[4]=mu2*x/(x+Kx2)*y2;// increase of consumer 2
   prob[5]=(dilution+d2*z/(z+Kz2))*y2;// decrease of consumer 2
   prob[6]=dilution*zmean;// increse of toxin with  mean supply
   prob[7]=dilution*z+d1/Yz1*z/(z+Kz1)*y1+d2/Yz2*z/(z+Kz2)*y2;//decreasing toxin
   prob[8]=nu;//switching resource
   sum=0.0;
    for(i=0;i<9;i++)
    {
        sum+=prob[i];
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
        //seasonal environmental change
        if(result[4]==4)
        {
            result[4]=1;
            
        }
        else
        {
            result[4]+=1;
        }
    }
    return result;
}
