#include<stdio.h>   
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"ran2.h"

#define TOPOLOGY 0                 //if =1 --> one-dimensional; otherwise, mean-field
#define L 2000                     //population's size

#define BETA_0 0.05                //infection probability per meeting
#define eta 0.15                   //noise (at every meeting, with probability eta BETA=1, otherwise =BETA_0
#define GAMMA 0.015                 //recovery probability
#define TAUs 50                   //average re-infection time
#define SIR 0                      //if =1 --> SIR (TAUs=\infty)

#define I0 0.05                    //initial density of infected (S0=1-I0; R0=0)

#define T_MAX 2000                 //single realization duration
#define S 20000                     //number of independent realizations

int main()
{
  FILE *fa,*fb,*fc,*fd,*fe;
  int s,i,j;
  unsigned long int t,tt;
  int vp[L],vm[L];
  int *type;
  unsigned long int *sus,*inf,*rec;
  int extr;
  double extr2,BETA,extr3;
  time_t Time;  
  long a2,*a;
  //char filen1[64],filen2[64],filen3[64],filen4[64],filen5[64];
  
  a2=-(unsigned)time(&Time);//-2573150310;//-53029820617;//
  a=&a2;

  sus=(unsigned long int*)calloc(T_MAX+1,sizeof(unsigned long int)); //number of susceptibles at time t;
  inf=(unsigned long int*)calloc(T_MAX+1,sizeof(unsigned long int)); //number of infected at time t;
  rec=(unsigned long int*)calloc(T_MAX+1,sizeof(unsigned long int)); //number of recovered at time t;

  type=(int*)calloc(L,sizeof(int)); //=0 --> susceptible; =1 --> infected; =2 --> recovered

  //GENERAL INITIALIZATIONS
  for (t=0;t<=T_MAX;t++){
    sus[t]=0;
    inf[t]=0;
    rec[t]=0;
  }
for (i=0;i<L;i++){
  if (i>0)
    vm[i]=i-1;
  else
    vm[i]=L;
  if (i<L-1)
    vp[i]=i+1;
  else
    vp[i]=0;
 }

  for (s=0;s<S;s++){

    //INIZIALIZZAZIONI SINGOLA REALIZZAZIONE
    for (i=0;i<L;i++){
      extr2=1.0*ran2(a);
      if (extr2<I0)
	type[i]=1;
      else
	type[i]=0;
    }
    for (i=0;i<L;i++){
      if (type[i]==0)
	sus[0]++;
      else {
	if (type[i]==1)
	  inf[0]++;
	else
	  rec[0]++;
      }
    }

    for (t=0;t<T_MAX;t++){
      //QUI DINAMICA
      for (tt=0;tt<L;tt++){

	i=L*ran2(a);
	if (TOPOLOGY!=1){
	paperino:
	  j=L*ran2(a);
	  if (j==i)
	    goto paperino;
	}
	else {
	  extr2=ran2(a);
	  if (extr2<0.5)
	    j=vm[i];
	  else
	    j=vp[i];
	}

	extr3=1.0*ran2(a);
	if (extr3<eta)
	  BETA=1.0;
	else
	  BETA=1.0*BETA_0;
	if (type[i]==0){
	  if (type[j]==1){
	    extr2=1.0*ran2(a);
	    if (extr2<BETA)
	      type[i]=1;
	  }
	}
	if (type[i]==1){
	  if (type[j]==0){
	    extr2=1.0*ran2(a);
	    if (extr2<BETA)
	      type[j]=1;
	  }
	  extr2=1.0*ran2(a);
	  if (extr2<GAMMA)
	    type[i]=2;
	}
	if (type[i]==2){
	  extr2=1.0*ran2(a);
	  if (SIR!=1 && extr2<1.0/TAUs)
	    type[i]=0;
	}
      }//fine ciclo su tt
      
      for (i=0;i<L;i++){
	if (type[i]==0)
	  sus[t+1]++;
	else {
	  if (type[i]==1)
	    inf[t+1]++;
	  else
	    rec[t+1]++;
	}
      }
      
    }//fine ciclo su t

    //plottigat

    if ((fa=fopen("STOCH_susceptibles.dat","w"))==NULL){
      printf ("Impossibile aprire il file-a \n");
      exit (1);
    }
    if ((fb=fopen("STOCH_infected.dat","w"))==NULL){
      printf ("Impossibile aprire il file-b \n");
      exit (1);
    }
    if ((fc=fopen("STOCH_recovered.dat","w"))==NULL){
      printf ("Impossibile aprire il file-c \n");
      exit (1);
    }
    if ((fd=fopen("STOCH_susc_vs_inf.dat","w"))==NULL){
      printf ("Impossibile aprire il file-d \n");
      exit (1);
    }
    if ((fe=fopen("STOCH_susc_vs_rec.dat","w"))==NULL){
      printf ("Impossibile aprire il file-e \n");
      exit (1);
    }

    fprintf(fa,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fa,"%lu\t%g\n",tt,1.0*sus[tt]/(L*(1.0+s)));

    fprintf(fb,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fb,"%lu\t%g\n",tt,1.0*inf[tt]/(L*(1.0+s)));

    fprintf(fc,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fc,"%lu\t%g\n",tt,1.0*rec[tt]/(L*(1.0+s)));

    fprintf(fd,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fd,"%g\t%g\n",1.0*sus[tt]/(L*(1.0+s)),1.0*inf[tt]/(L*(1.0+s)));

    fprintf(fe,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fe,"%g\t%g\n",1.0*sus[tt]/(L*(1.0+s)),1.0*rec[tt]/(L*(1.0+s)));

    fclose(fa);
    fclose(fb);
    fclose(fc);
    fclose(fd);
    fclose(fe);

  }//fine ciclo su s

}
