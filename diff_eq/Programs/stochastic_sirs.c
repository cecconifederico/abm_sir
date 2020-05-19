#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<time.h>
#include"ran2.h"

#define SIR 0                    //if SIR=1 --> SIR model (tau_s=\infty), otherwise SIRS 

#define T_MAX 5000               //dynamics duration time
#define delta 0.00005            //integration step
#define S 200                    //number of independent realizations

#define S0 0.99                  //initial fraction of susceptible people

#define BETA 0.50                //institutional infection rate
#define GAMMA 0.50               //recovering rate

#define uncompl 0.001            //probability to disobey the law

#define tau_s 100                //average time to become susceptible again
#define N 2000                   //maximum number of people which non-compliers can infect

int main()
{

  FILE *fa,*fb,*fc,*fd,*fe,*ff,*fg,*fh;
  unsigned long int t,tt,tau,END;
  int s;
  double susc,inf,rec,beta0,beta,gamma,omega,aux;
  double susc0,inf0,rec0,susc1,inf1,rec1;
  double su[T_MAX+1],in[T_MAX+1],re[T_MAX+1];
  time_t Time;  
  long a2,*a;

  a2=-(unsigned)time(&Time);
  a=&a2;

  END=1.0/(1.0*delta);
  beta0=1.0*BETA;
  gamma=1.0*GAMMA;
  if (SIR==1)
    omega=0.0;
  else
    omega=1.0/(1.0*tau_s);

  for (tt=0;tt<=T_MAX;tt++){
    su[tt]=0.0;
    in[tt]=0.0;
    re[tt]=0.0;
  }

  for (s=0;s<S;s++){
    susc=1.0*S0;
    inf=1.0-susc;
    rec=0.0;

    if (s==0){
      if ((fa=fopen("susceptibles_single.dat","w"))==NULL){
	printf ("Impossibile aprire il file-a \n");
	exit (1);
      }
      if ((fb=fopen("infected_single.dat","w"))==NULL){
	printf ("Impossibile aprire il file-b \n");
	exit (1);
      }
      if ((fc=fopen("recovered_single.dat","w"))==NULL){
	printf ("Impossibile aprire il file-c \n");
	exit (1);
      }
      if ((fd=fopen("susc_vs_inf_single.dat","w"))==NULL){
	printf ("Impossibile aprire il file-d \n");
	exit (1);
      }
  
      fprintf(fa,"0\t%g\n",susc);
      fprintf(fb,"0\t%g\n",inf);
      fprintf(fc,"0\t%g\n",rec);
      fprintf(fd,"%g\t%g\n",susc,inf);
    }
    su[0]=su[0]+1.0*susc;
    in[0]=in[0]+1.0*inf;
    re[0]=re[0]+1.0*rec;

    for (t=0;t<T_MAX;t++){
      for (tau=0;tau<END;tau++){

	aux=1.0*ran2(a);
	if (aux>=uncompl)
	  beta=1.0*beta0;
	else
	  beta=1.0*beta0+(N-1.0)*ran2(a);
      
	susc0=1.0*susc+delta*(omega*rec-beta*susc*inf);
	inf0=1.0*inf+delta*(beta*susc*inf-gamma*inf);
	rec0=1.0*rec+delta*(gamma*inf-omega*rec);

	susc1=susc+0.5*delta*(omega*rec0-beta*susc0*inf0);
	inf1=inf+0.5*delta*(beta*susc0*inf0-gamma*inf0);
	rec1=rec+0.5*delta*(gamma*inf0-omega*rec0);

	susc=1.0*susc1;
	inf=1.0*inf1;
	rec=1.0*rec1;
      }

      if (s==0){
	fprintf(fa,"%lu\t%g\n",t+1,susc);
	fprintf(fb,"%lu\t%g\n",t+1,inf);
	fprintf(fc,"%lu\t%g\n",t+1,rec);
	fprintf(fd,"%g\t%g\n",susc,inf);
      }

      su[t+1]=su[t+1]+1.0*susc;
      in[t+1]=in[t+1]+1.0*inf;
      re[t+1]=re[t+1]+1.0*rec;

    }

    if (s==0){
      fclose(fa);
      fclose(fb);
      fclose(fc);
      fclose(fd);
    }

    if ((fe=fopen("Susceptibles.dat","w"))==NULL){
      printf ("Impossibile aprire il file-e \n");
      exit (1);
    }
    if ((ff=fopen("Infected.dat","w"))==NULL){
      printf ("Impossibile aprire il file-f \n");
      exit (1);
    }
    if ((fg=fopen("Recovered.dat","w"))==NULL){
      printf ("Impossibile aprire il file-g \n");
      exit (1);
    }
    if ((fh=fopen("Susc_vs_inf.dat","w"))==NULL){
      printf ("Impossibile aprire il file-h \n");
      exit (1);
    }

    fprintf(fe,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fe,"%lu\t%g\n",tt,su[tt]/(1.0+s));

    fprintf(ff,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(ff,"%lu\t%g\n",tt,in[tt]/(1.0+s));

    fprintf(fg,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fg,"%lu\t%g\n",tt,re[tt]/(1.0+s));

    fprintf(fh,"#realizzazioni=%d\n\n",s+1);
    for (tt=0;tt<=T_MAX;tt++)
      fprintf(fh,"%g\t%g\n",su[tt]/(1.0+s),in[tt]/(1.0+s));

    fclose(fe);
    fclose(ff);
    fclose(fg);
    fclose(fh);
    
  } //end cycle on s

  
}
