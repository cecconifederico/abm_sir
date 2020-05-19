#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>

#define SIR 0                    //if SIR=1 --> SIR model (tau_s=\infty), otherwise SIRS 

#define T_MAX 5000               //dynamics duration time
#define delta 0.00005            //integration step

#define S0 0.99                  //initial fraction of susceptible people

#define BETA 0.50                //infection rate
#define GAMMA 0.50               //recovering rate

#define tau_s 100                //average time to become susceptible again

int main()
{

  FILE *fa,*fb,*fc,*fd;
  unsigned long int t,tau,END;
  double susc,inf,rec,beta,gamma,omega;
  double susc0,inf0,rec0,susc1,inf1,rec1;

  END=1.0/(1.0*delta);
  beta=1.0*BETA;
  gamma=1.0*GAMMA;
  if (SIR==1)
    omega=0.0;
  else
    omega=1.0/(1.0*tau_s);

  susc=1.0*S0;
  inf=1.0-susc;
  rec=0.0;
  
  if ((fa=fopen("0_susceptibles.dat","w"))==NULL){
    printf ("Impossibile aprire il file-a \n");
    exit (1);
  }
  if ((fb=fopen("0_infected.dat","w"))==NULL){
    printf ("Impossibile aprire il file-b \n");
    exit (1);
  }
  if ((fc=fopen("0_recovered.dat","w"))==NULL){
    printf ("Impossibile aprire il file-c \n");
    exit (1);
  }
  if ((fd=fopen("0_susc_vs_inf.dat","w"))==NULL){
    printf ("Impossibile aprire il file-d \n");
    exit (1);
  }
  
  fprintf(fa,"0\t%g\n",susc);
  fprintf(fb,"0\t%g\n",inf);
  fprintf(fc,"0\t%g\n",rec);
  fprintf(fd,"%g\t%g\n",susc,inf);
  for (t=0;t<T_MAX;t++){
    for (tau=0;tau<END;tau++){
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

    fprintf(fa,"%lu\t%g\n",t+1,susc);
    fprintf(fb,"%lu\t%g\n",t+1,inf);
    fprintf(fc,"%lu\t%g\n",t+1,rec);
    fprintf(fd,"%g\t%g\n",susc,inf);
    
  }

  fclose(fa);
  fclose(fb);
  fclose(fc);
  fclose(fd);
}
