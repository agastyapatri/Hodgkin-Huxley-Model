/****************************************************************************
 *                                                                          *
 * File    : main.c                                                         *
 *                                                                          *
 * Purpose : Console mode (command line) program.                           *
 *                                                                          *
 * History : Date      Reason                                               *
 *           00/00/00  Created                                              *
 *                                                                          *
 ****************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "plot.h"

double  GNa=120 ;
double	GK = 36;
double	GL=0.3;
double	R = 30;
double	VNa = 115;
double	VK = -12;
double	VL=10.5995;
double	c = 1;


//Function Declaration 

double KV( double I, double V,double m,double h,double n );
double Km( double V , double m );
double Kh( double V , double h );
double Kn( double V , double n );
double KIext ( double time , double w , double Iext0 );



int main()
{
		

	FILE *fV_vs_t_id;
	clock_t start, end;
	fV_vs_t_id = fopen("V_vs_t.txt", "w");
		
/*	FILE *out1;
	out1 = fopen("m.txt", "w");
		
	FILE *out2;
	out2 = fopen("h.txt", "w");

	FILE *out3;
	out3 = fopen("n.txt", "w");
*/		
	FILE *fI_vs_t_id;
	fI_vs_t_id = fopen("I_vs_t.txt", "w");
		
	double V,m,h,n,Iext;
  	double runTime;
	double dt,T,w,freq;
	double  time;
	double GNa,GK,GL,VNa,VK,VL,c,R;
	double alpham,betam,alphah,betah,alphan,betan ;
	double k1,m1,n1,h1,Iext1;
	double k2,m2,n2,h2,Iext2;
	double k3,m3,n3,h3,Iext3;
	double k4,m4,n4,h4;
	long  i,j, jf, jI0;
	long  n_max;
	double Iext0;
	char *forcing_param[30];

	
	start = clock();



// *********************** Input parameter values (read from file below)*****************************************************
#include</home/agastya123/PycharmProjects/ComputationalNeuroscience/HodgkinHuxleyModel/SOP/Codes/par_HH.h>




// *********************** Start the FREQUENCY loop **************************************************************		
	freq = freq_min;
	for(jf = 0 ; jf < Num_freq  ; jf++) 
	{

		// ************************** Calculate T, dt, n_max and w (all of which depend on freq - loop is over frequency) *********

		// Time period
		T=(1/freq); // in seconds
		printf("Time period (in seconds) = %3.10f \n \n",T);

		// Time step
		dt=(T/1000); // in seconds (time step is always chosen as (1/1000)th of the time period of ext current)
		printf("Time step (in seconds) = %3.10f \n \n",dt);

		// Number of time steps
		n_max = t_max*freq*1000; // n_max = t_max / dt , where dt = T/1000 (1 time period divided into 1000 time intervals), and T = 1/freq
		printf("Number of time steps  = %ld \n \n",n_max);

		// Anglular frequency
		 w = freq*2*M_PI; // angular frequency (radians per second)
		printf("Angular frequency (rad/second) = %3.10f \n \n",w);

	
		// **************************** convert parameters to milliseconds (from seconds) *********************************
	
		// Anglular frequency
		w=(w/1000); // angular frequency (radians per millisecond)
		printf("Angular frequency (rad/milli-second) = %3.10f \n \n",w);
	

		// Time step
		dt=dt*1000; // in milliseconds. Factor of 1000 converts dt to milliseconds.
		printf("Time step (in milli-second) = %3.10f \n \n",dt);


	
		// ***************************** Prepare simulation for a fresh value of forcing amplitude ***************************************
		Iext0 = Iext0_min;
		for(jI0 = 0 ; jI0 < Num_I0  ; jI0++)  // Loop over a range of current (external) amplitudes (Max number of current values = Num_I0)
		{	
			 printf("jI0= %ld Iext0= %.11f \n",jI0,Iext0);
		

			 //  ************* Initialize time, V[], m1, h1, n1 before starting RK4 *********************************
			 // Bring time back to zero
			 time=0;		 

			 // Bring back voltage and gate variables back to rest state
			 V = 0; 		 
			 alpham = (0.1*(25-V))/(exp((25 - V)/10) - 1);
			 betam = 4*exp(-V/18);
			 alphah = 0.07*exp(-V/20);
			 betah = 1/(exp((30-V)/10) + 1);
			 alphan = 0.01*(10-V)/(exp((10-V)/10) - 1);
			 betan = 0.125*exp(-V/80);
			 m = (alpham/(alpham + betam));
			 h = (alphah/(alphah + betah));
			 n = (alphan/(alphan + betan));

		
		          fprintf(fV_vs_t_id,"\n");
		          fprintf(fV_vs_t_id,"\n");
		          fprintf(fV_vs_t_id,"%s %.11f \n","#Iext0=",Iext0);
			 fprintf(fV_vs_t_id,"%s %.11f \n","#Freq=",freq);
			 fprintf(fV_vs_t_id,"%d %.11f %14.11f %14.11f %14.11f %14.11f \n",0,time*0.00,V,m,h,n);

			 fprintf(fI_vs_t_id,"\n");
		          fprintf(fI_vs_t_id,"\n");
		          fprintf(fI_vs_t_id,"%s %.11f \n","#Iext0=",Iext0);
			 fprintf(fI_vs_t_id,"%s %.11f \n","#Freq=",freq);

		         

			 // ***************************** START TIME STEPPING WITH RK4 ***************************************
			 for(i=0; i<n_max ; i++) // Loop over time for a chosen forcing parameter values
			 {
				 
				 time = (i+1)*dt;
			
				 Iext=Iext0*cos(w*time);

				 k1 = dt*KV(Iext,V,m,h,n);
	   			 m1 = dt*Km(V,m);
	    			 h1 = dt*Kh(V,h);
	    			 n1 = dt*Kn(V,n);
	    			 Iext1=dt*KIext(time,w,Iext0);

				 k2 = dt*KV(Iext+0.5*Iext1 ,V+(0.5*k1),m+(0.5*m1),h+(0.5*h1),n+(0.5*n1));
	    			 m2 = dt*Km(V+(0.5*k1),m+(0.5*m1));
	    			 h2 = dt*Kh(V+(0.5*k1),h+(0.5*h1));
	   			 n2 = dt*Kn(V+(0.5*k1),n+(0.5*n1));
	  			 Iext2=dt*KIext(time+0.5*dt,w,Iext0);

				 k3 = dt*KV(Iext+0.5*Iext2 ,V+(0.5*k2),m+(0.5*m2),h+(0.5*h2),n+(0.5*n2));
	    			 m3 = dt*Km(V+(0.5*k2),m+(0.5*m2));
	    			 h3 = dt*Kh(V+(0.5*k2),h+(0.5*h2));
	   			 n3 = dt*Kn(V+(0.5*k2),n+(0.5*n2));
	    			 Iext3=dt*KIext(time+0.5*dt,w,Iext0);

				 k4 = dt*KV(Iext+Iext3 ,V+k3,m+m3,h+h3,n+n3);
	    			 m4 = dt*Km(V+k3,m+m3);
	   		          h4 = dt*Kh(V+k3,h+h3);
	   		          n4 = dt*Kn(V+k3,n+n3);


				V = V + ((k1 + 2*k2 + 2*k3 + k4)/6);
	    			m = m + ((m1 + 2*m2 + 2*m3 + m4)/6);
	   			h = h + ((h1 + 2*h2 + 2*h3 + h4)/6);
	   			n = n + ((n1 + 2*n2 + 2*n3 +n4)/6);
					
				fprintf(fV_vs_t_id,"%ld %.11f %14.11f %14.11f %14.11f %14.11f \n",i+1,time*0.001,V,m,h,n);
				fprintf(fI_vs_t_id," %ld %.11f %.11f \n",i+1,time*0.001,Iext); 					
						 
			} // end of RK4 loop
			

			Iext0 = Iext0 + dIext0; // change the forcing amplitude value and simulate once again.
		} // end of amplitude loop		
  	
		freq = freq + dfreq; // change the forcing frequency value and simulate once again.
  	} // end of frequency loop


	fclose(fV_vs_t_id);
	//fclose(out1);
	//fclose(out2);
	//fclose(out3);
	fclose(fI_vs_t_id); 
	return 0; 
}


	double KV( double I, double V,double m,double h,double n)
	{
		double p ;
		p =  ((I - ((GNa*(m*m*m)*h*(V-VNa)) + (GK*(n*n*n*n)*(V-VK))+(GL*(V-VL))))/c);
		return (p) ;
	}

	double Km( double V , double m )
	{
		double q;
		q =  ((((0.1*(25-V))/(exp((25-V)/10) - 1))*(1 - m)) - (m*(4*exp(-V/18))));
		return (q);
	}

	double Kh( double V , double h )
	{
		double r;
		r = ((0.07*exp(-V/20))*(1-h)) - ((1/(exp((30-V)/10) + 1))*h);
		return (r);
	}

	double Kn( double V , double n )
	{
		double s;
		s = (0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)) ;
		return (s);
	}

	double KIext ( double time , double w , double Iext0 )
	{
		double t;
		t = (-w*Iext0*sin(w*time));
		return (t);

	} 

