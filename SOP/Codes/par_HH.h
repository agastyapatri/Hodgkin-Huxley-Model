
double	freq_min, freq_max, dfreq ;
long   Num_I0, Num_freq;
double Iext0_min, Iext0_max, dIext0;
double t_max;




// *********************** Input parameter values *****************************************************
	// Forcing amplitude
	Iext0_min=1.86346;  // least value of external current amplitude
        Iext0_max=1.86350;  // largest value of - do -
        dIext0 =  0.0000001;    // Difference between successive values of current amplitude
	Num_I0 = ((Iext0_max - Iext0_min)/dIext0) + 1; 
	printf("Number of forcing Amplitudes = %ld \n \n",Num_I0);

	// Forcing frequency (See HH code)
	freq_min = 50; // (in Hz - per second)
	freq_max=  50; // (in Hz - per second)
	dfreq =     2;      // (in Hz - per second)
	Num_freq = ((freq_max - freq_min)/dfreq) + 1; // number of different frequency values (of current) for which simulations will be carried
	printf("Number of different values of forcing frequencies= %ld \n \n",Num_freq);

	// Maximum time for which simulation will be run
	t_max =20;  // in seconds
	printf("Maximum time(in seconds) = %3.10f \n \n",t_max);



