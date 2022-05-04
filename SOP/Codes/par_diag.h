
float	freq_min, freq_max, dfreq ;
long   Num_I0, Num_freq;
//here
double Iext0_min, Iext0_max, dIext0;
double t_max;
long   n_one_cycle;
float T_start_map_constrn;
float T_start_spike_count;
float T_end_spike_count;
long Max_spike_per_grp;
long periodm;


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
	dfreq =     0;      // (in Hz - per second)
	Num_freq = 1; // number of different frequency values (of current) for which simulations will be carried
	printf("Number of different values of forcing frequencies= %ld \n \n",Num_freq);

	// Maximum time for which simulation will be run
	t_max =20;  // in seconds
	printf("Maximum time(in seconds) = %3.10f \n \n",t_max);

	// Number of time steps in one cycle of sinusoidal input current (later dt is chosen as dt = T/1000). If a change is made,1000 should be changed both places.
	n_one_cycle=1000;
	printf("Number of time steps in one cycle  = %ld \n \n ",n_one_cycle);

	// Adjust dimensions, if required
	long j_mode_lck[100][100];

	// Start Poincare map construction at the following time
	T_start_map_constrn = 15; // time at which poincare map construction starts (in seconds)
	printf("Poincare map constructruction starts at time (in seconds) = %f \n \n ",T_start_map_constrn);

	// Start Poincare map construction & also spike counting at the following time
	T_start_spike_count = 14; // time at which spike counting starts (in seconds)
	printf("Spike counting starts at time (in seconds) = %f \n \n ",T_start_spike_count);

	// End Poincare map construction & also spike counting at the following time
	T_end_spike_count = 19; // time at which spike counting ends (in seconds)
	printf("Spike counting ends at time (in seconds) = %f \n \n ",T_end_spike_count);

	// For the statistics of spike numbers per group we set an upper limit on its value. The statistics will be carried out only till this number
        Max_spike_per_grp = 30;

	// We will here look at statistics in the range from 1/periodm to 1/(periodm+1). We wish to find out whether periodm-2 and periodm+1 
        // spikes per group occur in this range.
	periodm = 3 ;
