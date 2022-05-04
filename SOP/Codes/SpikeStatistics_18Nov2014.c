#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(){
int Num_lengths = 61;
int lenarr[Num_lengths] = {0};
double Iext0;
double freq;
double T, dt, w;
float spikegroup;
float AveragingTime,FiringRate;
float tjunk;
double tmp1, tmp2;
float tmp3, Iext0tmp;


long kf, kI0;
long i, j, k;
long m, n, grpnum, nospike_grpnum;
long s, nos;
long ijunk;
long n_max;

long n_start_spike_count, n_end_spike_count;
long Nspikes_tot,Ngrp_tot, Ngrp_nospike_tot;
long z;
long ilocmax ;

char label1[6], label2[5];

double *V, *m1, *h1, *n1, *tspike, *tNospike;
float *lyapunov ;
//float *ISI;
long *spike_per_grp, *nospike_per_grp;
float ISI, AvgISI;


// ***** variables for repetitive sequence
long Start_grpnumL, End_grpnumL, Start_grpnumR, End_grpnumR, diff_grpnumLR;
long imatch;
long grpnumL, grpnumR;
long diff;
long Last_grp_repeating_seq, Num_grp_repeating_seq ;
long *Spikes_grp_repeating_seq;
long tempL, tempR;
long iexpand, ishift;
int temp_iter ;
long ss; 

// *********


//long spike_per_grp[100000];

FILE *fV_vs_t_id, *fjunkid, *fSpikeTimeid, *fSpike_per_grpid, *fISIid, *fRepeatSeqid , *fhistspike1, *fhistspike2, *fAvgISI ;
FILE *fNoSpike_per_grpid , *fLengthRepeatSeqid,*lenfreq;
FILE *fdebug;
//FILE *fLyap_id;
//FILE *fLyp_perdc_spike_grp, *fLyp_chaos_spike_grp, *fmultidots_per_grp ;

// *********************** Input parameter values (read from file below)***************************************************

#include</media/mohit/Seagate Backup Plus Drive/hhNeuron/Trials/t7/par_diag.h>

int count_of_seq[Num_I0] ; 
long ispike_grp, freq_of_spike_per_grp[Num_I0][Max_spike_per_grp] ;
double prob_of_spike_per_grp[Num_I0][Max_spike_per_grp];  // variables for calculating statistics of spike per group

// **************** Open Files ***************************************************************************

	fV_vs_t_id=fopen("V_vs_t.txt","r");
         fjunkid=fopen("Vjunk.txt","w");
	fSpikeTimeid = fopen("SpikeTime.txt","w");
	fSpike_per_grpid = fopen("Spike_per_grp.txt","w");
         fNoSpike_per_grpid = fopen("NoSpike_per_grp.txt","w");
	//fISIid = fopen("ISI.txt","w");
         fRepeatSeqid = fopen("RepeatSeq.txt","w");
	fLengthRepeatSeqid = fopen("LengthRepeatSeq.txt","w");
	fdebug = fopen("debug.txt","w");
	fhistspike1 = fopen("Freq_of_Spike_per_grp.txt" , "w");	
	fhistspike2 = fopen("Freq_of_Spike_per_grp2_vs_Iext0.txt" , "w");
	lenfreq = fopen("len_freq", "w");
	//fAvgISI = fopen("AvgISI.txt","w");	
	//fLyap_id = fopen("Lyapunov.txt", "r");
         //fLyp_perdc_spike_grp = fopen("Lyp_perdc_spike_grp.txt","w");
	//fLyp_chaos_spike_grp = fopen("Lyp_chaos_spike_grp.txt","w");
	//fmultidots_per_grp = fopen("Multidots_per_grp.txt","w");

// ********************** Initializing the histogram array ****************************************
        for(kI0 = 0 ; kI0 < Num_I0  ; kI0++) 
	{ 
		for(ispike_grp=1; ispike_grp <= Max_spike_per_grp; ispike_grp++)
		{
			freq_of_spike_per_grp[kI0][ispike_grp]=0;		
			prob_of_spike_per_grp[kI0][ispike_grp]=0;
		}
	}

// *********************** START LOOP OVER FORCING FREQUENCY ************************************************

	for(kf=0 ; kf<Num_freq; kf++){
		
		 printf("Reached here Kf = %ld \n", kf);
		//made changes here dete this or change if want to take more then 400 steps
		if(kf>1100)
		{
			goto label; 
		}
		// Above loop index kf is frequency index. Frequencies are read from a file


		// *********************** START LOOP OVER FORCING AMPLITUDES************************************************
		
		for(kI0 = 0 ; kI0 < Num_I0  ; kI0++) { 
			 printf("Reached here kI0 = %ld \n", kI0);
			// Above loop index kI0 is amplitude index. Amplitudes are read from a file

			// frequencies and amplitudes are read first. Frequency is needed to calculate dt and n_max needed later.
			fscanf(fV_vs_t_id,"\n");
			fscanf(fV_vs_t_id,"\n");
			fscanf(fV_vs_t_id,"%s %lf \n",label1,&Iext0);
			fscanf(fV_vs_t_id,"%s %lf \n",label2,&freq);
		        

			// print labels of frequency and Iext on every file that is created.
			fprintf(fjunkid,"#Iext0 = %.11f Freq = %.11f \n \n",Iext0, freq);
			fprintf(fSpikeTimeid,"\n \n");
			fprintf(fSpikeTimeid,"#Iext0 = %.11f Freq = %.11f \n",Iext0, freq);
			//fprintf(fISIid,"\n \n");
			//fprintf(fISIid,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			fprintf(fSpike_per_grpid,"\n \n");
			fprintf(fSpike_per_grpid,"#Iext0 = %.11f Freq = %.11f \n",Iext0, freq);
			fprintf(fRepeatSeqid,"\n \n");
			fprintf(fRepeatSeqid,"#Iext0 = %.11f Freq = %.11f \n",Iext0, freq);
			fprintf(fhistspike1,"\n \n");
			//fprintf(fhistspike1,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			//fprintf(fAvgISI,"\n \n");
			//fprintf(fAvgISI,"#Iext0 = %f Freq = %f \n",Iext0, freq);
			//fprintf(fmultidots_per_grp,"\n \n");
			//fprintf(fmultidots_per_grp,"#Iext0 = %f Freq = %f \n",Iext0, freq);

			


			//********************************  CALCULATE T, dt, n_max ********************************************
			// ************************dt needed for calculating n_start_map_constrn *******************************
			// **************n_max = number of voltage values generated by RK4 (depends on the forcing frequency ******
			
			// Time period
			T=(1/freq); // in seconds
			printf("Time period (in seconds) = %f \n \n",T);

			// Time step
			dt=(T/1000); // in seconds (time step is always chosen as (1/1000)th of the time period of ext current)
			printf("Time step (in seconds) = %f \n \n",dt);

			// Anglular frequency
			w = freq*2*M_PI; // angular frequency (radians per second)
			printf("Angular frequency (rad/second) = %f \n \n",w);

			n_max = t_max*freq*1000;  //n_max = t_max/dt; 
			printf("Number of time steps = %ld \n \n",n_max);

			// ****************** size allocation to arrays (depends on n_max. n_max depends on frequency)************* 

			V = (double *)malloc(n_max*sizeof(double));    // Voltage
			m1 = (double *)malloc(n_max*sizeof(double));  // m-gate
			h1 = (double *)malloc(n_max*sizeof(double));   // h-gate
		         n1 = (double *)malloc(n_max*sizeof(double));   // n-gate
			tspike = (double *)malloc((n_max/n_one_cycle)*sizeof(double));
			tNospike = (double *)malloc((n_max/n_one_cycle)*sizeof(double));
			//ISI = (float *)malloc((n_max/n_one_cycle)*sizeof(float));
			spike_per_grp = (long *)malloc((n_max/n_one_cycle)*sizeof(long));
			nospike_per_grp = (long *)malloc((n_max/n_one_cycle)*sizeof(long));
			Spikes_grp_repeating_seq = (long *)malloc(n_max*sizeof(long));
			lyapunov = (float *)malloc((Num_I0+2)*sizeof(float));
		
	
			//*******************************************************************************************************
			// *					 READ VOLTAGE TIME SERES FROM A FILE  			
			// ******************************************************************************************************
			printf("Start Reading voltage time series\n");
			for(i=0;i<=n_max;i++) 
       			{
			fscanf(fV_vs_t_id,"%ld %f %lf %lf %lf %lf \n",&ijunk,&tjunk,&V[i],&m1[i],&h1[i],&n1[i]);  
			}
			//fscanf(fLyap_id,"%f %f \n",&Iext0tmp,&lyapunov[kI0]);  
			

			// ***********************************************************************************************************
			//* CALCULATE - SPIKE TIMINGS & NUMBER OF SPIKES IN EACH GROUP  					
			//************************************************************************************************************
			

			n_start_spike_count = T_start_spike_count/dt + 1; 
			n_end_spike_count = T_end_spike_count/dt + 1 ; 
                        //Once started poincare map construction goes on till the end of the time series		

			i=n_start_spike_count;
			
			s=0; //n is the spike counter
			nos = 0; // m is the no spike counter


			printf("Number of cycles of input current = %ld \n",n_max/n_one_cycle);

			for(grpnum=0;grpnum<=n_max/n_one_cycle;grpnum++)  //Initialize arrange spike_count[] = 0. This array stores the number of spike in each group.
			{
			 	spike_per_grp[grpnum]=0;  //grp_spikecount[grpnum] is the number of spikes in the group grpnum
				nospike_per_grp[grpnum]=0;
			}

			for(ss=0;ss<=(n_max/n_one_cycle);ss++)
			{
				tspike[ss]=0;		// tspike[n] is the time at which nth spike (maxima in V) occurs
				tNospike[ss]=0;
			}

			grpnum = 0;
			nospike_grpnum = 0;
			ilocmax = 0;
			printf("Entering spikes-per-group counting loop \n");
			for(i=n_start_spike_count;i<=n_end_spike_count-1;i++)	  //Iterating on time steps
			{
				//printf("i= %ld \n",i);
				tmp1 = V[i]-V[i+1];	
				tmp2 = V[i]-V[i-1];
				if(tmp1 >= 0 && tmp2 >= 0)
				{
					//printf("i= %ld \n",i);
					if(V[i]>70)
					{
						s++ ; 
						//printf("s= %ld \n",s);                 
						tspike[s]=i*dt ; // tspike[n] is the time at which nth spike (maxima in V) occurs
						//printf("s= %ld \n",s);   
						if (ilocmax == -1)
                                                   		nospike_grpnum++ ;
						spike_per_grp[grpnum]++ ;   
						ilocmax = +1 ;
						//printf("s= %ld \n",s);
						//printf("i= %ld grpnum = %ld spike_per_grp = %ld\n",i, grpnum,spike_per_grp[grpnum]);
						fprintf (fSpikeTimeid,"%ld %f %ld %ld\n",s,tspike[s],grpnum, ilocmax);
						
					}
					else if(V[i]<70)
					{
						nos++ ;
						tNospike[nos] = i*dt ; // tNospike
						if (ilocmax == +1)
							grpnum++ ; // once a missing spike is detected increment grpnum by one.
						nospike_per_grp[nospike_grpnum]++ ;
						ilocmax = -1;
						//printf("grpnum = %ld \n",grpnum);

						fprintf (fSpikeTimeid,"%ld %f %ld %ld\n",nos,tNospike[nos],grpnum, ilocmax);
					
					}

				}
			}

			printf("Number of spike groups = %ld \n",grpnum);
			Ngrp_tot = grpnum+1; // total number of spike groups
			Ngrp_nospike_tot = nospike_grpnum ; // Total number of No-spike groups
			
		
			Nspikes_tot = s ; // number of spikes that were generated
			printf("Number of spikes = %ld \n",Nspikes_tot);
			printf("Number of spike groups = %ld \n",Ngrp_tot);

			//Print - number of spikes in each group

			printf("Writing in file - number of spikes per group \n");
			
			for(grpnum=1;grpnum<=Ngrp_tot;grpnum++) 
			// grpnum starts from '1' instead of '0'. no spike may appear before a spike. Then grpnum
			// would get incremented from 0 to 1 even befor a spike has appeared.
			{
			   fprintf (fSpike_per_grpid,"%ld %ld \n",grpnum,spike_per_grp[grpnum]);
			}

			printf("Writing in file - number of No spikes per group \n");

			// Frequency of multiple dots is print is printed later

			/*for(nospike_grpnum=1;nospike_grpnum<=Ngrp_nospike_tot;nospike_grpnum++) 
			{
			   fprintf (fNoSpike_per_grpid,"%ld %ld \n",nospike_grpnum,nospike_per_grp[nospike_grpnum]);
			   if(nospike_per_grp[nospike_grpnum] > 1)
			   fprintf(fmultidots_per_grp,"%f %f %ld \n",Iext0,lyapunov[kI0],nospike_per_grp[nospike_grpnum]);
			}* */

			

			// *****************************************************************************************************************
			// *		       Calculate & Print -  Interspike interval (ISI) 						   *
			// *****************************************************************************************************************
/*                        
			printf("Calculate ISI \n");
			ISI = 0;
			AvgISI = 0;
			
			for (j=1;j<=Nspikes_tot-1;j++){ // should it be j<=Nspikes_tot-1 or j<Nspikes_tot-1
			  	//ISI[j]=tspike[j+1]-tspike[j];
	                        ISI = tspike[j+1]-tspike[j];
 				AvgISI = AvgISI + ISI ;
			  	//fprintf(fISIid,"%ld %f \n",j,ISI[j]);
				fprintf(fISIid,"%ld %f \n",j,ISI);	
			}
			AvgISI = AvgISI/(Nspikes_tot-1) ;
			fprintf(fAvgISI,"%f %f \n",Iext0,AvgISI);

			printf("Finished ISI \n");
*/			// ***************************************************************************************************
			// **				Histogram of # of spikes in a group 
			// ****************************************************************************************************

			printf("Start calculation of Histogram \n");

	        		for(ispike_grp=1;ispike_grp<=Max_spike_per_grp;ispike_grp++) 
			// for finding frequency of number of spikes in a grp equal to nspike_grp
           		{
				for(grpnum=2;grpnum<=Ngrp_tot-2;grpnum++) //scan all the groups
				{
			   		if(spike_per_grp[grpnum]== ispike_grp) // if the number of spikes in that group = ispike_grp
					{

					freq_of_spike_per_grp[kI0][ispike_grp] = freq_of_spike_per_grp[kI0][ispike_grp] + 1 ; 
					// MISLEADING NAME: freq_of_spike_per_grp[kI0][ispike_grp] denote the number of groups with a particular number of spikes

					}
				}

				tmp3 = freq_of_spike_per_grp[kI0][ispike_grp] ;
				prob_of_spike_per_grp[kI0][ispike_grp] = tmp3/(Ngrp_tot-3) ;
			//printf("Probability = %ld  %f \n",Ngrp_tot-2,prob_of_spike_per_grp[kI0][ispike_grp]);	
			fprintf(fhistspike1,"%ld %ld %.11f\n",ispike_grp,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
			}	

			printf("Complete calculation and printing of Histogram \n");

			printf("periodm= %ld\n",periodm);

		

			// ****************************************************************************************************************
			// ********************* Average firing rate *****************************************************************
			// ****************************************************************************************************************
			//AveragingTime=(n_end_spike_count - n_start_spike_count)/dt; //.....is this formula correct or does -1 have to be subtracted from top?
			//FiringRate=Nspikes_tot/AveragingTime;
			//fprintf(............)

		// ********************************************************************************************************************			
			// ********************************************************************************************************************	
			// *   					IDENTIFY REPETITIVE SPIKE SEQUENCE					       *
		// *********************************************************************************************************************
			// *********************************************************************************************************************





			Start_grpnumL = 1; // The 1st group in all sequences that will be tested for repetivity
					   // (this value will not change hereafter. 1st group in sequence is always group no. 0)
			End_grpnumL = 1;   // the last group in the current sequence that will be tested for repetivity.
					   // if the current sequence is not repetitive, sequence size will be increased by 1.
					   // End_grpnumL will thus be increased by 1.
			iexpand = 0;
			ishift = 0;
			do // this loop expands the size of sequence till a repetitive sequence is found)
			{	
				 iexpand++ ;
				 fprintf(fdebug,"iexpand= %ld ishift = %ld \n", iexpand, ishift);
				// sequences labeled by L (left in time) will be matched against sequences labeled by R (right in time).
				// Matched 'against' sequence will start at group Start_grpnumR and end at group End_grpnumR.
				// The 1st grp in matched 'against' sequence will start just after the last group of the 'being' matched sequence.

				Start_grpnumR = End_grpnumL+1 ; // the group number on right (in time) where the matching starts right after
								// 'i' is incremented by 1 (that is after R sequence is shifted)
						

				diff_grpnumLR = Start_grpnumR - Start_grpnumL; // difference in group nos where left sequence and right sequence begin.

				imatch = +1;  // if L sequences matches with R sequence , imatch = +1 (imatch also initialized to +1)

				i=1; 
	
        				do   // this loop shifts the R sequence against which the L sequence will be matched.
       				{
					ishift++;
                                       	fprintf(fdebug,"ishift= %ld iexpand= %ld \n", ishift, iexpand);
					// grpnumL is the group label in left(L) sequence. This loop is over groups in a chosen L sequence 
					for(grpnumL=Start_grpnumL; grpnumL<=End_grpnumL; grpnumL++) // loop over all groups in L sequence
					{
	
						grpnumR =  grpnumL+i* diff_grpnumLR ;// (grpnumL will be matched 'against' grpnumR)
									
						//diff = spike_per_grp[grpnumL] - spike_per_grp[grpnumL+i*diff_grpnumLR];
							   // diff in number of spikes in grps being matched in left & right sequences
						tempL = spike_per_grp[grpnumL];
						tempR = spike_per_grp[grpnumL+i*diff_grpnumLR];

					 	fprintf(fdebug,"grpnumL= %ld grpnumR= %ld imatch=%ld",grpnumL, grpnumR,imatch);						
						fprintf(fdebug,"tempL= %ld tempR = %ld \n", tempL, tempR);
						
					 	if(tempL != tempR) 											
						//if (fabs(diff) != 0) //(if diff is non-zero sequence L does not repeat)
						{
							imatch = -1; // repetitive sequence not found (imatch = -1)
						}

					}  // loop continues till all groups in sequence L have been compared against the R sequence
					   // (even if imatch == -1 already. not efficient. can be rectified later)

					i++;  // 'i' locates the R sequence. Each R sequence has same length as L sequence.
					      //  But the R sequence can be 1 sequence away, 2 sequences away, or i sequence away from L sequence.

				 }while (imatch == +1 && End_grpnumL+i*diff_grpnumLR <= Ngrp_tot-2);

        				End_grpnumL++ ; // change the number of groups in L sequence by 1 (since the shorter sequence was not repetitive)
				
				fprintf(fdebug,"Now L sequence will be expanded , imatch= %ld \n",imatch);

			}while(imatch == -1 && End_grpnumL <= (Ngrp_tot-2)/2); 
								//if a sequence matches all  other same size sequences imatch=+1 and we
								 //have found a repetitive sequence. Job is done.On the other hand if a 
								//repetitive seq is not found imatch = false and we must go on.
			fprintf(fdebug,"Repetitive sequence found, imatch= %ld \n",imatch);
    
			Last_grp_repeating_seq = End_grpnumL -1 ; // This is the length of the repeating sequence
			
			Num_grp_repeating_seq = Last_grp_repeating_seq;

			if (imatch == +1) // Periodic sequence found
			{
				for(k = 1;k<= Last_grp_repeating_seq;k++) 
				{
			  		Spikes_grp_repeating_seq[k] = spike_per_grp[Start_grpnumL+k-1];
				 	fprintf(fRepeatSeqid,"%ld \n",Spikes_grp_repeating_seq[k]);
				}
				fprintf(fLengthRepeatSeqid,"%.11f %ld \n",Iext0,Last_grp_repeating_seq);
				lenarr[Last_grp_repeating_seq]++;
				printf("Repeating sequence is found for I0=%.11f \n--------------------------------------------------------\n\n", Iext0);
			}
			else if (imatch == -1) // No periodic sequence found
			{
				fprintf(fRepeatSeqid,"%d \n", 0);
				fprintf(fLengthRepeatSeqid,"%.11f %d \n", Iext0,0);
				lenarr[0]++;
				printf("No repeating sequence is found for I0=%.11f \n--------------------------------------------------------\n\n", Iext0);
			}

//************************************** Frequency of multiple dots in a group***************************************************
/*			for(nospike_grpnum=1;nospike_grpnum<=Ngrp_nospike_tot;nospike_grpnum++) 
			{
			  // fprintf (fNoSpike_per_grpid,"%ld %ld \n",nospike_grpnum,nospike_per_grp[nospike_grpnum]);
			   //if(nospike_per_grp[nospike_grpnum] > 1)
			   fprintf(fNoSpike_per_grpid,"%f %f %d %ld %ld \n",Iext0,lyapunov[kI0],imatch,nospike_grpnum,nospike_per_grp[nospike_grpnum]);
			   //fprintf(fMultidots_per_grp,"%f %f %d  %ld %ld \n",Iext0,lyapunov[kI0],imatch,nospike_grpnum,nospike_per_grp[nospike_grpnum]);
			}
*/
// ********************************* EXPLORE : GROUPS WITH HOW MANY SPIKES CAN OCCUR ********************************************
// ********************************* IN A PERIODIC SEQUENCE LYING BETWEEN period-m and period-(m+1) ****************************
			// Do periodic sequences between period-3 ({2.}) and period-4 ({3.}) region only display 2 spikes or 3 spikes in a
			// group. OR . can periodic sequneces also display 1 or 4,...spikes in a group 

			// Below we will explore the above question in region between period-(periodm) and period-(periodm+1).
			// The value of 'periodm' can be chosen by the user in file par_diag.h 
		
			// Below, we run through all possible values of 'number of spikes in a group'. If this number is non-zero
			// then find out if the Lyapunov exp is +ve or -ve. 

			// We expect the following :
			// For periodic sequence (Lyapunov exp -ve) in region between period-m ({m-1.})and period-(m+1) [{m.}]
			// f(i) = 0 for i=1,2,...,m-2, and
			// f(i) = 0 for i = m+1,.....
			// f(i) != 0 for i = m-1 and m.
			// However for chaotic sequence (Lyapunov exp +ve) f(i) could be non-zero even if i=1,2,..,m-2 and for i=m+1,...

			// Hence, for the chosen region (period-m and period-m+1) we wish to 
			// record for each Iext0 the Lyapunov exponent , i f(i) for i=1,...m-2 (Ask - f(i) = 0 for Lyapunov exp -ve ?)
			//                           Lyapunov exponent , i,f(i) for i = m+1,....(Ask - Is f(i) = 0 for Lyapunov exp -ve?)
			//        
	
			// find frequency f(n) corresponding to Iext0 for which Lyapunov -ve 
/*			for (ispike_grp = 1; ispike_grp <= periodm-2; ispike_grp++)//loop goes from 1 spike to (periodm-2) spikes in a grp
			{	
			if(freq_of_spike_per_grp[kI0][ispike_grp]!=0 && lyapunov[kI0] < 0) //note down only if f(n) is not zero.
			fprintf(fLyp_perdc_spike_grp,"%f%f%d%ld%ld%f\n",Iext0,lyapunov[kI0],imatch,ispike_grp,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
			}
			for (ispike_grp = periodm+1; ispike_grp <= Max_spike_per_grp; ispike_grp++)//go from (periodm+1) spike per grp onwards
			{	
			if(freq_of_spike_per_grp[kI0][ispike_grp]!=0 && lyapunov[kI0] < 0)
			fprintf(fLyp_perdc_spike_grp,"%f%f%d%ld%ld%f\n",Iext0,lyapunov[kI0],imatch,ispike_grp,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
			}
			// find f(n) for Lyapunov +ve
			for (ispike_grp = 1; ispike_grp <=  periodm-2; ispike_grp++)//loop goes from 1 spike to (periodm-2) spikes in a grp
			{	
			if(freq_of_spike_per_grp[kI0][ispike_grp]!=0 && lyapunov[kI0] > 0)
			fprintf(fLyp_chaos_spike_grp,"%f%f%d%ld%ld%f\n",Iext0,lyapunov[kI0],ispike_grp,imatch,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
			}
			for (ispike_grp = periodm+1; ispike_grp <= Max_spike_per_grp; ispike_grp++)//go from (periodm+1) spike per grp onwards
			{	
			if(freq_of_spike_per_grp[kI0][ispike_grp]!=0 && lyapunov[kI0] > 0)
			fprintf(fLyp_chaos_spike_grp,"%f%f%d%ld%ld%f\n",Iext0,lyapunov[kI0],imatch,ispike_grp,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
			}
			
			
*/

// **********************************************************************************************************************************
		} // amp loop
	} // freq loop

 // ****************** Histogram of # of spikes in a group. It was formed earlier in the code already and for each amplitude ****************************
// *******************************************the numer of spikes per group vs  frequency was also saved in a file ***************************************
// ************************************Here we will plot frequency vs amplitude for each value of 'number of spikes per group' ***************************


label:
for(ispike_grp=1;ispike_grp<=Max_spike_per_grp;ispike_grp++) // for finding frequency of number of spikes in a grp equal to nspike_grp
{
	fprintf(fhistspike2,"\n \n");
	fprintf(fhistspike2,"#Number of spike per group = %ld \n",ispike_grp);
	for(kf=0 ; kf<Num_freq; kf++)
	{
    		freq = kf*dfreq + freq_min;
		for(kI0 = 0 ; kI0 < Num_I0  ; kI0++) 
		{ 
                Iext0 = kI0*dIext0 + Iext0_min;
		fprintf(fhistspike2,"%.11f %ld %f\n",Iext0,freq_of_spike_per_grp[kI0][ispike_grp],prob_of_spike_per_grp[kI0][ispike_grp]);
		}
	}
}
for(int i = 0; i< Num_lengths; i++)
{ 
	fprintf(lenfreq,"%d \n",lenarr[i]);
}	
fclose(lenfreq);
fclose(fLengthRepeatSeqid);
fclose(fhistspike1);
fclose(fhistspike2);
} // main loop











