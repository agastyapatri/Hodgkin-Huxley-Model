

23/09/21 
\
There seems to be a tiny phase difference between the membrane spiking and the 
external current. To elucidate, the spike initiation is not 
taking place at exactly the same time during the Current cycle. 

1. How to quantify this?
2. Is a dot simply a result of phase difference adding up?


29/09/21
\
Continue investigation on the phase difference between the current and the membrane 
voltage.

1. Check if the spike is initiating at the same points on the current cycle 
2. Explore the time constants of the h gate
3. Explore the h gate of the reduced HH model via phase plane analysis

4/10/2021
\
Smaller values of Tau result in faster dynamics of the variable.
\
Correlate the different values of time constant to different types of dynamics 


-------------------------------------------------------------
* **Action Potential**, _Neuronal Dynamics_
\
When V = V_rest, all inward and outward currents balance each other so that 
the net current is 0.

    A small pulse of current applied via I(t) produces a small positive
perturbation of the membrane potential (depolarization), which results in a small net current that 
drives V back to the resting (Repolarization). 

    An intermediate size pulse current produces a perturbation 
that is amplified far more due to the conductances being dependant on V. 
\
This nonlinear amplification causes V to deviate significantly from V_rest, and is called 
_**Action Potential**_


* Strong depolarization increases activation variables _m, n_ and 
    decreases the inactivation variable _h_. 
Since Tau_m is relatively small, _m_ is relatively fast (_i.e the activation of Na+ conductance is fast_)
Fast activation of the Na+ conductance drives V toward E_Na, resulting in further depolarization and further activation of g_Na.
    This positive feedback loop results in an upstroke of V. 


* While V moves towards E_Na, _h, n_ catch up. _h_ tends to 0 and _n_ tends to 1, causing a slow activation of of the outward K+ current. 
* The K+ and the Leak currents repolarize the membrane potential towards V_rest
* Tau_m and Tau_h are relatively large, so the recovery of variables _n, h_ is slow. 
* The outward K+ current continues to be activated even after the action potential downstroke, resulting in V to go below V_rest towards E_k
    This is called **_After Hyperpolarization_**


----------------------------------------------------------------------------------
3/11/21

Recovering mix mode oscillations by keeping _m_ instantaneous and changing _Tau_m_ and _Tau_n_

-----------------------------------------------------------------------------------
24/1/22

Keeping m instantaneous, Tau_h can be tuned in such a way that there is only one missing spike.
~ for (factor_m, factor_h, factor_n) = (1/10e10, 4, 1) there is only one missing spike.

 
-----------------------------------------------------------------------------------
10/2/22

For the same values of (_factor_m, factor_h, factor_n_), check how the different values of the Amplitude effect the membrane potential.

For different current amplitudes, find what values of (_factor_m, factor_h, factor_n_) give the same qualitative behaviour.
    1. I_0 = 1.85 at (1/10e10, 4, 1) and I_0 = 1.55 at (1,1,1) show the same pattern 
    2. 

----------------------------------------------------------------------------------
15/2/22
REDUCING THE HH MODEL TO A 2 DIMENSIONS

The section 3.6 (pg 66) of "MATHEMATICAL FOUNDATIONS OF NEUROSCIENCE" describes two ways in which dynamical systems methods can be 
used to formally reduce the four dimensional HH model to a two dimensional system of equations. 

* METHOD 1 of REDUCTION: Rinzel developed a simple method based on two observations: 
    1. Tau_m is much smaller than Tau_h and Tau_n. This implies that m(t) is close to m0(V(t)). If m is replaced by m0 in the 
        voltage equation, then the HH model is reduced by one equation. 
    2. The second observation is that (n, h) lies nearly along the line n = b - r*h, where b and r are constants.
       (b, r) depend on the current, but that can be ignored. IF we replace n by b-rh in the voltage equation, we obtain 
        a 2 DIMENSIONAL MODEL.
    
        dVdt = (1/C)*(I_ext - (g_Na*m0**3*h*(V-V_Na) + g_k*(b-r*h)**4*(V-V_k) + g_l*(V -V_l)) )
        dhdt = alpha_h * (1-h) - beta_h * h  
-------------------------------------------------------------------------------------
    
21/2/22

PLOTTING A BIFURCATION DIAGRAM

For each value of the amplitude I_0, there is an array of membrane voltages of the same size as the timescale. 
How then can a bifurcation plot be done?

22/2/22
Plot the value of the Voltage whenever current is at the maximum value. 

describe what is understood in the paper and ask what is not understood.

---------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
29/3/22

Changed Params: Investigating the region between (n)-dot and (n+1)-dot. The Bifurcation Region. 
1. BETWEEN 2 DOT and 3 DOT
    ~ I_0 = 2.25: 2-dot.
    ~ I_0 = 2.275: 3-dot followed by 2-dot. 
    ~ I_0 = 2.29: 2 dot.
    ~ I_0 = 2.35 : 2-dot.
------------------------------------------------------------------------------------------------------------------------
11/4/22
I need to look into why the bifurcation diagrams in the smaller regions are yielding the wrong plots. 
Am i collecting the wrong values of the membrane potential for the step size?

What needs to be done: 
    re-write the loop to append the currents and the Voltages 
------------------------------------------------------------------------------------------------------------------------
12/04/22

The code seems to be fixed for the most part. 
What i need to do now is:
    1. Check if i can recover results from the Colab Code. 
    2. Go through the Bifurcation regions and see what is happening. 

------------------------------------------------------------------------------------------------------------------------
<<<<<<< HEAD
13/04/2022

Trying new things: 
    1. Replace m with steady state m to make sure the error is not from there 
    2. Increase the resolution of the time and re-run the simulations. 

------------------------------------------------------------------------------------------------------------------------
=======
19/04/2022
turns out that replacing m with m_infinity is done by introducing m_infinity as a function of V, instead of the steady 
value of m. Im a dunce. 
------------------------------------------------------------------------------------------------------------------------
24/04/2022

The time period being 0.02, choose the time step such that the sampling occurs once every peak. 
It s found that sampling wasnt optimal. for the time array t = np.linspace(0, 5000, x), x being the number of time steps,
sampling used to happen when t = 25, 125, 225.. 

But, with a change in the sampling frequency (every 400 time points), it is found that the voltage will be sampled once 
every Time Period.


I need to find a value of time steps, and the sampling frequency that will be suitable for the simulation.

(time_step, sampling_frequency) = (100000, 400)
------------------------------------------------------------------------------------------------------------------------
25/04/2022

Bifurcation: Sampling once every current peak is returning weird results; every current value is returning only two 
potential values, ~ 112 and ~ -4.2. The first one i understand, but the second is bizzare. 
I need to go through the voltage array manually and see what's happening at the relevant indices.

UPDATE: at the relevant indices (every [500 + n(400)]th index), the current reaches its amplitude and there are two 
values of the potential, no matter the current amplitude. 
------------------------------------------------------------------------------------------------------------------------
27/4/2022

Restructuring of the code: Why does the reduced model in colab give only 1-dot patterns for any value of the amplitude.

UPDATE: Good News(?): Code written on pycharm yields promising results about the reduced model. For I_0 - 2.75, there 
are 6 different values, which is what is expected according to the bifurcation diagram. 

Seems like it works. Now the only thing to do is to gather data files for a few different values of the amplitude and 
then finally go ahead with the bifurcation diagram.

Its Working!!

Note: To make the bifurcation function work, All the definitions of I_0 need to be changed
------------------------------------------------------------------------------------------------------------------------
4/05/2022 
Many Things have happened, including the deletion of the Git Repository associated with this project. I need to be more 
careful about the tools that I use to conduct my work. 

Figure out why i cannot push the changes to the remote repository.

The Problem seems to persist, where for whatever reason the value of I_0 defined in the beginning is the only one 
which is chosen for the bifurcation. How the fuck do I fix this?

UPDATE: Turns out that I_0 is declared globally (enabling the definition of other functions), and so declaring a 
"new" I_0 within other functions will not update the value of the variable. 
_A fix seems to the Bifurcation Problem seems to be to run the logic of the bifurcation outside the bounds of a function 
/ "Globally"_


FINAL BIFURCATION UPDATE: The Code finally works in its entirety. I'm not quite sure why i cannot write a function 
to do the bifurcation analysis, but i've been struggling with this for so long that I'm just relieved that I have 
a working solution.


Things are looking up. 
------------------------------------------------------------------------------------------------------------------------
Print the n and h data as well 





