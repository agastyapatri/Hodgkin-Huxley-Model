"""
    Simulating the Hodgkin-Huxley model of Action Potential under a sinusoidal input current.
    Written by Agastya Patri.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, ode




"""---------------------------------------------------------------------------------------------------------------------
    1. Defining the constants used in the model 
    ~ Cm: the membrane capacitance of the neuron in uF/cm^2
    ~ V_Na, V_K, V_L: the reversal potentials of the Na, K, L channels in mV
    ~ g_Na, g_K, g_L: the maximal conductances of the Na, K, L channels in mS/cm^2
    ~ t: time over which the integration is taking place
    ~ init_cond: the initial values of the parameters  
---------------------------------------------------------------------------------------------------------------------"""

# Capacitance
C_m = 1.0

# Nernst Potential
V_Na = 115.0
V_K = -12.0
V_L = 10.5995

# Maximal Conductances
g_Na = 120.0
g_K = 36.0
g_L = 0.3

# Temperature
T = 6.3




"""---------------------------------------------------------------------------------------------------------------------
    2. Defining Auxiliary functions for the quantities used in the model
---------------------------------------------------------------------------------------------------------------------"""



def Phi(T):
    """
    Temperature dependence function. Based on the Q10 temperature sensitivity and the Arrhenius Law.
    :param T:
    :return:
    """
    phi = 3**((T-6.3)/10)
    return phi

# Gating Currents-------------------------------------------------------------------------------------------------------
def I_Na(V, m, h):
    """
    Sodium Gating Current
    """
    current = g_Na * m ** 3 * h * (V - V_Na)
    return current

def I_K(V, n):
    """
    Potassium Gating Current
    """
    current = g_K * n**4 * (V - V_K)
    return current

def I_L(V):
    """
    Leakage Current. Mainly Chloride.
    """
    current = g_L * (V - V_L)
    return current

# Channel Gating Kinetics-----------------------------------------------------------------------------------------------
def alpha_m(V):
    a_m = Phi(T)*0.1*(25-V)/(np.exp((25-V) / 10.0) - 1)
    return a_m

def beta_m(V):
    b_m = Phi(T)*4.0*np.exp(- V / 18.0)
    return b_m

def alpha_h(V):
    a_h = Phi(T)*0.07*np.exp(-V / 20.0)
    return a_h

def beta_h(V):
    b_h = Phi(T)*1.0/(np.exp((30.0-V)/10.0) + 1.0 )
    return b_h

def alpha_n(V):
    a_n = Phi(T)*0.01*(10.0 - V)/( np.exp( (10.0 - V) / 10.0) - 1 )
    return a_n

def beta_n(V):
    b_n = Phi(T)*0.125*np.exp(-(V) / 80.0)
    return b_n


# External Current Stimulus---------------------------------------------------------------------------------------------
def I_ext(time, amplitude):
    """
    This function returns current in mA/cm^2. Not sure how exactly the dimensions match in the equations, but the \
    results corroborate the expectations.
    """

    current_sinusoidal = amplitude*np.sin(2*np.pi*0.05*time)
    current_constant = 10
    current_pulse = 10*(time>100) - 10*(time>200) + 35*(time>300) - 35*(time>400)
    return current_sinusoidal


# Time Constants--------------------------------------------------------------------------------------------------------
def Tau_m(V):
    # tau = (x0 - x) / (alpha_m(V)*(1-x) - beta_m(V)*(x) )
    # tau = 1 / (alpha_m(V) + beta_m(V))
    tau = Phi(T) / (alpha_m(V) + beta_m(V))
    return tau


def Tau_h(V):
    # tau = (x0 - x) / (alpha_h(V)*(1-x) - beta_h(V)*(x) )
    tau = 1 / (alpha_h(V) + beta_h(V))
    return tau

def Tau_n(V):
    # tau = (x0 - x) / (alpha_n(V)*(1-x) - beta_n(V)*(x) )
    tau = 1 / (alpha_n(V) + beta_n(V))
    return tau




"""---------------------------------------------------------------------------------------------------------------------
    3. Defining the Hodgkin-Huxley Model
---------------------------------------------------------------------------------------------------------------------"""

def steady_state(v0):
    """
    Function to return the steady state values of the Gating Variables
    """
    m0 = float( alpha_m(v0) / (alpha_m(v0) + beta_m(v0)) )
    h0 = float( alpha_h(v0) / (alpha_h(v0) + beta_h(v0)) )
    n0 = float( alpha_n(v0) / (alpha_n(v0) + beta_n(v0)) )
    return [v0, m0, h0, n0]

init_cond = steady_state(-25)

# Parameters
t = np.linspace(0, 5000, 50000)

I_0 = 2.18

factor_m = 1/10e10
factor_h = 4
factor_n = 1



def HHModel(X, t, I_0):
    """
    Re-casting the HH model so that time constant manipulation is easier.
    :param X: The array of variables
    :param t: time
    :return: The 4 Equations of the Hodgkin Huxley Model
    """

    V = X[0]
    m = X[1]
    h = X[2]
    n = X[3]



    dVdt = ( (I_ext(t, I_0) ) - I_Na(V, m, h) - I_K(V, n) - I_L(V)) / (C_m)

    # dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
    dmdt = ( alpha_m(V) * (1.0 - m) - beta_m(V) * m ) /( (Tau_m(V)*factor_m)*(alpha_m(V)  + beta_m(V))  )

    # dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h
    dhdt = (alpha_h(V) * (1.0 - h) - beta_h(V) * h) / ( (Tau_h(V)*factor_h) * (alpha_h(V) + beta_h(V)))

    # dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n
    dndt = (alpha_n(V) * (1.0 - n) - beta_n(V) * n) / (Tau_n(V)*factor_n * (alpha_n(V) + beta_n(V)))

    return [dVdt, dmdt, dhdt, dndt]



params = g_Na,g_K,g_L,V_Na,V_K,V_L





def RunModel():

    """
    A function to integrate the equations of the HH model with respect to time and the initial conditions specified.
    :return: The integrated equations, as well as the gating currents.
    """

    Y = odeint(HHModel, init_cond, t, args=(I_0, ))

    V = Y[:,0]
    m = Y[:,1]
    h = Y[:,2]
    n = Y[:,3]

    INa = I_Na(V, m, h)
    IK = I_K(V, n)
    IL = I_L(V)

    return V, m, h, n, INa, IK, IL

V, m, h, n, INa, IK, IL = RunModel()



"""---------------------------------------------------------------------------------------------------------------------
   4. Plotting the Results
---------------------------------------------------------------------------------------------------------------------"""


def PlotResults():

    """
    Function to plot the integrated functions. Uncomment the relevant lines of code to plot what is needed
    :return: Responses of the variables with time, voltage.
    """

    plt.figure()
    # plt.title("The Hodgkin Huxley Neuron under a constant current of 10 $\mu A / cm^2 $. ( $Tau_{m}$*" + str(factor_m) + ", $Tau_{h}$*" + str(factor_h) + ", $Tau_{n}$*" + str(factor_n) + ")." )
    # plt.title("The Hodgkin-Huxley neuron under a discrete current pulse stimulus.")
    # plt.title("( $Tau_{m}$/" + str(Phi(T)) + ", $Tau_{h}$/" + str(Phi(T)) + ", $Tau_{n}$/" + str(Phi(T)) + "). $I_{0}$ = " + str(I_0) + ". $T$ (Celcius) = " + str(T) )

    plt.title("( $Tau_{m}$*" + str(factor_m) + ", $Tau_{h}$*" + str(factor_h) + ", $Tau_{n}$*" + str(
         factor_n) + "). $I_{0}$ = " + str(I_0) + ". $T$ (Celcius) = " + str(T))

    plt.plot(t, V, label="$V$")
    plt.grid(b = True, which="major", color="b", linestyle="-")
    plt.plot(t, I_ext(t, I_0), "r-", label="$I_{ext}$")
    plt.legend()
    plt.show()

    amplitudes = [1.55, 1.65, 1.75, 1.85, 1.95]





def PlotGating():
    """
    Function to plot the gating variables. Created to make plotting easier.
    :return: Plots for the gating variables
    """
    plt.figure()
    plt.title("Plots of $m$, $h$, $n$")


    plt.subplot(4, 1, 1)
    plt.plot(t, m, "r-", label="$m$")
    plt.plot(t, h, "b-", label="$h$")
    plt.plot(t, n, "g-", label="$n$")
    plt.xlabel("time")
    plt.ylabel("Gating Variable")
    plt.legend()
    plt.grid()


    plt.subplot(4, 1, 2)
    plt.plot(V, m, "r-", label="$V-m$")
    plt.xlabel("V")
    plt.ylabel("m")
    plt.grid()

    plt.subplot(4, 1, 3)
    plt.plot(V, h, "b-", label="$V-h$")
    plt.xlabel("V")
    plt.ylabel("h")
    plt.grid()

    plt.subplot(4, 1, 4)
    plt.plot(V, n, "g-", label="$V-h$")
    plt.xlabel("V")
    plt.ylabel("n")
    plt.grid()

    plt.show()

def PlotGating2():
    """
    Function to plot m vs h vs n to check the nature of their relationship
    """
    plt.figure()
    plt.title("Plots of the Gating Variables")

    plt.subplot(3,1,1)
    plt.plot(m, h, "g-", label="$m-h$")
    plt.xlabel("m")
    plt.ylabel("h")
    plt.grid()


    plt.subplot(3,1,2)
    plt.plot(n, h, "r-", label="$n-h$")
    plt.xlabel("n")
    plt.ylabel("h")
    plt.grid()

    plt.subplot(3,1,3)
    plt.plot(m, n, "b-", label="$m-n$")
    plt.xlabel("m")
    plt.ylabel("n")
    plt.grid()
    plt.show()


def PlotCurrent():
    """
    Function to plot the gating Currents. Created to make plotting easier.
    :return: Plots for the gating Currents
    """
    plt.figure()

    plt.subplot(3, 1, 1)
    plt.title("Plots of $I_{Na}$, $I_{K}$, $I_{L}$")
    plt.plot(t, INa, "g-", label="$I_{Na}$")
    plt.grid()

    plt.subplot(3, 1, 2)
    plt.plot(t, IK, "r-", label="$I_{K}$")
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(t, IL, "b-", label="$I_{L}$")
    plt.grid()

    plt.show()

"""--------------------------------------------------------------------------------------------------------------------
    5. Plotting the Time Constants and the steady state variables  
--------------------------------------------------------------------------------------------------------------------"""
m0 = []
h0 = []
n0 = []
for i in V:
    temp = steady_state(i)
    # a is now the array of steady state m's for all values of V
    a = temp[1]
    b = temp[2]
    c = temp[3]
    m0.append(a)
    h0.append(b)
    n0.append(c)


def PlotTime():
    """
    Function to plot the time constants of the gating variables
    :return: Plots
    """

    tau_m = Tau_m(V)
    tau_h = Tau_h(V)
    tau_n = Tau_n(V)

    plt.figure()

    #Plotting the steady state gating variables
    plt.subplot(2, 1, 1)
    plt.title("$m_0$, $h_0$ and $n_0$ with V for T = " + str(T))
    plt.plot(V, m0, "r-", label="$m_0$")
    plt.plot(V, h0, "b-", label="$h_0$")
    plt.plot(V, n0, "g-", label="$n_0$")
    plt.xlabel("V $mV$")
    plt.ylabel("Steady State Gating Variables")
    plt.legend()
    plt.grid()

    #Potting the Time Constants of the gating variables
    plt.subplot(2, 1, 2)
    plt.title("$Tau_m$, $Tau_h$ and $Tau_n$ with V for T = " + str(T))
    plt.plot(V, tau_m, "r-", label="$Tau_m$")
    plt.plot(V, tau_h, "b-", label="$Tau_h$")
    plt.plot(V, tau_n, "g-", label="$Tau_n$")
    plt.xlabel("V $mV$")
    plt.ylabel("Time Constants")
    plt.legend()
    plt.grid()

    plt.show()




"""---------------------------------------------------------------------------------------------------------------------
Code Written by Agastya Patri. 
---------------------------------------------------------------------------------------------------------------------"""




