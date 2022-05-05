"""
    Simulating the Hodgkin-Huxley model of Action Potential under a sinusoidal input current.
    Written by Agastya Patri.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, ode
np.set_printoptions(threshold=np.inf)



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

init_cond = steady_state(-25)[0], steady_state(-25)[2], steady_state(-25)[3]

# Parameterskanye we
t = np.linspace(0, 5000, 100000)
factor_m = 1/10e10
I_0 = 2.75
factor_h = 4
factor_n = 1

def m_infinity(V):
    m_inf = float(alpha_m(V) / (alpha_m(V) + beta_m(V)))
    return m_inf


def HHModel_Reduced(X, t):
    """
    Hodgkin Huxley Model with instantaneous m dynamics
    """
    V = X[0]
    h = X[1]
    n = X[2]
    m_inf = m_infinity(V)

    I = lambda t, I_0: I_0 * np.sin(2 * np.pi * 0.05 * t)

    dVdt = ( I(t, I_0) - I_Na(V, m_inf, h) - I_K(V, n) - I_L(V) )/(C_m)

    dhdt = (alpha_h(V) * (1.0 - h) - beta_h(V) * h) / ( (Tau_h(V)*factor_h) * (alpha_h(V) + beta_h(V)))

    dndt = (alpha_n(V) * (1.0 - n) - beta_n(V) * n) / (Tau_n(V)*factor_n * (alpha_n(V) + beta_n(V)))

    return [dVdt, dhdt, dndt]



params = g_Na,g_K,g_L,V_Na,V_K,V_L



def RunModel():
    """
    A function to integrate the equations of the HH model which has been relieved of the m gating variable.
    :return: The integrated equations, as well as the gating currents.
    """
    Y = odeint(HHModel_Reduced, init_cond, t)

    V = Y[:,0]
    h = Y[:,1]
    n = Y[:,2]

    return V, h, n

V, h, n = RunModel()



# Code for the Bifurcation

def Bifurcation(amp_start, amp_end, step_size):
    """
    Function to plot return the Bifurcation dynamics of the Hodgkin Huxley Model.
    """

    I_0 = amp_start
    X = []
    Y = []

    while I_0 <= amp_end:
        I_arr = I_ext(t, I_0)
        Ys = RunModel()
        V = Ys[0][32000:]

        for i in range(100, len(V), 400):
            X.append(I_arr[i])
            Y.append(V[i])

        I_0 += step_size

    X = np.array(X)
    Y = np.array(Y)

    xdata = np.savetxt(str(amp_start) + "-to-" + str(amp_end) + "_x.csv", X)
    ydata = np.savetxt(str(amp_start) + "-to-" + str(amp_end) + "_y.csv", Y)
    # return X, Y





"""---------------------------------------------------------------------------------------------------------------------
   4. Plotting the Results
---------------------------------------------------------------------------------------------------------------------"""


def PlotResults(I_0):

    """
    Function to plot the integrated functions. Uncomment the relevant lines of code to plot what is needed
    :return: Responses of the variables with time, voltage.
    """
    V2 = V[32000:]
    t2 = t[32000:]
    I = I_ext(t, I_0)[32000:]


    plt.figure()
    # plt.title("The Hodgkin Huxley Neuron under a constant current of 10 $\mu A / cm^2 $. ( $Tau_{m}$*" + str(factor_m) + ", $Tau_{h}$*" + str(factor_h) + ", $Tau_{n}$*" + str(factor_n) + ")." )
    # plt.title("The Hodgkin-Huxley neuron under a discrete current pulse stimulus.")
    # plt.title("( $Tau_{m}$/" + str(Phi(T)) + ", $Tau_{h}$/" + str(Phi(T)) + ", $Tau_{n}$/" + str(Phi(T)) + "). $I_{0}$ = " + str(I_0) + ". $T$ (Celcius) = " + str(T) )

    plt.title("( $Tau_{m}$*" + str(factor_m) + ", $Tau_{h}$*" + str(factor_h) + ", $Tau_{n}$*" + str(
         factor_n) + "). $I_{0}$ = " + str(I_0) + ". $T$ (Celcius) = " + str(T))

    plt.plot(t2, V2, label="$V$")
    # plt.plot(t, I_ext(t, I_0), "r-", label="$I_{ext}$")
    plt.plot(t2, I, "r-", label="$I_{ext}$")
    plt.grid(b = True, which="major", color="b", linestyle="-")
    plt.grid(b=True, which="minor", color="b", linestyle="-", alpha=0.2)
    plt.minorticks_on()
    plt.legend()
    plt.show()




def PlotBifurcation(amp_start, amp_end):
    """
    Function to plot the data collected around the bifurcation region.
    """
    # PATH = "/home/agastya123/PycharmProjects/ComputationalNeuroscience/HodgkinHuxleyModel/Figures_and_Results/Bifurcation/Changed Params/Values_Around_Bifurcation/"
    # PATH = "/home/agastya123/Downloads/"
    PATH = "/home/agastya123/PycharmProjects/ComputationalNeuroscience/HodgkinHuxleyModel/Figures_and_Results/Bifurcation_Reduced_Model/Current_Ranges/"


    file_x = f"{amp_start}-to-{amp_end}_x.csv"
    file_y = f"{amp_start}-to-{amp_end}_y.csv"

    x = np.genfromtxt(PATH + file_x, delimiter=",")
    y = np.genfromtxt(PATH + file_y, delimiter=",")

    plt.figure(figsize=(10,10))
    plt.title(f"Bifurcation Diagram in the region {amp_start} to {amp_end} for the reduced model")
    plt.scatter(x, y, s=5)
    plt.grid(b=True, which="major", color="b", linestyle="-")
    plt.grid(b=True, which="minor", color="b", linestyle="-", alpha=0.2)
    plt.minorticks_on()
    plt.show()





"""---------------------------------------------------------------------------------------------------------------------
   6. Bifurcation for ranges of amplitude values
---------------------------------------------------------------------------------------------------------------------"""



amp_start = 2.801
amp_end = 3
I_0 = amp_start
X = []
Y = []

while I_0 <= amp_end:
    I_arr = I_ext(t, I_0)
    Ys = RunModel()
    V = Ys[0][32000:]
    for i in range(100, len(V), 400):
        X.append(I_arr[i])
        Y.append(V[i])

    I_0 += 0.001

X = np.array(X)
Y = np.array(Y)

xdata = np.savetxt(str(amp_start) + "-to-" + str(amp_end) + "_x.csv", X)
ydata = np.savetxt(str(amp_start) + "-to-" + str(amp_end) + "_y.csv", Y)

print(X)
print(Y)




"""---------------------------------------------------------------------------------------------------------------------
Code Written by Agastya Patri. 
---------------------------------------------------------------------------------------------------------------------"""


