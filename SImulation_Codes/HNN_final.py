"""
    Simulating the Hodgkin Huxley Model of action potential
    Written by Agastya Patri
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp

class HodgkinHuxley:
    def __init__(self, amplitude, temperature, frequency, factor_m, factor_h, factor_n):

        self.I_0 = amplitude
        self.T = temperature
        self.omega = frequency
        self.factor_m = factor_m
        self.factor_h = factor_h
        self.factor_n = factor_n


        self.C_m = 1.0
        self.V_Na = 115.0
        self.V_K = -12.0
        self.V_L = 10.5995
        self.g_Na = 120.0
        self.g_K = 36.0
        self.g_L = 0.3

        # self.t = np.arange(0, 450, 0.001)


    """-----------------------------------------------------------------------------------------------------------------
    1. Defining the Gating Currents and Temperature Dependence 
    -----------------------------------------------------------------------------------------------------------------"""

    def Phi(self):
        """
        Temperature dependance
        :return:
        """
        phi = 3**((self.T-6.3)/10)
        return phi

    def I_ext(self, time):
        current_sinusoidal = self.I_0 * np.sin(2 * np.pi * self.omega * time)
        current_constant = 10
        current_pulse = 10*(time>100) - 10*(time>200) + 35*(time>300) - 35*(time>400)
        return current_sinusoidal

    def gating_currents(self, V, h, b, r):
        """
        Function to define the gating currents used in the 2 Dimensional Model
        """
        v0, m0, h0, n0 = self.steady_state(-25)

        current_Na = self.g_Na* m0**3 * h * (V - self.V_Na)
        current_K = self.g_K* ((b-r*h)**4) * (V - self.V_K)
        current_L = self.g_L * (V - self.V_L)
        return [current_Na, current_K, current_L]


    """-----------------------------------------------------------------------------------------------------------------
    2. Defining Channel Gating Kinetics
    -----------------------------------------------------------------------------------------------------------------"""

    def alphas(self, V):
        a_m = self.Phi() * 0.1 * (25 - V) / (np.exp((25 - V) / 10.0) - 1)
        a_h = self.Phi() * 4.0 * np.exp(- V / 18.0)
        a_n = self.Phi() * 0.01 * (10.0 - V) / (np.exp((10.0 - V) / 10.0) - 1)
        return a_m, a_h, a_n

    def betas(self, V):
        b_m = self.Phi() * 4.0 * np.exp(- V / 18.0)
        b_h = self.Phi() * 1.0 / (np.exp((30.0 - V) / 10.0) + 1.0)
        b_n = self.Phi() * 0.125 * np.exp(-(V) / 80.0)
        return b_m, b_h, b_n

    def steady_state(self, v0):
        alpha_m, alpha_h, alpha_n = self.alphas(v0)
        beta_m, beta_h, beta_n = self.betas(v0)

        m0 = alpha_m / (alpha_m + beta_m)
        h0 = alpha_h / (alpha_h + beta_h)
        n0 = alpha_n / (alpha_n + beta_n)
        return [v0, m0, h0, n0]

    """-----------------------------------------------------------------------------------------------------------------
    3. Defining External Current and Time constants  
    -----------------------------------------------------------------------------------------------------------------"""

    def Taus(self, V):
        alpha_m, alpha_h, alpha_n = self.alphas(V)
        beta_m, beta_h, beta_n = self.betas(V)

        tau_m = 1 / (self.Phi() * (alpha_m + beta_m))
        tau_h = 1 / (self.Phi() * (alpha_h + beta_h))
        tau_n = 1 / (self.Phi() * (alpha_h + beta_h))

        return tau_m, tau_h, tau_n


    """-----------------------------------------------------------------------------------------------------------------
    3. Defining and Integrating the Hodgkin Huxley Model
    -----------------------------------------------------------------------------------------------------------------"""

<<<<<<< HEAD

    def HHModel(self, X, t):
=======
    def HHModel(self, t, X):
>>>>>>> origin/Reduced_Model_Bifurcation
        """
        Function to define the Hodgkin Huxley Model
        """
        V = X[0]
        m = X[1]
        h = X[2]
        n = X[3]

        # Calling the parameters

        alpha_m, alpha_h, alpha_n = self.alphas(V)
        beta_m, beta_h, beta_n = self.betas(V)

        Tau_m, Tau_h, Tau_n = self.Taus(V)

        I_Na, I_K, I_L = self.gating_currents(V, h, 0.95, 0.91)

        # Equations of the Model

        dVdt = (self.I_ext(t) - I_Na - I_K - I_L) * (1/self.C_m)

        dmdt = (alpha_m * (1.0 - m) - beta_m * m) / ((Tau_m * self.factor_m) * (alpha_m + beta_m))

        dhdt = (alpha_h * (1.0 - h) - beta_h * h) / ((Tau_h * self.factor_h) * (alpha_h + beta_h))

        dndt = (alpha_n * (1.0 - n) - beta_n * n) / ((Tau_n * self.factor_n) * (alpha_n + beta_n))

        return [dVdt, dmdt, dhdt, dndt]


    def RunModel(self):

        init_cond = self.steady_state(-25)

        # t = np.arange(0, 450, 0.001)

        Y = odeint(self.HHModel, init_cond, t)

        # INa, IK, IL = self.gating_currents(V, h, 0.95, 0.91)

<<<<<<< HEAD
        return Y
=======
        return Y.y
>>>>>>> origin/Reduced_Model_Bifurcation


    """-----------------------------------------------------------------------------------------------------------------
    3. Plotting the Results
    -----------------------------------------------------------------------------------------------------------------"""

    def PlotResults(self):
        """
        Function to plot the Action potential and the external current with time
        """
        V = self.RunModel()


        plt.figure(figsize=(10, 8))
<<<<<<< HEAD
        plt.title("2 Dimensional Hodgkin Huxley Model with Amplitude = " + str(self.I_0) )
        plt.plot(t, V[:,0], "g-", label="$t-V$")
        plt.plot(t, self.I_ext(time=t), "r-", label="$t-I_{ext}$")
=======
        plt.title("Hodgkin Huxley Model with Amplitude = " + str(self.I_0) )
        plt.plot(t, V, "g-", label="$t-V$")
        plt.plot(t, self.I_ext(t), "r-", label="$t-I_{ext}$")
>>>>>>> origin/Reduced_Model_Bifurcation
        plt.xlabel("time (ms)")
        plt.grid(b=True, which="major", color="b", linestyle="-")
        plt.grid(b=True, which="minor", color="b", linestyle="-", alpha=0.2)
        plt.minorticks_on()
        plt.legend()
        plt.show()


    def PlotGating(self):
        """
        Function to plot the gating variables with time
        """
        V, h, INa, IK, IL = self.RunModel()

        t = np.arange(0, 450, 0.001)
        plt.figure(figsize=(10, 8))
        plt.title("Gating Variables of the 2 D HH Model")
        plt.plot(h, n, "g-", label="$h-t$")
        # plt.plot(t, self.I_ext(time=t), "r-", label="$t-I_{ext}$")
        plt.xlabel("time (ms)")
        plt.grid(b=True, which="major", color="b", linestyle="-")
        plt.legend()
        plt.show()

    def PlotCurrent(self):
        """
        Function to plot the gating currents of the HH Model
        """
        V, h, INa, IK, IL = self.RunModel()

        t = np.arange(0, 450, 0.001)
        plt.figure(figsize=(10, 8))
        plt.title("Gating Currents")
        plt.plot(t, INa, "g-", label="$I_{Na}-t$")
        plt.plot(t, IK, "r-", label="$I_{K}-t$")
        plt.plot(t, IL, "b-", label="$I_{L}-t$")

        plt.xlabel("time (ms)")
        plt.ylabel("Current ($\mu A / cm^2$)")
        plt.grid(b=True, which="major", color="b", linestyle="-")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    t = np.linspace(0, 5000, 500000)
    neuron1 = HodgkinHuxley(amplitude=1.85, frequency=0.05, temperature=6.3, factor_m=1, factor_h=1, factor_n=1)


