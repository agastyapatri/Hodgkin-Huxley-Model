"""
    Rewriting the Hodgkin Huxley simulation code within the paradigm of the Object Oriented
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
np.set_printoptions(threshold=np.inf)


PATH = "/home/agastya123/PycharmProjects/ComputationalNeuroscience/HodgkinHuxleyModel/Figures_and_Results/Bifurcation_Reduced_Model/Current_Ranges/final_data/"

amp_array = np.genfromtxt("/home/agastya123/PycharmProjects/ComputationalNeuroscience/HodgkinHuxleyModel/Figures_and_Results/Bifurcation_Reduced_Model/Current_Ranges/final_data/2-to-3_x.csv")
relevant_array = amp_array[::170]

for i in range(len(relevant_array)):
    amp_array[i*170 : (i+1)*170] = relevant_array[i]


final_X = np.savetxt("2-to-3_x.csv", amp_array)







