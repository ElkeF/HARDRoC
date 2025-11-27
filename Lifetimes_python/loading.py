import constants
import numpy as np
import os
import os.path
import re

def get_SIP(SIP_file):
    #Reads in the Single Ionization Energies(SIP) from a file, 
    #with energies in the first column and polestrengths in the
    #second. Stores them in an array.

    SIP_data = np.loadtxt(SIP_file)

    SIP_energies = SIP_data[:,0] #[in eV]
    SIP_polestrengths = SIP_data[:,1] #in a.u.

    return SIP_energies

def get_DIP(DIP_file, SIP_file):
    #Reading in the Doubleionisation energies (DIP) from a file,
    #with energies in the first column and polestrengths in the second
    #and returns only those energies which are smaller than the SIP
    
    DIP_data = np.loadtxt(DIP_file)

    DIP_energies = DIP_data[:,0] #[in eV]
    DIP_pole_strength = DIP_data[:,1] #in a.u.

    filtered_DIP_energies = DIP_energies[DIP_energies < max(get_SIP(SIP_file))]

    return filtered_DIP_energies
