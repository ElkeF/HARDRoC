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

def load_PI_xs(xs_file):
    """
    Reads in the energies in eV and photoionization cross-section in Mb from a file
    and places these into a dictionary
    """
    xs_data = np.loadtxt(xs_file)

    xs_energies = xs_data[:,0] #in eV
    xs = xs_data[:,1] #in Mb

    xs_dict = dict(zip(xs_energies,xs), dtype= float)
    
    return xs_dict
