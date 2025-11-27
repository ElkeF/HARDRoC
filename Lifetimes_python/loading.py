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
