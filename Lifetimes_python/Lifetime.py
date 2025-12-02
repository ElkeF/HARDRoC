import numpy as np

class Lifetime:

  def __init__(self, SIP_energies, DIP_energies, xs_dict, rec_tau, n_in_val_orbs, n_out_val_orbs):
    """
    :SIP_energies: the single ionization energies (in eV)
    :DIP_energies: the double ionization energies (in eV)
    :xs_dict: a dictionary with the photoionization cross-sections (in Mb) 
    :rec_tau: an array with the reciprocal radiative lifetimes (in s**(-1))
    :n_in_val_orbs: number of inner valence orbitals
    :n_out_val_orbs: number of outer valence orbitals
    """

  def gamma_comp():
    """
    calculates the lifetime (gamma) for all SIPs

    gamma: lifetime of ICD for defined parameters in asymptotic approximation
    """
    return gamma

  def gamma_omega(self, omega, rec_tau_A, PI_xs_B):
    """
    calculates the lifetime for one specific omega

    :omega: energy difference between final and initial state
    :rec_tau_A: reciprocal radiative lifetime of the transition in monomer A
    :PI_xs_B: photoionization cross-section of monomer B for the specified omega
    """
    return gamma_omega
  
