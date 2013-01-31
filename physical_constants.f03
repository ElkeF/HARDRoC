MODULE physical_constants
!Purpose: This module holds physical constants and conversion factors
! from one unit to another one

IMPLICIT NONE

!Physical constants
REAL,PARAMETER :: c = 299792458 !Speed of light in m/s
REAL,PARAMETER :: h = 6.62606957E-34 !Plancks constant in SI
REAL,PARAMETER :: m_e = 9.10938215E-31 !Rest mass of electron
REAL,PARAMETER :: e = 1.602176565E-19 !Charge of the electron in SI
REAL,PARAMETER :: fs_const = 1/137.035999074 !alpa alias fine-structure constant
REAL,PARAMETER :: pi = 3.141592654 !pi halt
REAL,PARAMETER :: a_0 = h/(m_e*c*fs_const) !Bohr radius

! Conversion factors
!REAL,PARAMETER :: angstrom_to_bohr
!REAL,PARAMETER :: bohr_to_angstrom
! Energies
REAL,PARAMETER :: joule_to_ev = 1/e
REAL,PARAMETER :: ev_to_joule = e
REAL,PARAMETER :: joule_to_hartree = h*c*fs_const/(a_0*2*pi)
REAL,PARAMETER :: hartree_to_joule = (a_0*2*pi)/(h*c*fs_const)
REAL,PARAMETER :: hartree_to_ev = hartree_to_joule*joule_to_hartree
REAL,PARAMETER :: ev_to_hartree = ev_to_joule*joule_to_hartree

END MODULE physical_constants
