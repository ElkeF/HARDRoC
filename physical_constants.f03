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
REAL,PARAMETER :: a_0 = h/(2*pi*m_e*c*fs_const) !Bohr radius

! Conversion factors
REAL,PARAMETER :: angstrom_to_bohr = a_0*1E10
REAL,PARAMETER :: bohr_to_angstrom = 1/a_0 * 1E-10
! Energies
REAL,PARAMETER :: joule_to_ev = 1/e
REAL,PARAMETER :: ev_to_joule = e
REAL,PARAMETER :: hartree_to_joule = h*c*fs_const/(a_0*2*pi)
REAL,PARAMETER :: joule_to_hartree = (a_0*2*pi)/(h*c*fs_const)
REAL,PARAMETER :: hartree_to_ev = hartree_to_joule*joule_to_ev
REAL,PARAMETER :: ev_to_hartree = ev_to_joule*joule_to_hartree

CONTAINS
SUBROUTINE write_phys_const()
  WRITE(*,*) 'Physical constants used in this programme'
  WRITE(*,131) 'Speed of light', c
  131 FORMAT(' ',A30,3X,ES20.13)
  WRITE(*,131) "Planck's constant", h
  WRITE(*,131) 'Rest mass of the electron', m_e
  WRITE(*,131) 'Electron charge', e
  WRITE(*,131) 'Fine structure constant', fs_const
  WRITE(*,131) 'pi', pi
  WRITE(*,131) 'Bohr radius', a_0
  
  WRITE(*,*) ' '
  WRITE(*,*) 'Conversion factors used in this programme'
  WRITE(*,131) 'Angstrom to Bohr', angstrom_to_bohr
  WRITE(*,131) 'Bohr to Angstrom', bohr_to_angstrom
  WRITE(*,131) 'Joule to eV', joule_to_ev
  WRITE(*,131) 'eV to Joule', ev_to_joule
  WRITE(*,131) 'Joule to Hartree', joule_to_hartree
  WRITE(*,131) 'Hartree to Joule', hartree_to_joule
  WRITE(*,131) 'Hartree to eV', hartree_to_ev
  WRITE(*,131) 'eV to Hartree', ev_to_hartree
  WRITE(*,*) ' '
END SUBROUTINE write_phys_const

END MODULE physical_constants
