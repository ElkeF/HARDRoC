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
REAL,PARAMETER :: c_au = 1/fs_const
REAL,PARAMETER :: pi = 3.141592654 !pi halt
REAL,PARAMETER :: a_0 = h/(2*pi*m_e*c*fs_const) !Bohr radius
REAL,PARAMETER :: E_h = fs_const**2 * m_e * c**2

! Conversion factors
! Distances, Areas, Volume
REAL,PARAMETER :: bohr_to_angstrom = a_0*1E10
REAL,PARAMETER :: angstrom_to_bohr = 1/a_0 * 1E-10
REAL,PARAMETER :: bohr_to_meter = a_0
REAL,PARAMETER :: meter_to_bohr = 1/a_0
REAL,PARAMETER :: megabarn_to_sqmeter = 1E-22
REAL,PARAMETER :: sqmeter_to_megabarn = 1E22
! Energies
REAL,PARAMETER :: joule_to_ev = 1/e
REAL,PARAMETER :: ev_to_joule = e
REAL,PARAMETER :: hartree_to_joule = h*c*fs_const/(a_0*2*pi)
REAL,PARAMETER :: joule_to_hartree = (a_0*2*pi)/(h*c*fs_const)
REAL,PARAMETER :: hartree_to_ev = hartree_to_joule*joule_to_ev
REAL,PARAMETER :: ev_to_hartree = ev_to_joule*joule_to_hartree
! Time
REAL,PARAMETER :: second_to_atu = E_h*2*pi/h
REAL,PARAMETER :: atu_to_second = h/(2*pi*E_h)

CONTAINS

SUBROUTINE write_phys_const()

  use control

  WRITE(of,*) 'Physical constants used in this programme'
  WRITE(of,131) 'Speed of light', c
  131 FORMAT(' ',A30,3X,ES20.13)
  WRITE(of,131) "Planck's constant", h
  WRITE(of,131) 'Rest mass of the electron', m_e
  WRITE(of,131) 'Electron charge', e
  WRITE(of,131) 'Fine structure constant', fs_const
  WRITE(of,131) 'pi', pi
  WRITE(of,131) 'Bohr radius', a_0
  WRITE(of,131) 'Hartree energy', E_h
  
  WRITE(of,*) ' '
  WRITE(of,*) 'Conversion factors used in this programme'
  WRITE(of,131) 'Angstrom to Bohr', angstrom_to_bohr
  WRITE(of,131) 'Bohr to Angstrom', bohr_to_angstrom
  WRITE(of,131) 'Megabarn to square meter', megabarn_to_sqmeter
  WRITE(of,131) 'Second to a.t.u', second_to_atu
  WRITE(of,131) 'a.t.u. to Second', atu_to_second
  WRITE(of,131) 'Joule to eV', joule_to_ev
  WRITE(of,131) 'eV to Joule', ev_to_joule
  WRITE(of,131) 'Joule to Hartree', joule_to_hartree
  WRITE(of,131) 'Hartree to Joule', hartree_to_joule
  WRITE(of,131) 'Hartree to eV', hartree_to_ev
  WRITE(of,131) 'eV to Hartree', ev_to_hartree
  WRITE(of,*) ' '
END SUBROUTINE write_phys_const

END MODULE physical_constants
