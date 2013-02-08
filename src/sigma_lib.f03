MODULE sigma_library
! Purpose: This module holds parameters for the calculation of
! ionization cross sections fitted to
! f(x) =  quad*x^2 + lin*x + const + oneover/x
! Here x is of course the energy of the virtual photon in eV.

! Xe p, combined 3/2 and 1/2
! Atomic and Nuclear Data Tables 22, 103 (1978)
  REAL,PARAMETER :: Xe5p_quad    =     0.695058
  REAL,PARAMETER :: Xe5p_lin     =   -41.741611
  REAL,PARAMETER :: Xe5p_const   =   773.748197
  REAL,PARAMETER :: Xe5p_oneover = -3655.360124

END MODULE sigma_library
