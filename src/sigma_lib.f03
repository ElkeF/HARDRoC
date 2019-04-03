MODULE sigma_library
! Purpose: This module holds parameters for the calculation of
! ionization cross sections fitted to
! f(x) =  quad*x^2 + lin*x + const + oneover/x
! Here x is of course the energy of the virtual photon in eV.

! ionization after excitation, state: He 1s2p
! Atomic and Nuclear Data Tables ,  (1976)
  REAL,PARAMETER :: He1s2p_quad    =  -0.036486
  REAL,PARAMETER :: He1s2p_lin     =   1.659328
  REAL,PARAMETER :: He1s2p_const   = -23.078388
  REAL,PARAMETER :: He1s2p_oneover = 100.082905

!! ionization after excitation, state: He 1s3p
!! Atomic and Nuclear Data Tables ,  (1976)
!! better for small energy differences
!  REAL,PARAMETER :: He1s3p_quad    =  -0.008797
!  REAL,PARAMETER :: He1s3p_lin     =   0.385338
!  REAL,PARAMETER :: He1s3p_const   =  -5.203295
!  REAL,PARAMETER :: He1s3p_oneover =  22.567650

! ionization after excitation, state: He 1s3p
! Atomic and Nuclear Data Tables ,  (1976)
  REAL,PARAMETER :: He1s3p_quad    =  -0.006761
  REAL,PARAMETER :: He1s3p_lin     =   0.319469
  REAL,PARAMETER :: He1s3p_const   =  -4.641940
  REAL,PARAMETER :: He1s3p_oneover =  21.261983

! Na 3s
! webcross sections elettra
  REAL,PARAMETER :: Na3s_quad    =     0.000053
  REAL,PARAMETER :: Na3s_lin     =    -0.005798
  REAL,PARAMETER :: Na3s_const   =     0.191587
  REAL,PARAMETER :: Na3s_oneover =     0.147936

! Ne p, combined 3/2 and 1/2
! Atomic and Nuclear Data Tables ,  (1976)
  REAL,PARAMETER :: Ne2p_quad    =     0.023484
  REAL,PARAMETER :: Ne2p_lin     =   -2.3761661
  REAL,PARAMETER :: Ne2p_const   =    90.936913
  REAL,PARAMETER :: Ne2p_oneover =  -957.652355

! Ar p, combined 3/2 and 1/2
! Atomic and Nuclear Data Tables ,  (1976)
  REAL,PARAMETER :: Ar3p_quad    =     0.076643
  REAL,PARAMETER :: Ar3p_lin     =    -9.891266
  REAL,PARAMETER :: Ar3p_const   =   344.970801
  REAL,PARAMETER :: Ar3p_oneover = -2820.658029

! Kr p, combined 3/2 and 1/2
! Atomic and Nuclear Data Tables ,  (1976)
  REAL,PARAMETER :: Kr4p_quad    =     0.614407
  REAL,PARAMETER :: Kr4p_lin     =   -43.033045
  REAL,PARAMETER :: Kr4p_const   =   977.917266
  REAL,PARAMETER :: Kr4p_oneover = -6479.795897

! Xe p, combined 3/2 and 1/2
! Atomic and Nuclear Data Tables 22, 103 (1978)
  REAL,PARAMETER :: Xe5p_quad    =     0.695058
  REAL,PARAMETER :: Xe5p_lin     =   -41.741611
  REAL,PARAMETER :: Xe5p_const   =   773.748197
  REAL,PARAMETER :: Xe5p_oneover = -3655.360124

END MODULE sigma_library
