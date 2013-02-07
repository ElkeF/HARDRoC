MODULE trdm_library
! Purpose: This module holds transition dipole moments molecules
! to be read for the ETMD calculation
! The nomenclature is Accptor|Donor|J_A|M_A|J_D|M_D facotr, alpha or const
! for the case f(x) = factor * exp(-alpha*x) + const

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%% ArXe %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! J_A=1, M_A=1, J_D=3, M_D=3
REAL, PARAMETER :: ArXe_1_1_3_3_factor = 4.627300
REAL, PARAMETER :: ArXe_1_1_3_3_alpha  = 1.692111
REAL, PARAMETER :: ArXe_1_1_3_3_const  = -8E-6
! J_A=1, M_A=1, J_D=3, M_D=1
REAL, PARAMETER :: ArXe_1_1_3_1_factor = 36.228524
REAL, PARAMETER :: ArXe_1_1_3_1_alpha  = 1.728164
REAL, PARAMETER :: ArXe_1_1_3_1_const  = -2.6E-5
! J_A=1, M_A=1, J_D=3, M_D=-1
REAL, PARAMETER :: ArXe_1_1_3_m1_factor = 0.753757
REAL, PARAMETER :: ArXe_1_1_3_m1_alpha  = 1.759957
REAL, PARAMETER :: ArXe_1_1_3_m1_const  = -3E-6
! J_A=1, M_A=1, J_D=1, M_D=1
REAL, PARAMETER :: ArXe_1_1_1_1_factor = 46.461521
REAL, PARAMETER :: ArXe_1_1_1_1_alpha  = 1.734922
REAL, PARAMETER :: ArXe_1_1_1_1_const  = -5.8E-5
! J_A=1, M_A=1, J_D=1, M_D=-1
REAL, PARAMETER :: ArXe_1_1_1_m1_factor = 15.192572
REAL, PARAMETER :: ArXe_1_1_1_m1_alpha  = 1.825265
REAL, PARAMETER :: ArXe_1_1_1_m1_const  = -8E-6

END MODULE trdm_library
