MODULE trdm_library
! Purpose: This module holds transition dipole moments molecules
! to be read for the ETMD calculation
! The nomenclature is Accptor|Donor|J_A|M_A|J_D|M_D facotr, alpha or const
! for the case f(x) = factor * exp(-alpha*x) + const

! This is a test
REAL, PARAMETER :: ArXe_1_1_3_3_factor = 0.5
REAL, PARAMETER :: ArXe_1_1_3_3_alpha  = 0.4
REAL, PARAMETER :: ArXe_1_1_3_3_const  = 0.3

END MODULE trdm_library
