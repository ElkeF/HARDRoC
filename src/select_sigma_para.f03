MODULE select_sigma_para

CONTAINS

SUBROUTINE select_sigma_fit_para(quad,lin,const,oneover)
! Purpose: select the correct parameters of an exponential fit of
! the transition dipole moments of molecules.

  use control
  use sigma_library

  IMPLICIT NONE

! Data dictionary: output parameters
  REAL, INTENT(OUT) :: quad,lin,const,oneover

! Data dictionary: 

  inatom: SELECT CASE (fin_atom_type2)
  
    CASE('Xe')
      quad    = Xe5p_quad
      lin     = Xe5p_lin
      const   = Xe5p_const
      oneover = Xe5p_oneover
    CASE DEFAULT
      WRITE(of,*) 'No ionization cross section provided for this atom type'

  END SELECT inatom

END SUBROUTINE select_sigma_fit_para

END MODULE select_sigma_para
