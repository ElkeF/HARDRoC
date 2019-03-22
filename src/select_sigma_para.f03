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
  
    CASE('Hx')
      quad    = He1s2p_quad
      lin     = He1s2p_lin
      const   = He1s2p_const
      oneover = He1s2p_oneover
    CASE('Ne')
      quad    = Ne2p_quad
      lin     = Ne2p_lin
      const   = Ne2p_const
      oneover = Ne2p_oneover
    CASE('Ar')
      quad    = Ar3p_quad
      lin     = Ar3p_lin
      const   = Ar3p_const
      oneover = Ar3p_oneover
    CASE('Kr')
      quad    = Kr4p_quad
      lin     = Kr4p_lin
      const   = Kr4p_const
      oneover = Kr4p_oneover
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
