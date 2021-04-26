MODULE plume_rise_IS4FIRES
  !
  ! This is the module for the plume-rise model of IS4FIRES
  ! Adopted with minor changes from the SILAM code.
  ! 
  ! Author: M.Sofiev, 2012
  ! Adaptation: M.Sofiev, 2020
  ! Language: FORTRAN-90, well, close to that
  !
  
  implicit none
  
  private
  
  public IS4FIRES_vertical_profile
  
CONTAINS
  
  !*************************************************************
  
  subroutine IS4FIRES_vertical_profile(FRP, &
                                       & fABL_height, fBruntVaisFreq, &
                                       & ifOneStepProcedure, &
                                       & fLevBottom, fLevTop)
    implicit none

    ! Imported parameters
    real, intent(in) :: FRP, fABL_height, fBruntVaisFreq
    logical, intent(in) :: ifOneStepProcedure
    real, intent(out) :: fLevBottom, fLevTop ! plume bounds

    ! Local parameters
    real, parameter :: alpha_full=0.24, alpha_step_1 = 0.15, alpha_step_2 = 0.93, &
                     & beta_full=169,   beta_step_1 = 102.,  beta_step_2 = 298., &
                     & gamma_full=0.35, gamma_step_1 = 0.49, gamma_step_2 = 0.13, &
                     & delta_full=254.,                      delta_step_2 = 259., &
                     & FRP_scale=1.e-6
    ! Local variables
    real :: fABL_height_local
    
    fABL_height_local = fABL_height
    if(fABL_height < 10. .or. fABL_height > 6000.)then
      print *, 'Funny ABL height:', fABL_height
      fABL_height_local = max(10.,min(6000.,fABL_height))
      print *, 'Funny ABL height, force'
    endif

    if(ifOneStepProcedure)then
      !
      ! One-step procedure means a unified formula
      !
      fLevTop = alpha_full * fABL_height_local + &
              & beta_full * ((FRP * FRP_scale) **gamma_full) * exp(-delta_full * fBruntVaisFreq)
      fLevBottom = fLevTop / 3.
    else
      !
      ! For two-steps, first compute the separation value, which is to be compared with ABL
      ! Then - either FT- or unified formula (the later is also OK for ABL)
      !
      fLevTop = alpha_step_1 * fABL_height_local + beta_step_1 * ((FRP * FRP_scale) **gamma_step_1)
      if(fLevTop > fABL_height_local)then
        fLevTop = alpha_step_2 * fABL_height_local + &
                & beta_step_2 * ((FRP * FRP_scale) **gamma_step_2) * &
                & exp(-delta_step_2 * fBruntVaisFreq)
        fLevBottom = fLevTop / 2.
      else
        fLevTop = alpha_full * fABL_height_local + &
                & beta_full * ((FRP * FRP_scale) **gamma_full) * exp(-delta_full * fBruntVaisFreq)
        fLevBottom = fLevTop / 3.
      endif
    endif

  end subroutine IS4FIRES_vertical_profile

end MODULE plume_rise_IS4FIRES
  