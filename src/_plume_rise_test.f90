program test_plume_rise
  !
  ! Test program for verifying technical implementation of the IS4FIRES and PRM
  ! plume-rise algorithms after their reprogramming.
  ! Creates a few fires, corresponding meteorology and calls the plume rise routines
  ! via the standard itnerface
  !
  use plume_rise_driver
  use rconstants
  
  implicit none
  
  ! Local variables
  type(Tfires) :: fires
  type(Tmeteo_data) :: meteo_data
  type(Tgrid_lonlat) :: meteo_grid
  type(Tvertical_hybrid) :: meteo_vertical
  integer :: iTmp, status, ix, iy, iz, izInverse
  logical :: eof
  real :: a, b, fABL, pressure
  character(len=100) :: strTmp

  !
  ! Create and initialise a bunch of fires
  !
  fires%nFires=5
  allocate(fires%lon(fires%nFires), fires%lat(fires%nFires), fires%FRP(fires%nFires), &
         & fires%burnt_area(fires%nFires), fires%mdur(fires%nFires), fires%moist(fires%nFires), &
         & fires%ifInsideGrid(fires%nFires))
  allocate(fires%inj_IS4FIRES(2, fires%nFires), fires%inj_PRM(2, fires%nFires))
  do iTmp = 1, fires%nFires
    fires%lon(iTmp) = 30.0 + iTmp
    fires%lat(iTmp) = 50.0 + iTmp
    fires%FRP(iTmp) = 1e6 * iTmp ** 2      ! W
    fires%burnt_area(iTmp) = 1e6 * iTmp    ! m2  MODIS pixel is 1x1 km2 or more
    fires%mdur(iTmp) = 3600.0              ! fire duration, sec
    fires%moist(iTmp) = 0.8                ! fuel moisture, fraction
  end do
  !
  ! Metadata for the run
  !
  meteo_grid = Tgrid_lonlat('lonlat', &
                          & 20., 40., &    ! sw_corner_lon, sw_corner_lat 
                          & 20, 20, &      ! nlon, nlat 
                          & 1.0, 1.0)      ! dlon_deg, dlat_deg
  !
  ! Vertical definition. Let's take the 137 levels from operational IFS
  !
  meteo_vertical%vertical_type = 'hybrid_sigma_pressure'
  meteo_vertical%nbr_of_levels = 137
  allocate(meteo_vertical%a_interface(meteo_vertical%nbr_of_levels), &
         & meteo_vertical%b_interface(meteo_vertical%nbr_of_levels), &
         & meteo_vertical%z_interface(meteo_vertical%nbr_of_levels))
  open(10, file='../dat_prm/referenceLevel_ecmwf_137.csv',status='old') !, mode='r', status='old')
  eof = .false.
  do while(.not.eof)
    iTmp = -1
    read(unit = 10, fmt = '(A)', iostat = status) strTmp
    strTmp = adjustl(strTmp)
    if(status /= 0)then
      eof = .true.
      cycle
    endif
    if(strTmp(1:1) == 'l' .or. strTmp(1:1) == '')cycle
    read(unit=strTmp, fmt=*) iTmp, a, b
    if(iTmp > 0)then
      meteo_vertical%a_interface(iTmp) = a
      meteo_vertical%b_interface(iTmp) = b
    endif
  end do
  !
  ! Now, the meteodata. Should give reasonable values
  !       TE = TEMPERATURE (K)
  !       PE = PRESSURE (kPa)
  !       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
  !
  allocate(meteo_data%data3d(meteo_vertical%nbr_of_levels, meteo_grid%nLon, meteo_grid%nLat, 6), &
         & meteo_data%data2d(meteo_grid%nLon, meteo_grid%nLat, 2), &
         & meteo_data%quantities_3d(6), meteo_data%quantities_2d(2))
  meteo_data%quantities_3d(1:6) = (/chU_wind, &      ! 'u', upe
                                  & chV_wind, &      ! 'v', vpe
                                  & chTempr, &       ! 't', te
                                  & chTheta, &       ! 'theta', the
                                  & chPressure, &    ! 'p', pe
                                  & chUhmidity/)     ! 'q', rhe
  meteo_data%quantities_2d(1:2) = (/chABL_qName, &        ! 'ABL_height'
                                  & chBruntVaisalaFreq_qName /)   !  'Brunt_Vaisala_freq'
  do iy = 1, meteo_grid%nLat
    do ix = 1, meteo_grid%nLon
      fABL = 100.0 + ix*iy
      meteo_data%data2d(ix,iy,1) = fABL   ! ABL height, m
      meteo_data%data2d(ix,iy,2) = 1e-3            ! Brunt-Vaisala frequency
      do iz = 1, meteo_vertical%nbr_of_levels
        meteo_data%data3d(iz,ix,iy,1) = 1.0   ! u, m/s
        meteo_data%data3d(iz,ix,iy,2) = 1.3   ! v, m/s
        izInverse = meteo_vertical%nbr_of_levels - iz + 1
        pressure = (meteo_vertical%a_interface(izInverse) + &
                  & meteo_vertical%b_interface(izInverse) * 1e5)
        meteo_data%data3d(iz,ix,iy,5) = pressure  ! Pa, drop 10%/km
        call temperatures(iz, fABL, pressure, meteo_vertical, &
                        & meteo_data%data3d(iz,ix,iy,3), meteo_data%data3d(iz,ix,iy,4))
        meteo_data%data3d(iz,ix,iy,6) = fu_q_prc(iz, fABL, meteo_vertical)
      end do
    end do
  end do
  !
  ! The main call. Dump file will be replicated to many in case of OMP run
  !
  call compute_plume_rise(fires, &                ! fires, input and output
                        & meteo_data, &           ! meteo data, input
                        & meteo_grid, meteo_vertical, &
                        & '../output/Dump.log_')
  !
  ! Print the results
  !
  
  CONTAINS
  
  !======================================================================

  subroutine temperatures(iz, fABL, pressure, meteo_vertical, T, Theta)
    implicit none
    integer, intent(in) :: iz
    real, intent(in) :: fABL, pressure
    type(Tvertical_hybrid), intent(inout) :: meteo_vertical
    real, intent(out) :: T, Theta
    
    real :: layer_thickness
    integer :: izInv

    if(iz == 1)then
      Theta = 300
      T = 300
      meteo_vertical%z_interface(iz) = 0.
      return
    endif

    izInv = meteo_vertical%nbr_of_levels - iz + 1
    if(meteo_vertical%z_interface(iz-1) <= fABL)then
      Theta = 300.0             ! neutral
    else
      Theta = 300 + 0.01 * (meteo_vertical%z_interface(iz-1) - fABL)  ! stable
    endif
    T = Theta / ((1e5/pressure)**0.2854)  ! = gas_constant_dryair/specific_heat_dryair)
    layer_thickness = rgas * T * &
        & LOG((meteo_vertical%b_interface(izInv+1) * 1e5 + meteo_vertical%a_interface(izInv+1)) / &
            & (meteo_vertical%b_interface(izInv) * 1e5 + meteo_vertical%a_interface(izInv))) / g
    meteo_vertical%z_interface(iz) = meteo_vertical%z_interface(iz-1)  + layer_thickness
  end subroutine temperatures
    
  !======================================================================
  
  real function fu_q_prc(iz, fABL, meteo_vertical)    
    implicit none
    integer, intent(in) :: iz
    real, intent(in) :: fABL
    type(Tvertical_hybrid), intent(inout) :: meteo_vertical

    if(meteo_vertical%z_interface(iz) <= fABL)then
      fu_q_prc = 50
    else
      fu_q_prc = 10
    endif
  end function fu_q_prc
  
  
end program test_plume_rise
  
