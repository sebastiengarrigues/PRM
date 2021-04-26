module plume_rise_driver
  !
  ! This is the module for the plume-rise model of IS4FIRES
  ! Adopted with minor changes from the SILAM code.
  ! 
  ! Author: M.Sofiev, 2012
  ! Adaptation: M.Sofiev, 2020
  ! Language: FORTRAN-90, well, close to that
  ! Units: SI basic units, coordinates in degrees [-180,180] and decimals
  ! FRP is in [W], etc.
  !
  use plume_rise_IS4FIRES
  use plume_rise_PRMv1
  use rconstants

  !$ use omp_lib
  
  implicit none
  
  private

  ! Public subrojtines from this module
  public compute_plume_rise
  
  ! private subroutines of this module
  private fu_index
  private extract_meteo_data
  
  ! Public constants from this module
  character(len=20), parameter, public :: chABL_qName = 'ABL_height'
  character(len=20), parameter, public :: chBruntVaisalaFreq_qName = 'Brunt_Vaisala_freq'
  character(len=20), parameter, public :: chU_wind = 'u'
  character(len=20), parameter, public :: chV_wind = 'v'
  character(len=20), parameter, public :: chTempr = 't'
  character(len=20), parameter, public :: chTheta = 'theta'
  character(len=20), parameter, public :: chPressure = 'p'
  character(len=20), parameter, public :: chUhmidity = 'q'  
  
  !
  ! Type for the grid definition. Should be passed together
  ! with the numerical fields of the meteodata
  !
  type Tgrid_lonlat
    character (len=6)  :: grid_type = 'lonlat'  ! just to be sure
    real*8 :: sw_corner_lon, sw_corner_lat    ! lon-lat coordinates of the first gridpoint given in lon-lat system
    integer :: nlon, nlat                     ! number of points along longitude / latitude
    real*8 :: dlon_deg, dlat_deg              ! Grid distance in west-east and south-north, degrees and decimals
  end type Tgrid_lonlat
  public Tgrid_lonlat
  !
  ! Type for the vertical definition. 
  ! Similar to the grid, should be passed together with the meteodata fields
  !
  type Tvertical_hybrid
    character(len=21) :: vertical_type = 'hybrid_sigma_pressure'  ! just to be sure
    integer :: nbr_of_levels                           ! the number of layers in the vertical
    real, dimension(:), allocatable :: a_interface, b_interface  ! hybrid coefficients of the upper interface of the layers
    real, dimension(:), allocatable :: z_interface  ! level  height above surface, [m]
  end type Tvertical_hybrid
  public Tvertical_hybrid
  !
  ! Fire data structure
  !
  type Tfires
    integer :: nFires                     ! the number of fires in the set
    real, dimension(:), allocatable :: lon, lat, FRP, burnt_area, mdur, moist   ! (nFires)
    real, dimension(:,:), allocatable :: inj_IS4FIRES, inj_PRM  ! (2,nFires) injection bottom and top, [m]
    logical, dimension(:), allocatable :: ifInsideGrid   ! for the case of regional application
  end type Tfires
  public Tfires
  !
  ! Meteodata structure
  !
  type Tmeteo_data
    real, dimension(:,:,:,:), allocatable :: data3d     ! (nz, nx, ny, nQauntities)
    real, dimension(:,:,:), allocatable :: data2d       ! (nx, ny, nQuantities)
    character(len=20), dimension(:), allocatable :: quantities_3d, quantities_2d
  end type Tmeteo_data 
  public Tmeteo_data

CONTAINS
  
  !*************************************************************
  
  subroutine compute_plume_rise(fires, &                ! fires, input and output
                              & meteo_data, &           ! meteo data, input
                              & meteo_grid, meteo_vertical, & ! metadata, input
                              & chDumpTemplate)
    !
    ! This is the main subroutine driving the plume-rise odules.
    ! It should be called with a list of fires and meteorological fields
    ! For each of the given fire (or a grid cell with the fire)
    ! the sub will return the top and bottomn of the plume
    !
    implicit none

    ! Imported parameters
    type(Tfires), intent(inout) :: fires
    type(Tmeteo_data), intent(in) :: meteo_data
    type(Tgrid_lonlat), intent(in) :: meteo_grid
    type(Tvertical_hybrid), intent(in) :: meteo_vertical
    character(len=*), intent(in) :: chDumpTemplate
    
    ! Local parameters.
    !
    ! Here, I put several supposedly constants, which should be reviewed
    ! and set their actual values
    !
    integer, parameter :: nTimeSteps_PRM = 200
    real, parameter :: alpha=0.1, C_epsi = 0.55 * 10, C_delta = -10.e0 * 5, &
                     & C_wind = 0.5e0, C_delta_wind = 1, &
                     & FRP2CHF=1.0       ! no idea what it should actually be
    logical, parameter :: wind_eff=.false., micro_eff=.false.
    
    
    ! Local variables
    real, dimension(:,:,:), allocatable :: meteo4fires_column  ! (nLevels, nQuantities, nFires)
    real, dimension(:,:), allocatable :: meteo4fires_sfc       ! (nQuantities, nFires)
    integer :: iABL, iBVFreq, iFire, uDump=50, iThread
    character(len=5) :: str5
    !
    ! First of all, get the data for the fire locations
    !
    allocate(meteo4fires_column(meteo_vertical%nbr_of_levels, 6, fires%nFires), &
           & meteo4fires_sfc(2, fires%nFires))
    call extract_meteo_data(fires, meteo_data, meteo_grid, meteo_vertical, &
                          & meteo4fires_column, meteo4fires_sfc)
    !
    ! Plume rise from IS4FIRES
    !
    iABL = fu_index(meteo_data%quantities_2d, chABL_qName)
    iBVFreq = fu_index(meteo_data%quantities_2d, chBruntVaisalaFreq_qName)
    !
    ! Scan all fires calling the corresponding plume-rise routines.
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iThread, iFire, &
    !$OMP & uDump, str5)
    iThread = 0
    !$ iThread = OMP_GET_THREAD_NUM()
    write(str5,fmt='(i5)')iThread
    str5 = adjustl(str5)

    open(uDump, file=chDumpTemplate // trim(str5))
    
    !$OMP DO schedule (guided)  ! or dynamic but for comparable-size loops guided is better
    do iFire = 1, fires%nFires
      !
      ! Plume rise from IS4FIRES.
      ! The model works in SI units
      !
      fires%inj_IS4FIRES = -1
      call IS4FIRES_vertical_profile(fires%FRP(iFire), &
                                   & meteo4fires_sfc(iABL, iFire), &
                                   & meteo4fires_sfc(iBVFreq, iFire), &
                                   & .true., &
                                   & fires%inj_IS4FIRES(1,iFire), fires%inj_IS4FIRES(2,iFire))
      write(uDump,*)'>>>>>>>>>>>> IS4FIRES bottom=', fires%inj_IS4FIRES(1,iFire), &
                                           & 'top=', fires%inj_IS4FIRES(2,iFire)
      !
      ! Plume rise from PRM
      ! PRM works in randomly picked units, watchout the conversion
      !
      fires%inj_PRM = -1
      call plumerise_PRMv1(meteo4fires_column(:,:,iFire), meteo_data%quantities_3d, &
                         & meteo4fires_sfc(:,iFire), meteo_data%quantities_2d, &  
                         & meteo_vertical, &
                         & fires%burnt_area(iFire) * 1e-4, &  ! from m2 to Ha
                         & fires%frp(iFire) * 1e-6, &         ! from W to MW
                         & fires%mdur(iFire) / 60., &         ! from sec to min
                         & fires%moist(iFire), &        ! remains fraction
                         & alpha, C_epsi, C_delta, C_wind, C_delta_wind,  &  ! fire configuration 
                         & FRP2CHF,  &
                         & wind_eff, micro_eff, &
                         & .true., .true.,   &   ! ifDebugPrint, ifStoreDump, 
                         & uDump, & ! debug/dump print file
                         & fires%inj_PRM(1,iFire), fires%inj_PRM(2,iFire))
      write(uDump,*)'>>>>>>>>>>>> IS4FIRES again bottom=', fires%inj_IS4FIRES(1,iFire), &
                                                 & 'top=', fires%inj_IS4FIRES(2,iFire)
      write(uDump,*)'>>>>>>>>>>>>>>>>>>>>>>> PRM bottom=', fires%inj_PRM(1,iFire), &
                                                 & 'top=', fires%inj_PRM(2,iFire)
    end do ! cycle over files
    !$OMP END DO

    !$OMP END PARALLEL
    
  end subroutine compute_plume_rise


  !*********************************************************                            
                              
  integer function fu_index(q_names, qIn)
    !
    ! Finds the index of the specific variable in the list of variable names
    ! Basically, a simple strong-comparison operation
    !
    implicit none
    
    ! Imported parameters
    character(len=*), dimension(:), intent(in) :: q_names
    character(len=*), intent(in) :: qIn
    
    ! Local variables
    integer :: idx
  
    fu_index = -1
    do idx = 1, size(q_names)
      if(qIn == q_names(idx))then
        fu_index = idx 
        return
      endif
    end do
    ! If we are here, something went wrong
    print *, 'Cannot find this variable: ', qIn
    print *, 'Available names:'
    do idx = 1, size(q_names)
      if(q_names(idx) == '')return
      print *, q_names(idx)
    end do
    
  end function fu_index

  
  !***************************************************************
  
  subroutine extract_meteo_data(fires, meteo_data, meteo_grid, meteo_vertical, &
                              & meteo4fires_column, meteo4fires_sfc)
    !
    ! Prepares meteodata for the fire plume rise models
    ! - picks the data from the grid
    ! - interpolates columns to the model vertical (for PRM model)
    !
    implicit none
    
    ! imported parameters
    type(Tfires), intent(inout) :: fires
    type(Tmeteo_data), intent(in) :: meteo_data
    type(Tgrid_lonlat), intent(in) :: meteo_grid
    type(Tvertical_hybrid), intent(in) :: meteo_vertical
    real, dimension(:,:,:), intent(out) :: meteo4fires_column  ! (nLevels, nFires, nQuantities)
    real, dimension(:,:), intent(out) :: meteo4fires_sfc       ! (nFires, nQuantities)
    
    ! Local variables
    integer :: iMet, iFire, ix, iy, izMet
    real :: xNew, yNew
    !
    ! Scan the fires one by one and pick the data
    !
    do iFire = 1, fires%nFires
      !
      ! Get 2d coordinates
      !
      xNew = 1 + (fires%lon(iFire) - meteo_grid%sw_corner_lon) / meteo_grid%dlon_deg
      yNew = 1 + (fires%lat(iFire) - meteo_grid%sw_corner_lat) / meteo_grid%dlat_deg

      ! If longitude out of grid, try rotate the earth once
      if(xNew < 0.5)then
        if(xNew + 360. / meteo_grid%dlon_deg < meteo_grid%nlon + 0.5) &
                                             & xNew = xNew + 360. / meteo_grid%dlon_deg
      elseif(xNew > meteo_grid%nlon + 0.5)then
        if(xNew - 360. / meteo_grid%dlon_deg > 0.5) &
                                             & xNew = xNew - 360. / meteo_grid%dlon_deg
      endif
      ! Are we inside the grid ?
      ix = nint(xNew)
      iy = nint(yNew)
      if(ix < 1 .or. ix > meteo_grid%nlon .or. iy < 1 .or. iy > meteo_grid%nlat)then
        fires%ifInsideGrid(iFire) = .false.
        cycle
      else
        fires%ifInsideGrid(iFire) = .true.
      endif
      !
      ! Extract data from 2d fields
      !
      do iMet = 1, size(meteo_data%quantities_2d)
        if(meteo_data%quantities_2d(iMet) == '')exit
        meteo4fires_sfc(iMet, iFire) = meteo_data%data2d(ix,iy,iMet)
      end do
      !
      ! Now, deal with the column. All vriables later will be interpolated 
      ! to the vertical of the plume-rise model. Here we only extract the column and store the 
      ! height of the levels
      !
      do iMet = 1, size(meteo_data%quantities_3d)  ! number of meteo quantities
        if(meteo_data%quantities_3d(iMet) == '')exit
        meteo4fires_column(:,iMet,iFire) = meteo_data%data3d(:,ix,iy,iMet)
      end do
    end do  ! nFires

  end subroutine extract_meteo_data


  !************************************************************************

  subroutine plumerise_PRMv1(meteo4fire_column, quantities_3d, &
                           & meteo4fire_sfc, quantities_2d, & 
                           & meteo_vertical, &
                           & burnt_area, frp_fire, mdur, moist, & 
                           & alpha,C_epsi,C_delta,C_wind,C_delta_wind, &  ! fire configuration 
                           & FRP2CHF, &
                           & wind_eff, micro_eff, &
                           & ifDebugPrint, ifStoreDump, uDump, &
                           & ztopmin, ztopmax)
    !
    ! Deals with a single fire case
    !-------------------------------------------------------------------------------------------!
    !
    ! Plume rise model for vegetation fires (CPTEC/INPE 2005-2006,2009)			    !
    ! Refs.:										    !
    ! Freitas, S. R., K. M. Longo, R. Chatfield, D. Latham, M. A. F. Silva Dias, M. O. Andreae, !
    ! E. Prins, J. C. Santos, R. Gielow and J. A. Carvalho Jr.: Including the sub-grid scale    !
    ! plume rise of vegetation fires in low resolution atmospheric transport models. 	    !
    !  Atmospheric Chemistry and Physics,2007.				                    !
    !-											    !
    ! Freitas, S. R.; Longo, K. M.; M. Andreae. Impact of including the plume rise of vegetation! 
    ! fires in numerical simulations of associated atmospheric pollutants. Geophys. Res. Lett., !
    ! 33, L17808, doi:10.1029/2006GL026608, 2006.                                               !
    !                                                                                           !
    ! questions/bugs/comments: Saulo Freitas (saulo.freitas@cptec.inpe.br) 			    !
    !-------------------------------------------------------------------------------------------!
    !
    implicit none
    
    ! Imported parameters
    real, dimension(:,:), intent(in) :: meteo4fire_column  ! (nLevels, nQuantities) no nFires dim
    real, dimension(:), intent(in) :: meteo4fire_sfc       ! (nQuantities)  no nFires dim
    character(len=*), dimension(:), intent(in) :: quantities_2d, quantities_3d
    type(Tvertical_hybrid), intent(in) :: meteo_vertical
    integer,intent(in) :: uDump
    real,intent(in) :: burnt_area, frp_fire, mdur, moist ! burnt_area in Ha, as it seems
    real,intent(in) :: alpha,C_epsi,C_delta,C_wind,C_delta_wind
    real,intent(in) :: FRP2CHF ! FRP2TOTALH,CHF2TOTALH
    logical,intent(in) :: wind_eff,micro_eff
    logical,intent(in) :: ifDebugPrint, ifStoreDump
    real,intent(out) :: ztopmax, ztopmin

    ! Local variables
    real :: heat_fluxW
    real, dimension(:), allocatable :: profile_mass_detrained_int_out, &
                                     & profile_mass_entrained_int_out
    
    type(TPRM_data) :: PRM_data  ! main dataset for the specific fire

    !
    ! initialize the main structure (only need to be done at the 1st time)
    ! Basic fire features are needed immediately
    !
    heat_fluxW = FRP2CHF * frp_fire * 100 / burnt_area  ! W/m2 from MW/Ha   ! scale form litterature
    !
    ! Set rest of the fire properties in the data structure
    !
    call pack_fire_properties(mdur, moist, & 
                            & alpha, C_epsi, C_delta, C_wind, C_delta_wind, &  ! fire configuration 
                            & FRP2CHF, &                 !FRP2TOTALH_in,CHF2TOTALH_in, &
                            & wind_eff, micro_eff, & 
                            & heat_fluxW, burnt_area, &
                            & ifDebugPrint, &
                            & PRM_data)
    !
    ! Brings the input meteodata into the computation vertical
    !
    call interpolate_env_profile(meteo4fire_column, quantities_3d, & ! (nzMeteo, nQuantities, nFires)
                               & meteo4fire_sfc, quantities_2d, &    ! (nQuantities, nFires)
                               & meteo_vertical%z_interface, &  ! interface height above the ground
                               & PRM_data)
    !
    ! Initialize the remainning pieces
    !
    call INITIAL(heat_fluxW, burnt_area * 1.e4, PRM_data, ifDebugPrint, ifStoreDump, uDump)
    !
    ! prepare the output of plumerise()
    !
    allocate(profile_mass_detrained_int_out(nzPRM))
    allocate(profile_mass_entrained_int_out(nzPRM))
    profile_mass_detrained_int_out(:) = 0.e0
    profile_mass_entrained_int_out(:) = 0.e0
    !
    ! The main PRM injection height model
    !
    call makeplume(PRM_data, heat_fluxW, &
                 & ifDebugPrint, ifStoreDump, uDump, & ! debug/dump print & files where to put it
                 & ztopmax, &
                 & profile_mass_detrained_int_out(1:nzPRM), &
                 & profile_mass_entrained_int_out(1:nzPRM))    

    ztopmin = ztopmax / 2.0   ! pretty much anything.

    if (ifDebugPrint .or. ifStoreDump) then 
        print*,'----------------------'
!        print*,' stop after ',mintime,' min'
        print*,' burnt_area (ha) - top cloud(km)- power (kW m^2)' 
        print*,burnt_area, ztopmax * 1.e-3 ,heat_fluxW*1.e-3
!        open(23,file='output.txt')
        write(uDump,*) '#burnt_area(ha)      FRP(MW)         top_cloud(km)      CHF(kWm^2)'
        write(uDump,'(4(f14.6,3x))') burnt_area, frp_fire, ztopmax/1000.,heat_fluxW*1.e-3
!        close(23)
    endif
  end subroutine plumerise_PRMv1


  !*******************************************************************

  subroutine interpolate_env_profile(meteo4fire_column, quantities_3d, &  
                                   & meteo4fire_sfc, quantities_2d, &
                                   & zMeteo, &               ! height of meteo levels, [m]
                                   & dat)
!  , ucon,vcon,urcon,prcon, tmpcon,dncon,zcon,thtcon,dimz)
    !
    ! Projects environmental profiles to the vertical of PRM
    !
    implicit none

    ! Imported parameters
    real, dimension(:,:), intent(in) :: meteo4fire_column  ! (nLevels, nQuantities)
    real, dimension(:), intent(in) :: meteo4fire_sfc       ! ( nQuantities)
    character(len=*), dimension(:), intent(in) :: quantities_2d, quantities_3d
    real, dimension(:), intent(in) :: zMeteo
    type(TPRM_data), intent(inout) :: dat
    
    ! Local variables
    integer :: k, nzMet
    real :: es
!    real, dimension(nzMet) :: ucon, vcon, thtcon, tmpcon, dncon, prcon, zcon, urcon
    !
    ! set the PRM grid
    !
    call set_grid() ! define vertical grid of plume model
    !
    ! interpolation of the environmental parameters
    !
    nzMet = size(meteo4fire_column, 1)
    
    ! Wind u and v
!    call htint(nzMet,  ucon, zMeteo, nzPRM, dat%upe, dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'u')), zMeteo, &
             & nzPRM, dat%upe, dat%zt)
!    call htint(nzMet,  vcon, zMeteo, nzPRM, dat%vpe, dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'v')), zMeteo, &
             & nzPRM, dat%vpe, dat%zt)
    ! Temperature and potential temperature
!    call htint(nzMet,tmpcon, zMeteo, nzPRM, dat%te , dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'t')), zMeteo, &
             & nzPRM, dat%te, dat%zt)
!    call htint(nzMet,thtcon, zMeteo, nzPRM, dat%the, dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'theta')), zMeteo, &
             & nzPRM, dat%the, dat%zt)
    ! Pressure, in Pa here
!    call htint(nzMet, prcon, zMeteo, nzPRM, dat%pe , dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'p')), zMeteo, &
             & nzPRM, dat%pe, dat%zt)
    ! Humidity?
!    call htint(nzMet, urcon, zMeteo, nzPRM, dat%rhe, dat%zt)
    call htint(nzMet, meteo4fire_column(:,fu_index(quantities_3d,'q')), zMeteo, &
             & nzPRM, dat%rhe, dat%zt)

    do k = 1, nzPRM
      !  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
      ES = fu_satur_water_vapour_kPa (dat%TE(k)) !blob saturation vapor pressure, kPa
      dat%QSAT (k) = (.622 * ES) / (dat%PE(k) * 1e-3 - ES) !saturation lwc g/g
      !calcula qvenv
      dat%qvenv(k) = max(0.01 * dat%rhe(k) * dat%QSAT(k), 1e-8)
      !   print*,k,dat%QSAT (k),rhe(k),qvenv(k),rhe(k)*dat%QSAT (k)

      dat%thve(k) = dat%the(k) * (1. + .61*dat%qvenv(k)) ! virtual pot temperature
      dat%dne(k) = dat%pe(k) *1e3 / (rgas * dat%te(k)*(1. + .61*dat%qvenv(k))) !  dry air density (kg/m3)
      dat%vel_e(k) = sqrt(dat%upe(k)**2 + dat%vpe(k)**2)

      call thetae(dat%pe(k), dat%te(k), dat%qvenv(k), dat%thee(k), dat%tde(k))

      !--------- converte press de Pa para kPa para uso modelo de plumerise
      dat%pe(k) = dat%pe(k) * 1.e-3

      dat%SCE(k) = 0.e0 !++ rp intialize the scalar in the environment
    enddo

    !++ rp get the high of the mixing layer form the maximum of the relative humidity.
    dat%hBL = meteo4fire_sfc(fu_index(quantities_2d, chABL_qName))

    CONTAINS

      !===========================================================

      subroutine set_grid()

        implicit none
        integer :: k

        !dz=50. ! set constant grid spacing of plume grid model(meters)
        dat%dz= dzPRM ! set constant grid spacing of plume grid model(meters)

        dat%zsurf = 0.e0
        dat%zt(1) = dat%zsurf
        dat%zm(1) = dat%zsurf
        dat%zt(2) = dat%zt(1) + 0.5 * dat%dz
        dat%zm(2) = dat%zm(1) + dat%dz
        do k=3,nzPRM
          dat%zt(k) = dat%zt(k-1) + dat%dz ! thermo and water levels
          dat%zm(k) = dat%zm(k-1) + dat%dz ! dynamical levels	
        enddo
        !print*,zsurf
        !Print*,zt(:)
        do k = 1, nzPRM-1
          dat%dzm(k) = 1. / (dat%zt(k+1) - dat%zt(k))
        enddo 
        dat%dzm(nzPRM) = dat%dzm(nzPRM-1)

        do k = 2, nzPRM
          dat%dzt(k) = 1. / (dat%zm(k) - dat%zm(k-1))
        enddo
        dat%dzt(1) = dat%dzt(2) * dat%dzt(2) / dat%dzt(3)
        !   dzm(1) = 0.5/dz
        !   dzm(2:nzPRM) = 1./dz
        return
      end subroutine set_grid

      !===========================================================

      SUBROUTINE htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)
        IMPLICIT NONE
        INTEGER, INTENT(IN ) :: nzz1
        INTEGER, INTENT(IN ) :: nzz2
        REAL,    INTENT(IN ) :: vctra(nzz1)
        REAL,    INTENT(OUT) :: vctrb(nzz2)
        REAL,    INTENT(IN ) :: eleva(nzz1)
        REAL,    INTENT(IN ) :: elevb(nzz2)

        INTEGER :: l
        INTEGER :: k
        INTEGER :: kk
        REAL    :: wt

        l=1
        DO k=1,nzz2
          DO
            IF ( (elevb(k) <  eleva(1)) .OR. &
              ((elevb(k) >= eleva(l)) .AND. (elevb(k) <= eleva(l+1))) ) THEN
              wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
              vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
              EXIT
            ELSE IF ( elevb(k) >  eleva(nzz1))  THEN
              wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
              vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
              EXIT
            END IF
            l=l+1
            IF(l == nzz1) THEN
              PRINT *,'htint:nzz1',nzz1
              DO kk=1,l
                PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
              END DO
              STOP 'htint'
            END IF
          END DO
        END DO
      END SUBROUTINE htint
      
      !========================================================================

      real function fu_satur_water_vapour_kPa(temperature)
        ! Description:
        !  Returns the saturation vapour pressure for a temperature
        !  valid for RANGE -40 C to +50 C
        !
        ! Method:
        ! The original Goff-Gratch-formulation ( in its  1946 modification). 
        !
        ! All units: NOT SI
        !
        ! Language: ANSI Fortran 90
        !
        ! Author: Ilkka Valkama, FMI, adapted by M.Sofiev, FMI

        IMPLICIT NONE

        ! Imported parameters with intent(in):
        REAL, INTENT(in) :: temperature
 
        ! Local declarations :

        REAL local_temperature1, local_temperature2, local_temperature3
        REAL water_vapour,local_vapour1,local_vapour2

        IF(temperature.lt.100)THEN        ! Only kelvin is acceptable
          local_temperature1=temperature+273.16
        ELSE
          local_temperature1=temperature
        END IF

        local_temperature2= local_temperature1/273.16
        local_temperature3= 1./local_temperature2

        local_vapour1= 8.29692*(1.0 - local_temperature2)                   
        local_vapour2= 4.76955*(1.0 - local_temperature3)                 
        local_vapour1= 1.0 - 10.**local_vapour1
        local_vapour2= 10.**local_vapour2 - 1.

        water_vapour = 10.79574*(1.0 - local_temperature3) - 5.0280 &
	                 & *(LOG10(local_temperature2))
        water_vapour=water_vapour + 1.50475E-4*local_vapour1 + 0.42873E-3 *local_vapour2

        fu_satur_water_vapour_kPa = 0.1 * 10.**(water_vapour + 0.78614)  ! in kPa
        !
        ! ANother alternative (difference in 3-rd digit): 
        ! https://www.cs.helsinki.fi/u/ssmoland/physics/envphys/lecture_2.pdf
        !
        !fu_satur_water_vapour_kPa = 0.001 * exp(77.34 - &
        !                                      & 7235/local_temperature1 - &
        !                                      & 8.2*log(local_temperature1) + &
        !                                      & 0.005711 * local_temperature1)
        
      end function fu_satur_water_vapour_kPa
      
    end subroutine interpolate_env_profile

end module plume_rise_driver
