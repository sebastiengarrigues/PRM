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
  use eccodes

  !$ use omp_lib
  
  implicit none
  
  private

  ! Public subrojtines from this module
  public compute_plume_rise

  public acquire_fires
  public acquire_meteo
  public store_plume_rise
  
  ! private subroutines of this module
  private fu_index
  
  ! Public constants from this module

  ! GRIB short names
  character(len=20), dimension(4), parameter, public :: quantities_3D = [Character(len=20) :: 'u','v','t','q']
  integer, parameter, public :: indU = 1
  integer, parameter, public :: indV = 2
  integer, parameter, public :: indT = 3
  integer, parameter, public :: indQ = 4
  integer, parameter, public :: nMetQ3D = 4 ! size(quantities_3D)

  character(len=20), dimension(2), parameter, public :: quantities_2D = [Character(len=20) :: 'blh', 'sp' ]
  integer, parameter, public :: indBLH = 1 ! 
  integer, parameter, public :: indSP = 2
  integer, parameter, public :: nMetQ2D = 2 !size(quantities_2D)

  
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
  
  type(Tgrid_lonlat), parameter :: grid_missing = Tgrid_lonlat('',-1,-1,-1,-1,-1,-1)

  public Tgrid_lonlat
  !
  ! Fire data structure
  !
  type Tfires
    integer :: nFires                     ! the number of fires in the set
    real, dimension(:), allocatable :: lon, lat, FRP, burnt_area, mdur, moist   ! (nFires)
    real, dimension(:,:), allocatable :: inj_IS4FIRES, inj_PRM  ! (2,nFires) injection bottom and top, [m]
    integer :: valid_date, valid_time
  end type Tfires
  public Tfires
  !
  ! Meteodata structure
  !
  type Tmeteo_data
    real, dimension(:,:,:), allocatable :: data3d     ! (nz, nfires, nMetQ3D)
    real, dimension(:,:), allocatable :: data2d       ! (nfires, nMetQ2D)
    integer :: nbr_of_levels  ! hybrid_sigma_pressure 
    real, dimension(:), allocatable :: a_interface,  &  !! [a] = [Pa], [b] = [1]
                                      & b_interface  ! hybrid coefficients of the upper interface of the layers
    integer :: valid_date, valid_time
  end type Tmeteo_data 
  public Tmeteo_data

CONTAINS
  
  !*************************************************************
  
  subroutine compute_plume_rise(fires, &                ! fires, input and output
                              & meteo_data, &           ! meteo data, input
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
    character(len=*), intent(in) :: chDumpTemplate

    integer :: iFire
    real :: BVfreq
    
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
    
    
    ! Scan all fires calling the corresponding plume-rise routines.
    !
    
    fires%inj_IS4FIRES(:,:) = -1
    do iFire = 1, fires%nFires
      !
      ! Plume rise from IS4FIRES.
      ! The model works in SI units
      !
      call bvfreq_for_fire(meteo_data,iFire, BVFreq)
      call IS4FIRES_vertical_profile(fires%FRP(iFire), &
                                   & meteo_data%data2d(iFire, indBLH), &
                                   & BVFreq, &
                                   & .true., &
                                   & fires%inj_IS4FIRES(1,iFire), fires%inj_IS4FIRES(2,iFire))
!      write(*,*)'>>>>>>>>>>>> IS4FIRES bottom=', fires%inj_IS4FIRES(1,iFire), &
!                                           & 'top=', fires%inj_IS4FIRES(2,iFire)
      !
      ! Plume rise from PRM
      ! PRM works in randomly picked units, watchout the conversion
      !
      fires%inj_PRM = -1
!!!!      call plumerise_PRMv1(meteo_data, iFire, &
!!!!                         & fires%burnt_area(iFire), &  ! from m2 
!!!!                         & fires%frp(iFire) * 1e-6, &         ! from W to MW
!!!!                         & fires%mdur(iFire) / 60., &         ! from sec to min
!!!!                         & fires%moist(iFire), &        ! remains fraction
!!!!                         & alpha, C_epsi, C_delta, C_wind, C_delta_wind,  &  ! fire configuration 
!!!!                         & FRP2CHF,  &
!!!!                         & wind_eff, micro_eff, &
!!!!                         & .True., .False.,   &   ! ifDebugPrint, ifStoreDump, 
!!!!                         & 6, & ! debug/dump print file
!!!!                         & fires%inj_PRM(1,iFire), fires%inj_PRM(2,iFire))
      write(*,*)'>>>>>>>>>>>> IS4FIRES again bottom=', fires%inj_IS4FIRES(1,iFire), &
                                                 & 'top=', fires%inj_IS4FIRES(2,iFire)
   !   write(*,*)'>>>>>>>>>>>>>>>>>>>>>>> PRM bottom=', fires%inj_PRM(1,iFire), &
   !                                              & 'top=', fires%inj_PRM(2,iFire)
    end do ! cycle over files
    
  end subroutine compute_plume_rise

 !*******************************************************

 subroutine bvfreq_for_fire(meteodata,iFire, BVFreq)
      implicit none
      type (Tmeteo_data), intent(in) :: meteodata
      integer, intent(in) :: iFire
      real, intent(out) :: BVFreq

      !FIXME
      ! Some smarter thing should go here
      BVFreq =  1e-3

 end subroutine bvfreq_for_fire



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

  

  !************************************************************************

  subroutine plumerise_PRMv1( meteo_data, iFire, &
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
    type(Tmeteo_data), intent(in) :: meteo_data
    integer,intent(in) :: ifire, uDump
    real,intent(in) :: burnt_area, frp_fire, mdur, moist ! burnt_area in M2
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
    call interpolate_env_profile(meteo_data, meteo_data%nbr_of_levels, iFire,  PRM_data)
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

  subroutine interpolate_env_profile(meteo_data, nz_meteo, iFire, dat)
    !
    ! Generates the vertical of PRM and
    ! Projects environmental profiles to it
    ! Calculating derived parameters
    !
    implicit none

    ! Imported parameters
    type (Tmeteo_data), intent(in) :: meteo_data
    integer, intent(in) :: nz_meteo, iFire
    type(TPRM_data), intent(inout) :: dat
    
    real, dimension(nz_meteo):: zMeteo, pMeteo
    ! Local variables
    integer :: k, iLevMet, iLevPRM
    real :: es, dz,p,t, pBot, pTop, ps, zbot, ztop
    real, parameter ::  RAir =   8.314 / 0.02897 ! Rgas / mu_air : [J/K/mol] /  [kg/mol] 
    real, parameter ::  g  = 9.8

    ! set the PRM grid
    call set_grid(dat) ! define vertical grid of plume model 
    !
    ! interpolation of the environmental parameters
    !

    dat%hBL = meteo_data%data2D(iFire, indBLH)
    ps = meteo_data%data2D(iFire, indSP)

    zBot = 0
    pBot =      meteo_data%a_interface(meteo_data%nbr_of_levels+1) &
           &  + meteo_data%b_interface(meteo_data%nbr_of_levels+1)  * ps

    iLevPRM = 1         

M:  do iLevMet = meteo_data%nbr_of_levels,1,-1
      pTop =  meteo_data%a_interface(iLevMet) +  meteo_data%b_interface(iLevMet)*ps
      T = meteo_data%data3D(iLevMet, iFire, indT)
      p = 0.5* (pBot + pTop)
      dz  =  (pBot - pTop) / p * T * RAir / g
      zTop = zBot + dz
      do while (dat%zt(iLevPRM) < zTop) 
        dat%upe(iLevPRM) = meteo_data%data3D(iLevMet, iFire, indU)
        dat%upe(iLevPRM) = meteo_data%data3D(iLevMet, iFire, indV)
        dat%rhe(iLevPRM) = meteo_data%data3D(iLevMet, iFire, indQ)
        dat%te(iLevPRM)  = T
        dat%pe(iLevPRM)  =  p
        dat%the(iLevPRM)  = T *  (1e5/p)**0.2854


        !  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
        ES = fu_satur_water_vapour_kPa (dat%TE(iLevPRM)) !blob saturation vapor pressure, kPa
        dat%QSAT (iLevPRM) = (.622 * ES) / (dat%PE(iLevPRM) * 1e-3 - ES) !saturation lwc g/g
        !calcula qvenv
        dat%qvenv(iLevPRM) = max(0.01 * dat%rhe(iLevPRM) * dat%QSAT(iLevPRM), 1e-8)
        !   print*,iLevPRM,dat%QSAT (iLevPRM),rhe(iLevPRM),qvenv(iLevPRM),rhe(iLevPRM)*dat%QSAT (iLevPRM)

        dat%thve(iLevPRM) = dat%the(iLevPRM) * (1. + .61*dat%qvenv(iLevPRM)) ! virtual pot temperature
        dat%dne(iLevPRM) = dat%pe(iLevPRM) *1e3 / (rgas * dat%te(iLevPRM)*(1. + .61*dat%qvenv(iLevPRM))) !  dry air density (kg/m3)
        dat%vel_e(iLevPRM) = sqrt(dat%upe(iLevPRM)**2 + dat%vpe(iLevPRM)**2)

        call thetae(dat%pe(iLevPRM), dat%te(iLevPRM), dat%qvenv(iLevPRM), dat%thee(iLevPRM), dat%tde(iLevPRM))

        !--------- converte press de Pa para kPa para uso modelo de plumerise
        dat%pe(iLevPRM) = dat%pe(iLevPRM) * 1.e-3

        dat%SCE(iLevPRM) = 0.e0 !++ rp intialize the scalar in the environment
        iLevPRM = iLevPRM + 1

        if (iLevPRM > nzPRM) then
              exit M
        endif
      enddo !! over iLevPRM
      zBot = zTop
      pBot = pTop
    enddo M
    if (iLevMet < 1) then
      print *, "iLevMet = ", iLevMet, ' zBot, zTop, p: ', zBot, zTop, p, "iLevPRM, nzPRM, dat%zt(iLevPRM) ", iLevPRM, &
                 & nzPRM, dat%zt(iLevPRM)
      print *, "This should not happen, iLevMet =", iLevMet 
      stop 1
    endif




    CONTAINS

      !===========================================================

      subroutine set_grid(dat)

        implicit none
        type(TPRM_data), intent(inout) :: dat
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


     !*************************************************************


    subroutine grid_from_grib(igrib_in, grid)
        integer, intent(in) :: igrib_in
        type( Tgrid_lonlat), intent(out) :: grid

        character(len=256) :: strTmp
        real*8 :: dTmp
        integer :: iTmp, jTmp, iStat
        character(len=*), parameter :: sub_name = 'grid_from_grib'


        call grib_get(igrib_in, 'gridType', strTmp, status = iStat)
        call codes_check(iStat, sub_name, 'gridType')
        if (strTmp /= 'regular_ll') then
          print *, 'Grid ', strTmp, 'not imlemented, only regular_ll so far...'
          stop 1
        endif 
        grid%grid_type = 'lonlat'

        call codes_get(igrib_in,'Ni', grid%nlon, status = iStat)
        call codes_check(iStat, sub_name, 'Ni')
        call codes_get(igrib_in,'Nj', grid%nlat, status = iStat)
        call codes_check(iStat, sub_name, 'Nj')

        call codes_get(igrib_in,'Ni', grid%nlon, status = iStat)
        call codes_check(iStat, sub_name, 'Ni')
        call codes_get(igrib_in,'Nj', grid%nlat, status = iStat)
        call codes_check(iStat, sub_name, 'Nj')



        !!! Here we do not flip the input array, but rather flip the origin and the sign of the increment
        !! Longitude
        call codes_get(igrib_in,'iDirectionIncrementInDegrees', grid%dlon_deg, status = iStat)
        call codes_check(iStat, sub_name, 'iDirectionIncrementInDegrees')
        call grib_get(igrib_in, 'iScansNegatively', iTmp, status = iStat)
        call codes_check(iStat, sub_name, 'iScansNegatively')

        call codes_get(igrib_in,'longitudeOfFirstGridPointInDegrees', grid%sw_corner_lon, status = iStat)
        call codes_check(iStat, sub_name, 'longitudeOfFirstGridPointInDegrees')
        if (iTmp == 0) then !!Negative scan order 
          grid%dlat_deg = -1 * grid%dlat_deg
        endif

        !! Latitude
        call codes_get(igrib_in,'jDirectionIncrementInDegrees', grid%dlat_deg, status = iStat)
        call codes_check(iStat, sub_name, 'jDirectionIncrementInDegrees')
        call grib_get(igrib_in, 'jScansPositively', iTmp, status = iStat)
        call codes_check(iStat, sub_name, 'jScanspositively')

        call codes_get(igrib_in,'latitudeOfFirstGridPointInDegrees', grid%sw_corner_lat, status = iStat)
        call codes_check(iStat, sub_name, 'latitudeOfFirstGridPointInDegrees')
        if (iTmp == 0) then !!Negative scan order
          grid%dlat_deg = -1 * grid%dlat_deg
        endif


    
    end subroutine grid_from_grib

     !*************************************************************


     subroutine datetime_from_grib(igrib_in, dataDate, dataTime)
        implicit none
        integer, intent(in) :: igrib_in
        integer, intent(out) :: dataDate, dataTime
        character(len=*), parameter :: sub_name = 'datetime_from_grib'
        integer :: iStat

        call codes_get(igrib_in,'dataDate', dataDate, status=iStat)
        call codes_check(iStat, sub_name, 'dataDate')
        call codes_get(igrib_in,'dataTime', dataTime, status=iStat)
        call codes_check(iStat, sub_name, 'dataTime')


     end subroutine datetime_from_grib

     !*************************************************************

    

     subroutine acquire_fires(fires, fire_grib)

       implicit none
        type(Tfires), intent(out) :: fires
        character(len=*), intent(in) :: fire_grib

        integer :: iUnit, igrib_in, nFires, iTmp, jTmp, iStat
        real :: fTmp
        character(len=200) :: shortname
        real, dimension(:), allocatable   :: values
        type( Tgrid_lonlat) :: grid
        character(len = *), parameter :: frpname = 'frpfire'  ! ParamID 210099
        character(len=*), parameter :: sub_name = 'acquire_fires'

        call codes_open_file(iUnit,fire_grib,'r')

        igrib_in = 1
        do while (igrib_in > 0) 
          call codes_grib_new_from_file(iUnit,igrib_in)
          call codes_get(igrib_in,'shortName',shortname, status=iStat)
          call codes_check(iStat, sub_name, 'shortName')
          if (shortname == frpname) exit
          call codes_release(igrib_in)
        enddo

        if (igrib_in < 0) then
          print *, "Failed to find frpfire in ", fire_grib 
          stop
        endif
        call grid_from_grib(igrib_in, grid)
        print *, "Fire grid:", grid
        allocate(values(grid%nlon * grid%nlat))
        call codes_get(igrib_in,'values', values, status=iStat) !! Comes as W/m2
        call codes_check(iStat, sub_name, 'values')

        call datetime_from_grib(igrib_in, fires%valid_date, fires%valid_time)

        ! Init fires
        fires%nFires= count(values>0)
        allocate(fires%lon(fires%nFires), fires%lat(fires%nFires), fires%FRP(fires%nFires), &
               & fires%burnt_area(fires%nFires), fires%mdur(fires%nFires), fires%moist(fires%nFires))
        fires%burnt_area(:) = F_NAN
        fires%mdur(:)  = F_NAN
        fires%moist(:) =  F_NAN

        allocate(fires%inj_IS4FIRES(2, fires%nFires), fires%inj_PRM(2, fires%nFires))
        fires%inj_IS4FIRES(:,:) = F_NAN
        fires%inj_PRM(:,:) = F_NAN

        ! Put fires to the structure
        jTmp = 1
        fTmp = earth_radius * degrees_to_radians !! degrees to meters
        fTmp =  fTmp * fTmp !! square of it
        do iTmp = 1, Size(values)
          if (values(iTmp) <= 0.) cycle
          
          fires%lon(jTmp) = grid%sw_corner_lon + modulo(iTmp-1,grid%nlon) * grid%dlon_deg
          fires%lat(jTmp) = grid%sw_corner_lat + ((iTmp - 1)  / grid%nlon)  * grid%dlat_deg
          if (fires%lat(jTmp) < -90) then
            print *, 'Ops'
          endif

          ! incoming frp * cell_area
          fires%FRP(jTmp) = values(iTmp) * abs(grid%dlon_deg*grid%dlat_deg) * fTmp &
                        & * cos(fires%lat(jTmp)*degrees_to_radians)

          jTmp = jTmp + 1
        enddo

        !FIXME   Put something to make Plume-rise model happpy
        fires%burnt_area(:) = 1e6     ! m2  MODIS pixel is 1x1 km2 or more
        fires%mdur(:) = 3600.0              ! fire duration, sec
        fires%moist(:) = 0.8                ! fuel moisture, fraction


        if (jTmp - 1 /= fires%nFires) then
          print *, 'Oooops: Fire number mismatch: ', jTmp -1,  fires%nFires
          stop 3
        endif
        call codes_release(igrib_in)
        call codes_close_file(iUnit) 

    end subroutine acquire_fires

    
  !***************************************************************
  
  subroutine fires_to_meteo(fires,  meteo_grid, idxMeteo)
    !
    ! Assigns meteo index to fires and sets ifInsideGrid(iFire)

    implicit none
    
    ! imported parameters
    type(Tfires), intent(in) :: fires
    type(Tgrid_lonlat), intent(in) :: meteo_grid
    integer, dimension(:), allocatable, intent(inout) ::  idxMeteo ! Can be allocated or not
    
    ! Local variables
    integer :: iMet, iFire, ix, iy, izMet
    real :: xNew, yNew


    if (.not. allocated(idxMeteo)) then
        allocate(idxMeteo(fires%nFires))
    endif
    !
    ! Scan the fires one by one and assign coordinates
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
        idxMeteo(iFire) = -1
      else
        idxMeteo(iFire) = ix + (iy - 1) * meteo_grid%nlon
      endif
    end do  ! nFires
    print *, "Projected ", count(idxMeteo>0), ' fires of total ', fires%nFires

  end subroutine fires_to_meteo


  !************************************************************************


    subroutine acquire_meteo(meteo, fires, meteo_grib_files)

      !
      !Store meteo profiles for each fire. 
      ! Extraction of fire columns from GRIB files
      ! 
     
 
      implicit none
      type(Tmeteo_data), intent(out) :: meteo
      character(len=*), dimension(:), intent(in) :: meteo_grib_files
      type(Tfires), intent(in) :: fires

      integer :: iUnit, igrib_in, nFires, iTmp, jTmp, iStat, iLev, iVal2D, iVal3D, iFile
      integer :: dataDate, dataTime
      real :: fTmp
      character(len=200) :: typeOfLevel, shortname
      real, dimension(:), allocatable   :: values, pv
      integer, dimension(:), allocatable   :: idxMeteo
      type( Tgrid_lonlat) :: grid, meteo_grid
      character(len =* ), parameter :: frpname = 'frpfire'
      character(len=*), parameter :: sub_name = 'acquire_meteo'


      
      allocate( meteo%data2d(fires%nFires, nMetQ2D))
      meteo%data2d(:, :) = F_NAN

      meteo_grid = grid_missing
      meteo%valid_date = fires%valid_date
      meteo%valid_time = fires%valid_time
     

      do iFile = 1, size(meteo_grib_files)

        call codes_open_file(iUnit,meteo_grib_files(iFile),'r', status=iStat)
        call codes_check(iStat, sub_name, 'open')

        do while (.True.) 
            call codes_grib_new_from_file(iUnit, igrib_in)
            if (igrib_in < 1) exit

            call datetime_from_grib(igrib_in, dataDate, dataTime)
            call codes_get(igrib_in,'typeOfLevel', typeOfLevel, status=iStat)
            call codes_check(iStat, sub_name, 'typeOfLevel')
            call grid_from_grib(igrib_in, grid)

            !!if ((dataDate /= meteo%valid_date) .or. (dataTime /= meteo%valid_time)) then
            if ((dataDate /= meteo%valid_date) ) then  !!! Proper FIXME time matching needed
                call codes_release(igrib_in)

                print *, "skipping message ", dataDate, dataTime
                print *, "Fires date and time", meteo%valid_date, meteo%valid_time
                cycle
            endif
            
            if (  meteo_grid%nlat /= grid%nlat .or. meteo_grid%nlon /= grid%nlon ) then ! Poor-man's "/="
              meteo_grid = grid
              print *, "Meteo grid:", grid
              call  fires_to_meteo(fires,  meteo_grid, idxMeteo)
              if (.not. allocated(values)) allocate(values(meteo_grid%nlat * meteo_grid%nlon))
            endif

            call codes_get(igrib_in,'shortName',shortName, status=iStat)
            call codes_check(iStat, sub_name, 'shortName')

            ! Not found yet
            iVal3D = -1
            iVal2D = -1

            if (typeOfLevel == 'hybrid' ) then
                 ! Do we ned it at all?
                 do iVal3D = 1, nMetQ3D
                    if (quantities_3d(iVal3D) == shortName) then
                      call codes_get(igrib_in,'level', iLev, status=iStat)
                      call codes_check(iStat, sub_name, 'level')
                      exit
                   endif
                 enddo

                 if (ival3D > nmetQ3D) then
                    ival3D = -1
                    call codes_release(igrib_in)
                !    print *, "skipping message 3D"
                    cycle 
                 endif

                 ! 
                 ! iquantity and level is here now, make sure that the vertical is okay

                 call codes_get(igrib_in,'NV', iTmp, status=iStat)
                 call codes_check(iStat, sub_name, 'NV')
                 if (.not. allocated(meteo%a_interface)) then
                      meteo%nbr_of_levels = iTmp/2 - 1 
                      allocate(pv(iTmp), meteo%a_interface(iTmp/2), meteo%b_interface(iTmp/2))
                      allocate(meteo%data3d(meteo%nbr_of_levels,fires%nFires, nMetQ3D))
                      meteo%a_interface(:) = F_NAN
                      meteo%b_interface(:) = F_NAN
                      meteo%data3d(:,:,:) = F_NAN
                 endif
                 if (meteo%nbr_of_levels + 1 /= iTmp/2) then
                    print *, "Vertical size mismatch"
                    stop 1
                 endif
                 call codes_get(igrib_in,'pv', pv, status=iStat)
                 call codes_check(iStat, sub_name, 'pv')
                 if (meteo%a_interface(1) /= meteo%a_interface(1)) then ! Verical not yet defined
                    meteo%a_interface(:) =  pv(1:iTmp/2)  !!Pa
                    meteo%b_interface(:) =  pv(iTmp/2+1:iTmp) ! 
                 else
                     if (any(meteo%a_interface /= pv(1:iTmp/2))) then
                        print *, "Vertical values mismatch"
                        stop 1
                     endif
                 endif

            elseif (typeOfLevel == 'surface' ) then
                 do iVal2D = 1, nMetQ2D
                    if (quantities_2D(iVal2D) == shortName) exit
                 enddo

                 if (ival2D > nMetQ2D) then
                    ival2D = -1
                    call codes_release(igrib_in)
               !     print *, "skipping message 2D "
                    cycle 
                 endif
            else
                call codes_release(igrib_in)
                print *, "skipping message (leveltype)"
                cycle 
            endif

            !
            ! Now we have either iVal2D > 0 or (iVal3D > 0 and reasonable iLevel)

            call codes_get(igrib_in,'values',values, status=iStat)
            call codes_check(iStat, sub_name, 'values')
            if (iVal2D > 0 ) then
              print *, "Store 2D ", trim(shortName), ', lev', iLev
              do iTmp = 1,fires%nFires
                  if (idxMeteo(iTmp)>0) then
                    meteo%data2d(iTmp, iVal2D) = values(idxMeteo(iTmp))
                  endif
              enddo
            elseif (iVal3D > 0 ) then
              !print *, "Store 3D ", trim(shortName), ', lev', iLev
              do iTmp = 1,fires%nFires
                  if (idxMeteo(iTmp)>0) then
                    meteo%data3d(iLev,iTmp, iVal3D) = values(idxMeteo(iTmp))
                  endif
              enddo
            else
              print *, "This should not happen ever!", sub_name
              stop 2
           endif
           call codes_release(igrib_in)
        enddo
        call codes_close_file(iUnit)
      enddo

    end subroutine acquire_meteo

  !************************************************************************

     subroutine store_plume_rise(fires, frpgrib, plumerisegrib) 

       !
       ! Stores output by substituting parmeterID and alues in the input messages
       !

       implicit none
        type(Tfires), intent(in) :: fires
        character(len=*), intent(in) :: frpgrib, plumerisegrib

        integer :: iUnit, indGrib,  indGribOut, iCell, iFire, iStat
        real :: fTmp
        character(len=200) :: shortname
        real, dimension(:), allocatable   :: values
        integer, dimension(:), allocatable   :: iOutFire ! Index in the output grid for plume
        type( Tgrid_lonlat) :: grid
        character(len = *), parameter :: frpname = 'frpfire'
        character(len=*), parameter :: sub_name = 'store_plume_rise'
        integer, parameter :: ParamID_I4Finjh = 210060  !	Injection height (from IS4FIRES)  210060 
        integer, parameter :: ParamID_I4Ftop = 210119 ! Mean altitude of maximum injection 210119 !! Metres above sea level

        call codes_open_file(iUnit, frpgrib ,'r')

        indGrib = 1
        do while (indGrib > 0) 
          call codes_grib_new_from_file(iUnit,indGrib)
          call codes_get(indGrib,'shortName',shortname, status=iStat)
          call codes_check(iStat, sub_name, 'shortName')
          if (shortname == frpname) exit
          call codes_release(indGrib)
        enddo

        if (indGrib < 0) then
          print *, "Failed to find frpfire in ", frpgrib
          stop
        endif
        call grid_from_grib(indGrib, grid)
        print *, "Fire grid:", grid
        allocate(values(grid%nlon * grid%nlat), iOutFire(fires%nfires))
        call codes_get(indGrib,'values', values, status=iStat) !! Comes as W/m2
        call codes_check(iStat, sub_name, 'values')
        call codes_close_file(iUnit) 


        !Here we rely on the fact that the frp file is exactly the same as one used for 
        iFire = 0
        do iCell = 1, Size(values)
          if (values(iCell) <= 0.) cycle
          iFire = iFire + 1
          if (iFire > fires%nfires) then
            print *, "Fire index ", iFire, " > fires%nfires = ", fires%nfires
            stop 1
          endif
          iOutFire(iFire) = iCell
        enddo
        if (iFire /= fires%nFires) then
          print *, 'Oooops: Fire number mismatch in ', sub_name, " : ", iFire,  fires%nFires
          stop 3
        endif


        ! Prepare array of "Injection height (from IS4FIRES)" from "plume bottom" and "plume top"
        ! with simple average
        do iFire = 1,fires%nfires
            iCell = iOutFire(iFire)
            fTmp = values(iCell)
            values(iCell) = 0.5*(fires%inj_IS4FIRES(1,iFire) + fires%inj_IS4FIRES(2,iFire))
        enddo 

        call codes_open_file(iUnit, plumerisegrib, 'w')
        
        call codes_set(indGrib, 'paramId' ,ParamID_I4Finjh, status=iStat)
        call codes_check(iStat, sub_name, 'Set paramId')

        call codes_set(indGrib, 'packingType' , 'grid_simple', status=iStat)
        call codes_check(iStat, sub_name, 'set packingType')

        call codes_set(indGrib, 'values' , values, status=iStat)
        call codes_check(iStat, sub_name, 'Set values')

        call codes_write(indGrib, iUnit, status=iStat)
        call codes_check(iStat, sub_name, 'codes_write')
        call codes_release(indGrib)
        call codes_close_file(iUnit) 
        deallocate (values, iOutFire)

    end  subroutine store_plume_rise

    
  !***************************************************************
  

end module plume_rise_driver
