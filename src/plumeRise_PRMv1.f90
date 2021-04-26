MODULE plume_rise_PRMv1
  !
  ! This is the module for the plume-rise model PRM v.1.
  ! Adopted with minor changes from the GFAS code.
  ! 
  ! Adaptatins: M.Sofiev, 2020
  ! Language: FORTRAN-90, well, close to that
  !
  use rconstants
!  use plumegen_coms
!  use var_out   !, only : profile_mass_detrained_int_out, profile_mass_entrained_int_out
  
  IMPLICIT NONE
  
  private     ! close all module variables and functions

  public MAKEPLUME     ! actual plume rise computations, the only public item here
  public INITIAL       ! initialization of the PRMv1 instance
  public pack_fire_properties  ! storing the fire features in the data structure
  public thetae

  ! Private subroutines used in this module
  private tend0_plumerise ! tendency for the plume-rise model
  private vel_advectc_plumerise  ! advection contribution to W tendency
  private visc_W
  private update_plumerise  ! actual timestep: new state from the previous one and tendencies
  private predict_plumerise ! something similar as well
  private buoyancy_plumerise
  private ENTRAINMENT       ! exchange with the outside environment
  private scl_misc
  private damp_grav_wave !friction
  private fallpart
  private entrainment_coeff
  private esat_hPa
  private LBOUND_MTT
  private BURN
  private WATERBAL
  private ESAT_PR
  private printout

  ! Module parameters
  !
  ! THe main dimension of the model
  ! - vertical dimension
  ! - temporal dimension
  !
  integer, parameter, public :: nzPRM = 200, ntimePRM = 4000   ! vertical and temporal dimension
  real, parameter, public :: dzPRM = 100.  ! thickness of PRM vertical layers

  !====================================================================================
  !
  ! All variables needed for the plume rise are encapsulated in the derived type
  !
  type TPRM_data
!    private
    real, dimension(nzPRM) :: w=0, t=0, theq=0, qv=0, qc=0, qh=0, qi=0, sc=0,  &  ! blob
                          & vth=0, vti=0, rho=0, txs=0, qt=0,                &
                          & est=0, qsat=0, qpas=0, qtotal=0, td=0, vel_p=0, rad_p=0, rho_e=0, thv=0
    real, dimension(nzPRM) :: wc=0, wt=0, tt=0, qvt=0, qct=0, qht=0, qit=0, sct=0, vel_t=0, rad_t=0
    real, dimension(nzPRM) :: dzm=0, dzt=0, zm=0, zt=0,      &
                          & vt3dc=0, vt3df=0, vt3dk=0, vt3dg=0, scr1=0,     &
                          & vt3dj=0, vt3dn=0, vt3do=0, vt3da=0, scr2=0,     &
                          & vt3db=0, vt3dd=0, vt3dl=0, vt3dm=0, vt3di=0,    &
                          & vt3dh=0, vt3de=0, rbuoy=0, dwdt_entr_dyn=0
    real, dimension(nzPRM) :: rad_pC=0, profile_mass_detrained=0, profile_mass_entrained=0, mass_sc=0
    real, dimension(3,nzPRM) :: hs_rad_p_noDe=0
    real, dimension(4,nzPRM) :: dwdt_entr=0
    ! environment at plume grid
    real, dimension(nzPRM) :: pke=0, the=0, thve=0, thee=0, pe=0, te=0, qvenv=0, &
                          & rhe=0, dne=0, sce=0, tde=0, upe=0, vpe=0, vel_e=0, VISC=0, CVH=0, CVI=0
    real :: ztop_max_global=0
    integer :: entrainment_model=0
    integer :: MTT=0
    real :: DZ=0, VISCOSITY=0, TSTPF=0
    integer :: N=0
    real :: ADVW=0, ADVT=0, ADVV=0, ADVC=0, ADVH=0, ADVI=0, ADIABAT=0, WBAR=0, VHREL=0, VIREL=0  ! advection
    real :: ZSURF=0, ZBASE=0, ZTOP=0
    integer :: LBASE=0
    ! entrainment
    real :: AREA=0, RSURF=0, ALPHA=0, ALPHA_dyn=0
    real, dimension(nzPRM) :: RADIUS=0, radius_mdt=0
    !
    real, dimension(ntimePRM) :: HEATING=0
    real :: FMOIST=0
    !
    real :: DT=0, dt_mdt=0
    integer :: MINTIME=0, MDUR=0, MAXTIME=0

    logical :: wind_eff=.false., micro_eff=.false.
    ! ++ rp height of the mixing layer
    real :: hBL=0
    real :: W_Treshold=0
    real :: C_epsi=0, C_delta=0, C_wind=0, C_delta_wind=0
    real :: FRP2TOTALH=0, CHF2TOTALH=0
  end type TPRM_data
  public TPRM_data
  
  
  CONTAINS
    
    !
    SUBROUTINE MAKEPLUME (dat, heat_fluxW, & ! input
                        & ifDebugPrint, ifStoreDump, uDump, & ! debug/dump print & files where to put it
                        & ztopmax, profile_mass_detrained_int, profile_mass_entrained_int)  ! output
      !
      ! *********************************************************************
      !
      !    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
      !    Alan Weinstein, JAS V.27 pp 246-255. (W),
      !    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
      !    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
      !    Originally developed by: Don Latham (USFS)
      !
      !
      ! ************************ VARIABLE ID ********************************
      !
      !     DT=COMPUTING TIME INCREMENT (SEC)
      !     DZ=VERTICAL INCREMENT (M)
      !     LBASE=LEVEL ,CLOUD BASE
      !
      !     CONSTANTS:
      !       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
      !       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
      !       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
      !       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
      !       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
      !       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
      !       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
      !       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
      !       TFREEZE = FREEZING TEMPERATURE (K)
      !
      !
      !     PARCEL VALUES:
      !       T = TEMPERATURE (K)
      !       TXS = TEMPERATURE EXCESS (K)
      !       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
      !       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
      !       QC = WATER CONTENT (G/G DRY AIR)
      !       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
      !       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
      !       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
      !       ES = SATURATION VAPOR PRESSURE (kPa)
      !
      !     ENVIRONMENT VALUES:
      !       TE = TEMPERATURE (K)
      !       PE = PRESSURE (kPa)
      !       QVENV = WATER VAPOR (G/G)
      !       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
      !       DNE = dry air density (kg/m^3)
      !
      !     HEAT VALUES:
      !       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
      !       MDUR = DURATION OF BURN, MINUTES
      !
      !       W = VERTICAL VELOCITY (M/S)
      !       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
      !	RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
      !	ALPHA = ENTRAINMENT CONSTANT
      !       MAXTIME = TERMINATION TIME (MIN)
      !
      !     Dump and output files
      !       ifDebugPrint - controls test messages
      !       ifStoreDump - controls additional output from intermediate variables and time steps
      !       All dump files are now in a single file uDump, all contrlled by ifStoreDump
      !       uDump - the main dump file, many variables and time steps
    
      !if (ifStoreDump) then
      !  open(unit=1000,file='massConservation.txt')
      !  open(unit=1004,file='massConservation_sc.txt')
      !  open(unit=1002,file='eflux.txt')
      !endif


      !
      !
      !**********************************************************************
      !**********************************************************************               

      implicit none 

      type(TPRM_data), intent(inout) :: dat
      integer,intent(in) :: uDump
      real,intent(in) :: heat_fluxW
      logical, intent(in) :: ifDebugPrint, ifStoreDump  ! massive print to screen; writing detailed dump to file
      real,intent(out) :: ztopmax
      real,dimension(:),intent(inout) :: profile_mass_detrained_int, profile_mass_entrained_int

      !logical :: endspace  
      character (len=10) :: varn
      integer ::  iconv,  itime, k, kk, kkmax, deltak, nrectotal,i_micro,n_sub_step
      real ::  vc, TIME=0,  DQSDZ, &
               tmelt,  heatsubl,  heatfus,  heatcond, tfreeze, wmax, es, heat,dt_save
      !
      integer :: ii, L
      real :: mass_ejected,mass_in_plume,mass_ejected_int
      real :: mass_detrained,mass_detrained_int
      real :: mass_entrained,mass_entrained_int
      ! for the scalar
      real :: mass_in_plume_sc,mass_ejected_int_sc
      real :: mass_detrained_sc,mass_detrained_int_sc
      real :: mass_entrained_sc,mass_entrained_int_sc

      real :: v_ground,minsv,maxsv
      logical :: stationary
      real :: kappa,Rg,var_mass_in_plume,mass_in_plume_mdt,d_mass_in_plume,mintime_sav
      integer, parameter :: N_sav = 10 
      real,dimension(N_sav) :: var_mass_in_plume_ar, d_mass_in_plume_ar
      integer :: iii, NM1
      !
      ! initialization of a few internal variables
      !
      mintime_sav = 0.
      var_mass_in_plume_ar = 0.
      d_mass_in_plume_ar = 0.
      dat%tstpf = 2.0  	!- timestep factor

      nrectotal = nzPRM
      !
      !*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
      dat%mintime = 1  
      ztopmax = 0. 
      dat%ztop    = 0. 
      time = 0.  
      dat%dt = 1.
      wmax = 1. 
      kkmax   = 10
      deltaK  = 20
      dat%viscosity = 1500. !- viscosity constant (original value: 0.001)

      print *,'heat_fluxW', heat_fluxW

      v_ground = dat%w(1)
      mass_ejected_int = 0 
      mass_ejected_int_sc = 0 
      mass_detrained_int = 0.e0 
      mass_entrained_int = 0.e0
      mass_detrained_int_sc = 0.e0 
      mass_entrained_int_sc = 0.e0

      NM1 = nzPRM ! max(nzPRM, kkmax + deltak)

      dat%W_treshold = 1.e0  ! when w id less than this value the entrainment rate epsi is
                             ! switch to environement value.
      itime = 1
      stationary = .FALSE.

      ! ******************* model evolution ******************************
      !
      mass_in_plume_mdt = 0.
      mass_in_plume     = 0. 

      DO WHILE (TIME <= dat%MAXTIME)  !beginning of time loop

          !-- set model top integration
          !++ rp remove 
          !nm1 = min(nzPRM, kkmax + deltak)

          !-- set timestep
          !dt = (zm(2)-zm(1)) / (tstpf * wmax)  
          dat%dt = min(1.,(dat%zm(2)-dat%zm(1)) / (dat%tstpf * wmax))
                                
          !-- elapsed time, sec
          time = time + dat%dt

          !-- elapsed time, minutes                                      
          dat%mintime = 1 + int (time) / 60     
          wmax = 1.  !no zeroes allowed.


          !****************** BEGIN SPACE LOOP ******************

          !-- zerout all model tendencies
          if (ifDebugPrint) write(uDump,*) 'tend0_plumerise', dat%wc(1:5)
          call tend0_plumerise(dat%wt, dat%tt, dat%qvt, dat%qct, dat%qht, &
                             & dat%qit, dat%vel_t, dat%rad_t, dat%sct)

          !-- surface bounday conditions (k=1)
          if (ifDebugPrint) write(uDump,*) 'Lbound', dat%wc(1:5)

          call lbound_mtt(heat_fluxW, time, dat, ifDebugPrint, ifStoreDump, uDump)
          !call lbound_rio(heat_fluxW,time)

          if (v_ground .lt. dat%w(1) ) v_ground = dat%w(1)

          !-- dynamics for the level k>1 

          !-- W advection 
          if (ifDebugPrint) write(uDump,*) 'W advection', dat%wc(1:5)
      !   call vel_advectc_plumerise(NM1,WC,WT,DNE,DZM)
          call vel_advectc_plumerise(NM1, dat%WC, dat%WT, dat%RHO, dat%DZM)
  
          !-- scalars advection 1
          if (ifDebugPrint) write(uDump,*) 'scl advection', dat%wc(1:5)
          call scl_advectc_plumerise('SC',NM1,ztopmax, &
                                   & dat%dt, dat%T, dat%QV, dat%QC, dat%QI, dat%QH, dat%VEL_P, dat%SC, &
                                   & dat%TT, dat%QVT, dat%QCT, dat%QIT, dat%QHT, dat%VEL_T, dat%SCT, &
                                   & dat%w, dat%wc, dat%rho, dat%dzm, dat%dzt, dat%zt, dat%zm, &
                                   & dat%vt3dc, dat%vt3df, dat%vt3dg, dat%vt3dk)
    
          !-- scalars advection 2
          !call scl_advectc_plumerise2('SC',NM1)

          !-- Buoyancy (first call for the entrainment coefficient calculation)
          if (ifDebugPrint) write(uDump,*) 'buoyancy', dat%wc(1:5)
          call buoyancy_plumerise(NM1, dat%T, dat%TE, dat%QV, dat%QVENV, dat%QH, &
                                & dat%QI, dat%QC, dat%WT, dat%rbuoy)

          !-- scalars entrainment, adiabatic
          if (ifDebugPrint) write(uDump,*) 'scl entrainment', dat%wc(1:5)
          
          call scl_misc(NM1, dat%entrainment_model, dat%wind_eff, dat%W_treshold, dat%hBl, &
                      & dat%T, dat%TE, dat%QV, dat%QC, dat%QH, dat%QI,  dat%QVENV, dat%SC, dat%VEL_P, dat%VEL_E, &
                      & dat%TT, dat%QVT, dat%QCT, dat%QHT, dat%QIT, dat%SCT, dat%VEL_T, &
                      & dat%radius, dat%radius_mdt, dat%zm, dat%zt, dat%w, dat%rho, dat%rho_e, dat%dt_mdt, dat%rbuoy, &
                      & dat%C_epsi, dat%C_delta, dat%C_wind, dat%C_delta_wind, dat%WBAR, dat%dwdt_entr)

          !-- scalars dynamic entrainment
          !if (ifDebugPrint) print*, 'scl dyn'
          !if(wind_eff) call  scl_dyn_entrain(NM1)
    
          !-- gravity wave damping using Rayleigh friction layer fot T
          ! Direct call. no need for intermediate
          if (ifDebugPrint) write(uDump,*) 'damp grav', dat%wc(1:5)
!          call friction(1, nm1, deltak, dat%dt, dat%zt, dat%zm, dat%T, dat%TT, dat%TE) 
          
          call damp_grav_wave(1, nm1, deltak, dat%dt, dat%zt, dat%zm, dat%w, dat%t, dat%tt, &
                            & dat%qv, dat%qh, dat%qi, dat%qc, dat%te, dat%pe, dat%qvenv, &
                            & dat%vel_p, dat%vel_t, dat%vel_e, dat%rad_p, dat%rad_t)
          ! and for radius
          !call damp_grav_wave(3,nm1,deltak,dt,zt,zm,w,rad_p,rad_t,qv,qh,qi,qc,te,pe,qvenv&
          !                   ,vel_p,vel_t,vel_e,rad_p,rad_t)


          if (dat%micro_eff) then ! bypass microphysics
              !-- microphysics
              if (ifDebugPrint) write(uDump,*) ' micro'

              dt_save = dat%dt    ! going substeps in below calls, have to save this one
              n_sub_step=3
              dat%dt = dat%dt / float(n_sub_step)

              do i_micro=1,n_sub_step
                  !-- sedim ?
                  if (ifDebugPrint) write(uDump,*) ' micro: sedim' 
                  call fallpart(NM1, dat%RHO, dat%W, dat%CVH, dat%CVI, dat%VTH, dat%VHREL, &
                                   & dat%VTI, dat%VIREL, dat%QH, dat%QI, dat%ZM, dat%QHT, dat%QIT)
                  !-- microphysics
                  do L=2,nm1-1 
                      dat%WBAR    = 0.5*(dat%W(L)+dat%W(L-1))
                      es      = 0.1*esat_hPa (dat%t(L))            !blob saturation vapor pressure, em kPa
                      dat%qsat(L) = (eps * es) / (dat%pe(L) - es)  !blob saturation lwc g/g dry air
                      dat%EST (L) = ES  
                      dat%RHO (L) = 3483.8 * dat%PE(L) / dat%T(L)    ! AIR PARCEL DENSITY , G/M**3
                      !srf18jun2005
                      !	IF (W(L) .ge. 0.) DQSDZ = (QSAT(L  ) - QSAT(L-1)) / (ZT(L  ) -ZT(L-1))
                      !	IF (W(L) .lt. 0.) DQSDZ = (QSAT(L+1) - QSAT(L  )) / (ZT(L+1) -ZT(L  ))
                      DQSDZ = (dat%QSAT(L+1) - dat%QSAT(L-1)) / (dat%ZT(L+1 )-dat%ZT(L-1))

                      if (ifDebugPrint) print*, ' micro: waterbal', L , dat%w(1)         
                      
                      call waterbal(dat%T(L), dat%RHO(L), dat%QC(L), dat%QH(L), dat%QI(L), &
                                  & dat%QV(L), dat%QSAT(L), dat%WBAR, dat%DT, DQSDZ, dat%EST(L), &
                                  & dat%CVI(L))
                  enddo
              enddo
              dat%dt = dt_save
          endif
    
          !-- W-viscosity for stability 
          if (ifDebugPrint) print*, 'viscosity'
          call visc_W(nm1, deltak, ztopmax, &
                    & dat%zt, dat%zm, dat%W, dat%T, dat%QC, dat%QV, dat%QH, dat%QI, dat%SC, dat%VEL_P, &
                                    & dat%WT, dat%TT, dat%QCT, dat%QVT, dat%QHT, dat%QIT, dat%SCT, dat%VEL_T,  &
                    & dat%RAD_P, dat%visc)

          !-- update scalars
          if (ifDebugPrint) print*, 'update scalar'
          call update_plumerise(nm1,'S', dat%DT, &
                              & dat%WT, dat%TT, dat%QVT, dat%QCT, dat%QHT, dat%QIT, dat%VEL_T, dat%SCT, &
                              & dat%W, dat%T, dat%QV, dat%QC, dat%QH, dat%QI, dat%VEL_P, dat%SC)
             !print*,'wi apos update=',w(1:nm1)
             !print*,'Ti apos update=',T(1:nm1)


          ! hadvance_plumerise is meaningless, it just calls predict_plume_rise. The order of parameters:
          ! call hadvance_plumerise        (1, nm1, dat%dt, dat%WC, dat%WT, dat%W, dat%mintime)          
          !  subroutine hadvance_plumerise( iac, m1,   dt,     wc,     wt,     wp,  mintime)
          !call        predict_plumerise( m1,  wc,    wp, wt, dummy, iac, 2.*dt, eps)           
          !
          !call hadvance_plumerise(1, nm1, dat%dt, dat%WC, dat%WT, dat%W, dat%mintime)
          call predict_plumerise(nm1, dat%wc, dat%w, dat%wt, 1, &
                               & 2.*dat%dt, max(0.2, 1.5-dat%mintime)) ! eps=0.5 first, then 0.2
    
          !call hadvance_plumerise(1,nm1,dt,rad_pC,rad_t,rad_p,mintime)
          !call hadvance_plumerise(2,nm1,dt,rad_pC,rad_t,rad_p,mintime)

          !-- Buoyancy
          if (ifDebugPrint) print*, 'buoyancy'
          call buoyancy_plumerise(NM1, dat%T, dat%TE, dat%QV, dat%QVENV, dat%QH, dat%QI, dat%QC, &
                                & dat%WT, dat%rbuoy)
 
          !-- Entrainment
          if (ifDebugPrint) print*, 'entrainment'
          call entrainment(NM1, dat%W, dat%WT,dat%RADIUS, dat%radius_mdt, &
                         & dat%vel_p, dat%vel_e, dat%zm, dat%entrainment_model, dat%zt, dat%hBL, &
                         & dat%W_treshold, dat%dt_mdt, dat%wind_eff, dat%rho, dat%rho_e, dat%rbuoy, &
                         & dat%C_epsi, dat%C_delta, dat%C_wind, dat%C_delta_wind)

          !-- update W
          if (ifDebugPrint) print*, 'update W'
          call update_plumerise(nm1,'W', dat%DT, &
                              & dat%WT, dat%TT, dat%QVT, dat%QCT, dat%QHT, dat%QIT, dat%VEL_T, dat%SCT, &
                              & dat%W, dat%T, dat%QV, dat%QC, dat%QH, dat%QI, dat%VEL_P, dat%SC)
          
!          call hadvance_plumerise(2,nm1,dat%dt, dat%WC, dat%WT, dat%W, dat%mintime) 
          call predict_plumerise(nm1, dat%wc, dat%w, dat%wt, 2, &
                               & 2.*dat%dt, max(0.2, 1.5-dat%mintime)) ! eps=0.5 first, then 0.2
    
          !call min_max_1D(sc,minsv,maxsv,'noprint')
          !if(minsv .lt. 0.e0) then 
          !    print*,'scalar is negative', minsv    
          !endif

          !-- misc
          if (ifDebugPrint) print*, 'misc'
          do k=2,nm1
              !    pe esta em kpa  - esat_hPa do rams esta em mbar = 100 Pa = 0.1 kpa
              es       = 0.1*esat_hPa(dat%t(k)) !blob saturation vapor pressure, em kPa
              !    rotina do plumegen calcula em kPa
              !    es       = esat_pr (t(k))  !blob saturation vapor pressure, em kPa
              dat%qsat(k) = (eps * es) / (dat%pe(k) - es)  !blob saturation lwc g/g dry air
              dat%est (k) = es 
              dat%txs (k) = dat%t(k) - dat%te(k)
              dat%rho (k) = 3483.8 * dat%pe(k) / dat%t(k) ! air parcel density , g/m**3
                                             ! no pressure diff with radius
              if((abs(dat%wc(k))).gt.wmax) wmax = abs(dat%wc(k)) ! keep wmax largest w

              !srf-27082005
              !     if((abs(wt(k))).gt.wtmax) wtmax = abs(wt(k)) ! keep wmax largest w
 
          enddo  

          ! Gravity wave damping using Rayleigh friction layer for W
          ! Call directly friction, no need for intermediate
          if (ifDebugPrint) print*, 'gravity damp'
!          call friction(2, nm1, deltak, dat%dt, dat%zt, dat%zm, dat%w) 
          call damp_grav_wave(2,nm1, deltak, dat%dt, dat%zt, dat%zm, dat%w, dat%t, &
                            & dat%tt, dat%qv, dat%qh, dat%qi, dat%qc, dat%te, dat%pe, dat%qvenv, &
                            & dat%vel_p, dat%vel_t, dat%vel_e, dat%rad_p, dat%rad_t)

          !++ rp 
          !- update radius for mass conservation
          call mass_conservation(NM1, dat%entrainment_model, dat%rad_p, dat%rad_t, dat%radius, &
                               & dat%radius_mdt, dat%vel_p, dat%zm, dat%zt, dat%w, dat%wc, & 
                               & dat%rho, dat%sc, &
                               & dat%W_treshold, dat%hBl, dat%dzm, dat%dt, dat%hs_rad_p_noDe, dat%VISC, &
                               & dat%dt_mdt, dat%wind_eff, &
                               & dat%rho_e, dat%vel_e, dat%rbuoy, &
                               & dat%C_epsi, dat%C_delta, dat%C_wind, dat%C_delta_wind, &
                               & dat%profile_mass_detrained, dat%profile_mass_entrained)

          ! compute the virtual potential temperature
          Rg = 287.e0
          kappa = Rg/cp
          do k=1,nzPRM
              dat%thv(k) = dat%T(k) / ((dat%pe(k) / dat%pe(1))**kappa) *(1.+.61*dat%qv(k)) ! virtual pot temperature 
              dat%thve(k)= dat%Te(k)/((dat%pe(k) / dat%pe(1))**kappa) *(1.+.61*dat%qvenv(k)) ! virtual pot temperature 
          enddo

          if (ifDebugPrint) print*, 'update radius'
          !do k=1,10 
          !    write(*,'(a,3x,i6.6,4(3x,f12.6))') 'rad', k,rad_p(k), sc(k),rho(k),w(k)
          !enddo
          dat%radius_mdt(:) = dat%radius(:)
          dat%dt_mdt = dat%dt
          !if (flag_print) write(*,'(a,i6.6,a,i2.2)') 'time ', int(time),'_',int(100*(time-int(time)))
          !if (flag_print) write(*,'(3a)')  '         k          zm               rad_p             sc              rho',&
          !                                '              w             hs_rad_p_noDe(1:3,k)                          ',&
          !                                '       radius            thv'
          do k=2,nm1
              if (dat%rad_p(k) .lt. 0.e0 ) then
                  dat%radius(K) = 0.e0
                  dat%rad_p(k) = 0.e0
              endif
              if (dat%radius(k-1) .gt. 1.e0 ) then 
          !        if (flag_print) write(*,'(a,3x,i6.6,10(3x,g14.6))') '    ', k,zm(k), rad_p(k), sc(k),rho(k),w(k),& 
          !                                                             hs_rad_p_noDe(1:3,k),radius(k),thv(k)-thve(k)
                  dat%radius(k) = sqrt(dat%rad_p(k) / ( dat%rho(k) ) )
              else 
                  dat%radius(K) = 0.e0
                  dat%rad_p(k) = 0.e0
              endif
          enddo
  
          !-- try to find the plume top (above surface height)
          if (ifDebugPrint) print*, 'get ztop'
          kk = 1
          do while (dat%w(kk) .gt. 1.)  
              kk = kk + 1  
          enddo 
          kk = max(kk-1, 1)
          dat%ztop =  dat%zm(kk) 
    
          ztopmax = max (dat%ztop, ztopmax)
          dat%ztop_max_global = dat%ztop
          kkmax   = max (kk  , kkmax  ) 

          ! expected mass if no detrainement
          mass_ejected = 1.e-3 * dat%rho(1) * dat%w(1) * pi * dat%radius(1)**2 * dat%dt
          mass_ejected_int = mass_ejected_int + mass_ejected
          mass_ejected_int_sc = mass_ejected_int_sc + mass_ejected * dat%sc(1)

          call min_max_1D(dat%radius, minsv, maxsv, ifDebugPrint)
          if(minsv .lt. 0.e0) then 
              print*,'radius is negative', minsv
              !stop
          endif
          call min_max_1D(dat%sc, minsv, maxsv,ifDebugPrint)
          if(minsv .lt. 0.e0) then 
              print*,'scalar is negative', minsv
              !stop
          endif


          ! mass in the plume
          mass_in_plume = 0.e0
          mass_in_plume_sc = 0.e0
          do k = 2, nm1-1
              mass_in_plume = mass_in_plume + 1.e-3*dat%rho(k)*pi*dat%radius(k)**2*(dat%zm(k)-dat%zm(k-1))
              !write(*,*) 'rad', k,sc(k),radius(k),dat%radius_mdt(k)
              mass_in_plume_sc = mass_in_plume_sc + 1.e-3*dat%rho(k)*pi*dat%radius(k)**2*(dat%zm(k)-dat%zm(k-1)) * dat%sc(k)
          enddo

          if ( (dat%mintime-mintime_sav) .gt. 1.e0 ) then
              var_mass_in_plume = (mass_in_plume - mass_in_plume_mdt) / (mass_in_plume_mdt+1e-10)
              d_mass_in_plume = (mass_in_plume - mass_in_plume_mdt) / dat%dt
              mintime_sav = dat%mintime

              do iii = N_sav,2,-1
                  var_mass_in_plume_ar(iii) = var_mass_in_plume_ar(iii-1)
                  d_mass_in_plume_ar(iii) = d_mass_in_plume_ar(iii-1)
              enddo
              var_mass_in_plume_ar(1) = var_mass_in_plume
              d_mass_in_plume_ar(1) = d_mass_in_plume

          !else
              !var_mass_in_plume = 999
              !d_mass_in_plume = 999
          endif

          var_mass_in_plume = 0.e0
          d_mass_in_plume = 0.e0
          do iii=1,N_sav
              var_mass_in_plume = var_mass_in_plume + var_mass_in_plume_ar(iii)/N_sav
              d_mass_in_plume = d_mass_in_plume + d_mass_in_plume_ar(iii)/N_sav
          enddo
    
          if ( abs(var_mass_in_plume) .le. 2.e-5 .and. dat%mintime .gt. 10 ) then 
              !print*,'============================'
              !print*, mass_in_plume, mass_in_plume_mdt, dt
              !print*, var_mass_in_plume, d_mass_in_plume
              ! stationary = .TRUE.
          endif

          mass_in_plume_mdt = mass_in_plume ! for next time step 

          ! compute entrained and detrained mass
          mass_detrained = 0.e0
          mass_entrained = 0.e0    
          mass_detrained_sc = 0.e0

          ! mass balance at time t
          do k = 2, nm1-1    

             ! SR protection
             if (dat%profile_mass_detrained(k) .le. 0.0) then
               dat% profile_mass_detrained(k)=-1.0*dat%profile_mass_detrained(k)
             endif

             if (dat%profile_mass_entrained(k) .le. 0.0) then
                dat%profile_mass_entrained(k) = -1.0 * dat%profile_mass_entrained(k)
             endif
             mass_detrained = mass_detrained + dat%profile_mass_detrained(k)
             mass_entrained = mass_entrained + dat%profile_mass_entrained(k)
             mass_detrained_sc = mass_detrained_sc + dat%profile_mass_detrained(k) * dat%sc(k) / dat%rho(k) 
             !mass_entrained_sc = mass_entrained_sc + profile_mass_entrained(k)
         enddo
         !print*, mass_detrained_sc
          ! integrate over time
          mass_detrained_int = mass_detrained_int + mass_detrained
          mass_entrained_int = mass_entrained_int + mass_entrained
          mass_detrained_int_sc = mass_detrained_int_sc + mass_detrained_sc
          !mass_entrained_int_sc = mass_entrained_int_sc + mass_entrained_sc
    

          ! integrate profile
          do k = 2, nm1-1
              profile_mass_detrained_int(k) = profile_mass_detrained_int(k) + dat%profile_mass_detrained(k)
              profile_mass_entrained_int(k) = profile_mass_entrained_int(k) + dat%profile_mass_entrained(k)
          enddo

          if(ifDebugPrint) then 
              print*, '---'
              print*, 't = ', dat%mintime, itime
              write(*,*) 'mass Conservation', mass_ejected_int, mass_in_plume,mass_detrained_int, mass_entrained_int
              write(*,*) 'mass in the plume = ', mass_in_plume, var_mass_in_plume, d_mass_in_plume
        
              print*, '---'
          endif
           !if (ifStoreDump) write(1000,'(5(3x,d20.10))') time,mass_ejected_int,mass_in_plume,mass_detrained_int,mass_entrained_int
           if (ifStoreDump) write(1000,*) time,mass_ejected_int,mass_in_plume,mass_detrained_int,mass_entrained_int
           if (ifStoreDump) write(1004,*) time,mass_ejected_int_sc,mass_in_plume_sc,mass_detrained_int_sc

          if (ifDebugPrint)then
            print*, 'call printout'
            call printout (dat, nrectotal, time, uDump)
          endif

          if (ifDebugPrint) print*, 'end loop', dat%w(1)

          itime = itime + 1

          ! stop the time loop if stationary
          if(stationary) exit
 
      ENDDO   !do next timestep
      !++ rp
      ! get the ratio of the mass entrained/detrained to the mass injected by the fire
      do k = 2, nm1-1
          profile_mass_detrained_int(k) = profile_mass_detrained_int(k) / mass_ejected_int
          profile_mass_entrained_int(k) = profile_mass_entrained_int(k) / mass_ejected_int
      enddo



      !
      !++ rp 
      if (ifStoreDump)then
        close(1000)
        close(1004)
        close(1002)
      endif

      !print * ,' ztopmax=',ztopmax,'m',mintime,'mn '
      !print*,'======================================================='
      !if (ifStoreDump) then  
      !    open(1000,file='./mass_profile_en_detrainment.txt')
      !    write(1000,*) '# height   detrain   entrain'
      !    do k = 2, nm1-1
      !        write(*,*) zm(k), profile_mass_detrained_int(k), profile_mass_entrained_int(k)
      !    enddo
      !    close(1000)
      !endif

      !the last printout
      if (ifDebugPrint) then
        call printout (dat, nrectotal, time, uDump)
        close (2) !; close (19)           
      endif

      !--------------------------------------
      !++rp work on the profile of the scalar
      !--------------------------------------
      if (ifStoreDump) then 
          print*,'mass ejected - mass in the plume - mass detrained'
          print*, mass_ejected_int, mass_in_plume,mass_detrained_int, mass_entrained_int
      endif

      RETURN  
      
    END SUBROUTINE MAKEPLUME
    !-------------------------------------------------------------------------------

    !==================================================================================
    !==================================================================================
    !
    ! A bunch of technical subroutines called by the main MAKEPLUME
    !
    !==================================================================================
    !==================================================================================
    
    !***************************************************************************
    !
    SUBROUTINE INITIAL (heat_fluxW, burnt_area, dat, ifDebugPrint, ifStoreDump, uDump)  
      !
      ! Sets initial conditions for the problem.
      ! Zeroying is not needed: it is done automatically at the definition moment
      !
      implicit none 
      real,intent(in) :: heat_fluxW, burnt_area
      type(TPRM_data), intent(inout) :: dat
      logical, intent(in) :: ifDebugPrint, ifStoreDump
      integer, intent(in) :: uDump
      
      real, parameter :: tfreeze = 269.3
      integer ::  isub,  k,  n1,  n2,  n3,  lbuoy,  itmp,  isubm1
      real ::     xn1,  xi,  es
      !
      dat%N = nzPRM
      dat%dt = 1.
      ! initialize temperature structure, to the end of equal spaced sounding,
!      dat%TXS (1:nzPRM) = 0.0  
!      dat%W (1:nzPRM) = 0.0             
      dat%T (1:nzPRM) = dat%TE(:) !blob set to environment		  
!      dat%WC(1:nzPRM) = 0.0
!      dat%WT(1:nzPRM) = 0.0
      dat%QV(1:nzPRM) = dat%QVENV (1:nzPRM)
!      dat%VTH(1:nzPRM) = 0.      !initial rain velocity = 0	                     
!      dat%VTI(1:nzPRM) = 0.      !initial ice  velocity = 0	                     
!      dat%QH(1:nzPRM) = 0.       !no rain				  
!      dat%QI(1:nzPRM) = 0.       !no ice
!      dat%QC(1:nzPRM) = 0.       !no cloud drops	       
      dat%RHO  (1:nzPRM) = 3483.8 * dat%PE(1:nzPRM) / dat%T(1:nzPRM) 	!dry air density g/m**3
      dat%VEL_P(1:nzPRM) = dat%VEL_e(1:nzPRM)
!      dat%sc(1:nzPRM) = 0.e0
!      dat%mass_sc(1:nzPRM) = 0.e0
      dat%rho_e(1:nzPRM) = 3483.8 * dat%PE(1:nzPRM) / dat%TE(1:nzPRM)
      do k = 1, nzPRM
        ES  = 0.1 * ESAT_hPa (dat%T(k)) !blob saturation vapor pressure, kPa
        dat%EST(k) = ES
        dat%QSAT (k) = .622 * ES / (dat%PE(k) - ES) !saturation lwc g/g
      end do
  
      ! Initialize the entrainment radius, Turner-style plume
      dat%radius(1) = dat%rsurf
      do k=2, dat%N
        if (dat%entrainment_model == 1 ) then 
          dat%radius(k) = 0.e0 
          dat%rad_p(k)  = dat%radius(k)**2 * dat%rho(k)
        else
          dat%radius(k) = dat%radius(k-1)+(6./5.)*dat%alpha*(dat%zt(k)-dat%zt(k-1))
          dat%rad_p(k)  = dat%radius(k)
        endif
      enddo
   
      dat%radius_mdt(:) = dat%radius(:)

      !  Initialize the viscosity
      dat%VISC (1) = dat%VISCOSITY
      do k=2, dat%N
        dat%VISC (k) = max(1.e-3, dat%visc(k-1) - 1.* dat%VISCOSITY/float(nzPRM))
      enddo
      !print*,'initial', heat_fluxW
      call LBOUND_MTT(heat_fluxW, 0.e0, dat, ifDebugPrint, ifStoreDump, uDump)
      !call LBOUND_rio(heat_fluxW,0.e0)

      RETURN  
    END SUBROUTINE INITIAL


    !***************************************************************************
 
    subroutine pack_fire_properties(mdur, moist, & 
                                  & alpha, C_epsi, C_delta, C_wind, C_delta_wind, &  ! fire configuration 
                                  & FRP2CHF, &                 !FRP2TOTALH_in,CHF2TOTALH_in, &
                                  & wind_eff, micro_eff, & 
                                  & heat_fluxW, burnt_area, &
                                  & ifPrintDebug, &
                                  & dat)
      !
      ! Just packs the fire information in the data structure defined in the plume_rise_PRMv1
      ! module
      !
      implicit none

      real,intent(in) :: mdur, moist
      real,intent(in) :: alpha, C_epsi, C_delta, C_wind, C_delta_wind
      real,intent(in) :: FRP2CHF                  !FRP2TOTALH_in,CHF2TOTALH_in
      logical,intent(in) :: wind_eff, micro_eff, ifPrintDebug
      real,intent(in) :: heat_fluxW, burnt_area    ! modify by ronan it is now the CHF per m^2
      type(TPRM_data), intent(inout) :: dat

      ! Local variable
      real :: hinc

      dat%mdur     = mdur
      dat%fmoist   = moist * 0.01  ! fuel moisture, fraction
      dat%AREA = burnt_area
      dat%alpha    = alpha
      dat%C_epsi   = C_epsi
      dat%C_delta  = C_delta
      dat%C_wind   = C_wind
      dat%C_delta_wind = C_delta_wind
!      FRP2TOTALH   = FRP2TOTALH_in 
!      CHF2TOTALH   = CHF2TOTALH_in  
      dat%FRP2TOTALH = -999
      dat%CHF2TOTALH = -999  
      dat%wind_eff   = wind_eff
      dat%micro_eff  = micro_eff
      dat%RSURF = SQRT (dat%area / Pi) !- entrainment surface radius (m)
      dat%alpha_dyn = -999
      dat%entrainment_model = 1 ! to be changed in the code at one point

      dat%maxtime = mdur-1  ! model time, min
      IF (MOD (dat%MAXTIME, 2) .NE.0) dat%MAXTIME = dat%MAXTIME+1  !make maxtime even
      if(dat%maxtime > size(dat%heating))then
        print*, dat%maxtime, size(dat%heating)
        STOP 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'
      endif
      dat%MAXTIME = dat%MAXTIME * 60  ! and put in seconds
      !
      ! calculate the energy flux and water content at lboundary.
      ! fills heating() on a minute basis. input has to be adjusted to a one
      ! minute timescale.
      !
      !- make sure of energy release
      dat%HEATING(:) = 0.0001  !- avoid possible divide by 0
      dat%HEATING (1:nint(mdur)) = heat_fluxW       ! W/m**2 (0.55 converte para energia convectiva)

!    bfract = 1.             !- combustion factor
!++ rp EFFLOAD = BLOAD * BFRACT  !- patchy burning
!    TDUR = MDUR * 60.       !- number of seconds in the burn
!    HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2 
!    HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2 

      ! ramp for 5 minutes
      HINC = dat%HEATING (1) / 4.  
      dat%HEATING (1) = 0.1  
      dat%HEATING (2) = HINC  
      dat%HEATING (3) = 2. * HINC  
      dat%HEATING (4) = 3. * HINC 

      if (ifPrintDebug) then  
        PRINT*,'======================================================='
        print * , ' FIRE BOUNDARY CONDITION   :'  
        print * , ' DURATION OF BURN, MINUTES    =',dat%MDUR  
        print * , ' AREA OF BURN, HA             =',dat%AREA*1.e-4
        print * , ' CONVECTIVE HEAT FLUX, kW/m^2 =', heat_fluxW*1.e-3
        !++ rp print * , ' TOTAL LOADING, KG/M**2    =',BLOAD  
        print * , ' FUEL MOISTURE, %             =',dat%FMOIST !average fuel moisture,percent dry
        !print * , ' MODEL TIME, MIN.	      =',MAXTIME  
        !
        print*, 'alpha        = ', dat%alpha
        print*, 'C_epsi       = ', dat%C_epsi
        print*, 'C_delta      = ', dat%C_delta
        print*, 'C_wind       = ',  dat%C_wind
        print*, 'C_delta_wind = ', dat%C_delta_wind
        print*, 'FRP2TOTALH,CHF2TOTALH =', dat%FRP2TOTALH  !, dat%CHF2TOTALH
!        print*, 'FRP2CHF      = ', dat%FRP2CHF
        PRINT*,'======================================================='
      endif

      return
    end subroutine pack_fire_properties


    !******************************************************************    
    
    subroutine tend0_plumerise(wt, tt, qvt, qct, qht, qit, vel_t, rad_t, sct)
      ! Zero the model tendencies
      implicit none
      real, dimension(:), intent(out) :: wt, tt, qvt, qct, qht, qit, vel_t, rad_t, sct
    
      wt(1:nzPRM)  = 0.
      tt(1:nzPRM)  = 0.
      qvt(1:nzPRM)  = 0.
      qct(1:nzPRM)  = 0.
      qht(1:nzPRM)  = 0.
      qit(1:nzPRM)  = 0.
      vel_t(1:nzPRM)  = 0.
      rad_t(1:nzPRM)  = 0.
      sct(1:nzPRM)  = 0.
    end subroutine tend0_plumerise
    
    
    !***************************************************************************
    
    subroutine vel_advectc_plumerise(m1,wc,wt,rho,dzm)
      !
      ! Compute advection contribution to W tendency
      !
      implicit none
      integer, intent(in) :: m1
      real, dimension(:), intent(in) :: wc,dzm,rho
      real, dimension(:), intent(inout) :: wt
      
      ! local variables
      real, dimension(m1) :: flxw, dn0 ! var local
      real :: c1z
      integer :: k

      dn0(1:m1) = rho(1:m1) * 1.e-3 ! converte de cgs para mks
      flxw(1) = wc(1) * dn0(1) 
      do k = 2,m1-1
         flxw(k) = wc(k) * .5 * (dn0(k) + dn0(k+1))
      enddo
      ! Compute advection contribution to W tendency
      c1z = .5 
      do k = 2,m1-2
         wt(k) = wt(k) + c1z * dzm(k) / (dn0(k) + dn0(k+1)) * &
	                   & ((flxw(k) + flxw(k-1))  * (wc(k) + wc(k-1)) - &
                        & (flxw(k) + flxw(k+1))  * (wc(k) + wc(k+1)) + &
                        & (flxw(k+1) - flxw(k-1)) * 2.* wc(k)       )
      enddo
      return
    end subroutine vel_advectc_plumerise 
    
    
    !***********************************************************************************
    !
    subroutine scl_advectc_plumerise(varn, mzp,ztopmax, dt, T, QV, QC, QI, QH, VEL_P, SC, &
                                                        & TT, QVT, QCT, QIT, QHT, VEL_T, SCT, &
                                    & w, wc, rho, dzm, dzt, zt, zm, vt3dc, vt3df, vt3dg, vt3dk)
      !
      ! Scalar-advection for the plume
      !
      implicit none
      
      integer, intent(in) :: mzp
      character(len=*), intent(in) :: varn
      real,intent(in) :: ztopmax, dt
      real, dimension(:), intent(in) :: T, QV, QC, QI, QH, VEL_P, SC, w, wc, rho, dzm, dzt, zt, zm
      real, dimension(:), intent(inout) :: vt3dc, vt3df, vt3dg, vt3dk, &
                                         & TT, QVT, QCT, QIT, QHT, VEL_T, SCT

      ! local variables
      real :: dtlto2
      integer :: k
      real, dimension(mzp) :: vctr1, vctr2, scr1

      !  wp => w
      !- Advect  scalars
      dtlto2   = .5 * dt
      !  vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * dne(1)
      vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * rho(1)*1.e-3!converte de CGS p/ MKS
      vt3df(1) = .5 * (w(1) + wc(1)) * dtlto2 * dzm(1)

      do k = 2,mzp
      !     vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (dne(k) + dne(k+1))
      !++ rp
        if (k.lt. mzp ) then 
          vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (rho(k) + rho(k+1))*1.e-3
        else
          vt3dc(k) =  (w(k) + wc(k)) * dtlto2 * (rho(k))*1.e-3
        endif
      !--rp
        vt3df(k) =  (w(k) + wc(k)) * dtlto2 *.5 *  dzm(k)
           !print*,'vt3df-vt3dc',k,vt3dc(k),vt3df(k)
      enddo

       !  do k = 1,mzp-1
      do k = 1,mzp
        if (k .lt. mzp) then 
          vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
        else
          vctr1(k) = 0.e0
        endif
        vctr2(k) = (zm(k)   - zt(k)) * dzm(k)
      !    vt3dk(k) = dzt(k) / dne(k)
        vt3dk(k) = dzt(k) /(rho(k)*1.e-3)
            !print*,'VT3dk',k,dzt(k) , dne(k)
      enddo

      !      scalarp => scalar_tab(n,ngrid)%var_p
      !      scalart => scalar_tab(n,ngrid)%var_t
      !- temp advection tendency (TT)
      scr1=T
      call fa_zc_plumerise(mzp, T, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      call advtndc_plumerise(mzp,T,scr1,TT,dt)

      !- water vapor advection tendency (QVT)
      scr1=QV
      call fa_zc_plumerise(mzp,QV, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      call advtndc_plumerise(mzp,QV,scr1,QVT,dt)

      !- liquid advection tendency (QCT)
      scr1=QC
      call fa_zc_plumerise(mzp, QC, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      call advtndc_plumerise(mzp,QC,scr1,QCT,dt)

      !- ice advection tendency (QIT)
      scr1=QI
      call fa_zc_plumerise(mzp, QI,scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      call advtndc_plumerise(mzp,QI,scr1,QIT,dt)

      !- hail/rain advection tendency (QHT)
      !   if(ak1 > 0. .or. ak2 > 0.) then
      scr1=QH
      call fa_zc_plumerise(mzp, QH,scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      call advtndc_plumerise(mzp,QH,scr1,QHT,dt)
      !   endif

      !- horizontal wind advection tendency (VEL_T)
      scr1=VEL_P
      call fa_zc_plumerise(mzp, VEL_P, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

            !call advtndc_plumerise(mzp,VEL_P,scr1,VEL_T,dt)
            !advtndc_plumerise is just down to be able to control the clip on the
            !outside plume where we do not solbe the equation for vel_p
      do k = 2,mzp-1
        if (zt(k) .lt. ztopmax) vel_t(k) = vel_t(k) + (scr1(k)-vel_p(k)) * 1.e0/dt
      enddo

      !- vertical radius transport
      !      scr1=rad_p
      !      call fa_zc_plumerise(mzp                  &
      !             	          ,rad_p     ,scr1    &
      !             	          ,vt3dc  ,vt3df   &
      !             	          ,vt3dg  ,vt3dk   &
      !             	          ,vctr1,vctr2	       )
      !
      !      call advtndc_plumerise(mzp,rad_p,scr1,rad_t,dt)
      !   return

      !- gas/particle advection tendency (SCT)
      !    if(varn == 'SC')return
      scr1=SC
      call fa_zc_plumerise(mzp, SC, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)
   
      call advtndc_plumerise(mzp,SC,scr1,SCT,dt)

      return

        CONTAINS

          !------------------------------------------------------------------------

          subroutine fa_zc_plumerise(m1, scp, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)
            implicit none
            integer, intent(in) :: m1
            real, dimension(:), intent(in) :: scp,vt3dc,vt3df,vt3dk
            real, dimension(:), intent(inout) :: vctr1,vctr2,scr1,vt3dg

            ! Local variables
            integer :: k
            real, parameter :: dfact = 0.5
            
            ! Compute scalar flux VT3DG
            do k = 1,m1-1
               vt3dg(k) = vt3dc(k) * (vctr1(k) * scr1(k) + &
                                    & vctr2(k) * scr1(k+1) +   &
                                    & vt3df(k) * (scr1(k) - scr1(k+1)))
            enddo
      
            ! Modify fluxes to retain positive-definiteness on scalar quantities.
            !    If a flux will remove 1/2 quantity during a timestep,
            !    reduce to first order flux. This will remain positive-definite
            !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
            !    both fluxes are evacuating the box.

            do k = 1,m1-1
              if (vt3dc(k) > 0.) then
                if (vt3dg(k) * vt3dk(k) > dfact * scr1(k)) vt3dg(k) = vt3dc(k) * scr1(k)
              else
                if (-vt3dg(k) * vt3dk(k+1) > dfact * scr1(k+1)) vt3dg(k) = vt3dc(k) * scr1(k+1)
              endif
            enddo

            ! Compute flux divergence
            do k = 2,m1-1
              scr1(k) = scr1(k) + vt3dk(k) * (vt3dg(k-1) - vt3dg(k) + scp(k) * (vt3dc(k) - vt3dc(k-1)))
            enddo
            return
          end subroutine fa_zc_plumerise
          
          !-------------------------------------------------------------------------------
 
          subroutine advtndc_plumerise(m1, scp, sca, sct, dtl)
            implicit none
            integer, intent(in) :: m1
            real, intent(in) :: dtl
            real, dimension(:), intent(in) :: scp,sca
            real, dimension(:), intent(inout) :: sct

            sct(2:m1-1) = sct(2:m1-1) + (sca(2:m1-1)-scp(2:m1-1)) / dtl
            
            return
          end subroutine advtndc_plumerise
 
    end subroutine scl_advectc_plumerise
 
    !******************************************************************************

    
    subroutine visc_W(m1,deltak,ztopmax, zt, zm, W, T, QC, QV, QH, QI, SC, VEL_P, &
                                                  & WT, TT, QCT, QVT, QHT, QIT, SCT, VEL_T,  &
                                                  & RAD_P, visc)
      implicit none
      
      integer, intent(in) :: m1,deltak
      real, intent(in) :: ztopmax
      real, dimension(:), intent(in) :: zt, zm, W, T, QC, QV, QH, QI, SC, VEL_P, RAD_P, visc
      real, dimension(:), intent(inout) :: WT, TT, QCT, QVT, QHT, QIT, SCT, VEL_T
      
      ! Local variables
      integer :: k,m2
      real :: dz1t,dz1m,dz2t,dz2m,d2wdz,d2tdz  ,d2qvdz ,d2qhdz ,d2qcdz ,d2qidz ,d2scdz, &
            & d2vel_pdz,d2rad_dz
      integer :: mI,mF,deltaM

      mI = 2
      mF = min(m1,nzPRM-1)
      deltaM = 1

      do k=mI,mF,deltaM !v2
      !do k=2,m2-1 !orig

          DZ1T   = 0.5*(ZT(K+1)-ZT(K-1))
          DZ2T   = VISC (k) / (DZ1T * DZ1T)  
          DZ1M   = 0.5*(ZM(K+1)-ZM(K-1))
          DZ2M   = VISC (k) / (DZ1M * DZ1M)  
    
          D2WDZ  = (W  (k + 1) - 2 * W  (k) + W  (k - 1) ) * DZ2M  
          D2TDZ  = (T  (k + 1) - 2 * T  (k) + T  (k - 1) ) * DZ2T  
          D2QVDZ = (QV (k + 1) - 2 * QV (k) + QV (k - 1) ) * DZ2T  
          D2QHDZ = (QH (k + 1) - 2 * QH (k) + QH (k - 1) ) * DZ2T 
          D2QCDZ = (QC (k + 1) - 2 * QC (k) + QC (k - 1) ) * DZ2T  
          D2QIDZ = (QI (k + 1) - 2 * QI (k) + QI (k - 1) ) * DZ2T  
          D2SCDZ = (SC (k + 1) - 2 * SC (k) + SC (k - 1) ) * DZ2T 
          d2vel_pdz=(vel_P  (k + 1) - 2 * vel_P  (k) + vel_P  (k - 1) ) * DZ2T
          d2rad_dz =(rad_p  (k + 1) - 2 * rad_p  (k) + rad_p  (k - 1) ) * DZ2T 
    
          WT(k) =   WT(k) + D2WDZ
          TT(k) =   TT(k) + D2TDZ
                            
          QVT(k) =  QVT(k) + D2QVDZ 
          QCT(k) =  QCT(k) + D2QCDZ
          QHT(k) =  QHT(k) + D2QHDZ 
          QIT(k) =  QIT(k) + D2QIDZ     
    
          if (zt(k) .lt. ztopmax) VEL_T(k) =   VEL_T(k) + d2vel_pdz 
          ! apply viscosity only if we are in the plume
    
          !rad_t(k) =   rad_t(k) + d2rad_dz

          SCT(k) =  SCT(k) + D2SCDZ
          !print*,'W-VISC=',k,D2WDZ
      enddo  
    end subroutine visc_W
        
        
    !************************************************************************************
        
    subroutine update_plumerise(m1,varn, DT, WT, TT, QVT, QCT, QHT, QIT, VEL_T, SCT, W, T, QV, QC, QH, QI, VEL_P, SC)
      !
      ! Computes new variable values form their tendencies
      !
      implicit none

      integer :: m1
      real, intent(in) :: DT
      character(len=*), intent(in) :: varn
      real, dimension(:), intent(inout) :: W, T, QV, QC, QH, QI, VEL_P, SC
      real, dimension(:), intent(in) :: WT, TT, QVT, QCT, QHT, QIT, VEL_T, SCT
      
      ! local variables
      integer :: k

      if(varn == 'W') then
        W(2:m1-1) = W(2:m1-1) +  WT(2:m1-1) * DT  
      else
        do k=2,m1-1
            T(k) =  T(k) +  TT(k) * DT  
            QV(k) = max(0., QV(k) + QVT(k) * DT)
            QC(k) = max(0., QC(k) + QCT(k) * DT) !cloud drops travel with air 
            QH(k) = max(0., QH(k) + QHT(k) * DT)
            QI(k) = max(0., QI(k) + QIT(k) * DT)
            VEL_P(k) =  VEL_P(k) + VEL_T(k) * DT  

            !rad_p(k) =  rad_p(k) + rad_t(k) * DT
            !++ rp : it is now done at the end of the temproal loop
            !rhs_sct_max =   -1.e0 * ( sc(k)/dt )  
            !clip the scalar to avoid to remove too much 
            !SCT(k) = max(  SCT(k), rhs_sct_max) 

            SC(k)  =  SC(k) +  SCT(k) * dt
        enddo
      endif
    end subroutine update_plumerise 
        

    !************************************************************************

    subroutine predict_plumerise(npts, ac, ap, fa, iac, dtlp, epsu)
      implicit none
      integer :: npts,iac,m
      real :: epsu,dtlp
      real, dimension(*), intent(inout) :: ac,ap,fa
      
      ! local variables
      real, dimension(npts) :: af

      !     For IAC=3, this routine moves the arrays AC and AP forward by
      !     1 time level by adding in the prescribed tendency. It also
      !     applies the Asselin filter given by:
      !              {AC} = AC + EPS * (AP - 2 * AC + AF)
      !     where AP,AC,AF are the past, current and future time levels of A.
      !     All IAC=1 does is to perform the {AC} calculation without the AF
      !     term present.  IAC=2 completes the calculation of {AC} by adding
      !     the AF term only, and advances AC by filling it with input AP
      !     values which were already updated in ACOUSTC.
      !
      if (iac .eq. 1) then
         do m = 1,npts
            ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
         enddo
         return
      elseif (iac .eq. 2) then
         do m = 1,npts
            af(m) = ap(m)
            ap(m) = ac(m) + epsu * af(m)
         enddo
      !elseif (iac .eq. 3) then
      !   do m = 1,npts
      !      af(m) = ap(m) + dtlp * fa(m)
      !   enddo
      !   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic(nzp,nxp,nyp,af,'T')
      !   do m = 1,npts
      !      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
      !   enddo
      endif
      ac(1:npts) = af(1:npts)
      return
    end subroutine predict_plumerise     


    !***************************************************************************

    subroutine  buoyancy_plumerise(m1, T, TE, QV, QVENV, QH, QI, QC, WT, rbuoy)
      implicit none
      integer, intent(in) :: m1
      real, dimension(:), intent(inout) :: T, TE, QV, QVENV, QH, QI, QC, WT, rbuoy
      
      real :: TV,TVE,QWTOTL,umgamai
      real, parameter :: mu = 0.15 
      integer :: k
      real, dimension(m1) :: scr1

      !- orig
      umgamai = 1./(1.+ 0.5) !gama) ! ompensates for the lack of the acceleration term associated with
                                    ! non-hydrostatic disturbances in the pressure field
      !- new                 ! Siesbema et al, 2004
      !umgamai = 1./(1.-2.*mu)

      do k = 2,m1-1
        TV =   T(k) * (1. + (QV(k)   /EPS))/(1. + QV(k)   )  !blob virtual temp.                                        	   
        TVE = TE(k) * (1. + (QVENV(k)/EPS))/(1. + QVENV(k))  !and environment

        QWTOTL = QH(k) + QI(k) + QC(k)                       ! QWTOTL*G is drag
        scr1(k)= G*  umgamai*( (TV - TVE) / TVE   - QWTOTL) 
      enddo

      rbuoy(:) = 0 
      do k = 2,m1-2
          !srf- just for output
          rbuoy(k)=0.5*(scr1(k)+scr1(k+1))
          wt(k) = wt(k)+0.5*(scr1(k)+scr1(k+1))
      enddo

      !srf- just for output
      rbuoy(1)=rbuoy(2)
    
    end subroutine  buoyancy_plumerise
  

    !***************************************************************************

    subroutine ENTRAINMENT(m1, w, wt, radius, radius_mdt, vel_p, vel_e, zm, entrainment_model, zt, hBL, &
                         & W_treshold, dt_mdt, wind_eff, rho, rho_e, rbuoy, &
                         & C_epsi, C_delta, C_wind,C_delta_wind)
      !
      ! COmputes entrainment of freh ir into the plume
      !
      implicit none

      real, dimension(:), intent(in) :: w, vel_p, vel_e, zm, zt, rho, rho_e, radius,radius_mdt, rbuoy
      real, dimension(:), intent(inout) :: wt
      integer, intent(in) :: m1, entrainment_model
      real, intent(in) :: hBL, W_treshold, dt_mdt, C_epsi, C_delta, C_wind,C_delta_wind
      logical, intent(in) :: wind_eff
      
      ! local varibles
      real, parameter :: mu = 0.15 ,gama = 0.5 ! mass virtual coeff.
      integer :: k
      real :: WBAR, epsi, delta, delta_wind, epsi_wind, epsi_tot

      do k=2,m1-1
        WBAR=W(k)          
        call entrainment_coeff (entrainment_model, wind_eff, 'W', W_treshold, hBl, m1, &
                                & radius,radius_mdt,zm,zt,w,rho,rho_e, vel_p, vel_e, & 
                                & dt_mdt,k, rbuoy,  &
                                & epsi, delta, epsi_wind, delta_wind, &
                                & C_epsi, C_delta, C_wind,C_delta_wind)
        if (wind_eff) then 
          epsi_tot = epsi + epsi_wind
          delta = delta + delta_wind
        else
          epsi_tot = epsi 
        endif
        wt(k) = wt(k)  - epsi_tot * WBAR
      enddo
    end subroutine  ENTRAINMENT

                         
    !******************************************************************************************************
                         
    subroutine scl_misc(m1, entrainment_model, wind_eff, W_treshold, hBl, &
                      & T, TE, QV, QC, QH, QI, QVENV, SC, VEL_P, VEL_E, &
                      & TT, QVT, QCT, QHT, QIT, SCT, VEL_T, &
                      & radius, radius_mdt, zm, zt, w, rho, rho_e, dt_mdt, rbuoy, &
                      & C_epsi, C_delta, C_wind,C_delta_wind, WBAR, dwdt_entr)

      implicit none

      integer,intent(in) :: m1, entrainment_model
      logical, intent(in) :: wind_eff
      real, intent(in) :: W_treshold,hBl,dt_mdt, C_epsi, C_delta, C_wind,C_delta_wind
      real, dimension(:), intent(in) :: T, TE, QV, QC, QH, QI, QVENV, SC, VEL_P, VEL_E
      real, dimension(:), intent(inout) :: TT, QVT, QCT, QHT, QIT, SCT, VEL_T
      real, dimension(:), intent(in) :: radius,radius_mdt,zm,zt,w,rho,rho_e,rbuoy
      real, intent(out) :: WBAR
      real, dimension(:,:), intent(out) :: dwdt_entr
      
      integer :: k
      real dmdtm, ADIABAT
      !++ rp 
      real :: epsi,delta, delta_wind, source_vel,epsi_tot,epsi_wind

      do k=2,m1-1
  
          WBAR    = 0.5*(W(k)+W(k-1))  
          !-- dry adiabat
          ADIABAT = - WBAR * G / CP 

          !-- entrainment    
          call entrainment_coeff (entrainment_model,wind_eff,'S',W_treshold,hBl,m1, &
                                  radius,radius_mdt,zm,zt,w,rho,rho_e,vel_p, vel_e, & 
                                  dt_mdt,k,  rbuoy,  &
                                  epsi,delta, epsi_wind, delta_wind, &
                                & C_epsi, C_delta, C_wind,C_delta_wind)
          dwdt_entr(1,k) = epsi 
          dwdt_entr(2,k) = delta
          dwdt_entr(3,k) = epsi_wind 
          dwdt_entr(4,k) = delta_wind

          !if(k == 2)print*, epsi_wind

          if (wind_eff) then 
              epsi_tot = epsi + epsi_wind
              delta = delta + delta_wind
          else
              epsi_tot = epsi 
          endif

          !-- tendency temperature = adv + adiab + entrainment
          TT(k) = TT(K) + ADIABAT - epsi_tot * ( T  (k) -    TE (k) ) 

          !-- tendency water vapor = adv  + entrainment
          QVT(K) = QVT(K)         - epsi_tot * ( QV (k) - QVENV (k) ) 

          QCT(K) = QCT(K)	      - epsi_tot * ( QC (k)  ) 
          QHT(K) = QHT(K)	      - epsi_tot * ( QH (k)  ) 
          QIT(K) = QIT(K)	      - epsi_tot * ( QI (k)  )

          !-- tendency horizontal speed = adv  + entrainment
          source_vel = 0.e0
          VEL_T(K) = VEL_T(K)  - epsi_tot * ( VEL_P (k) - VEL_E (k) )  + source_vel 

          !-- tendency gas/particle = adv  + entrainment
          SCT(K) = SCT(K) - epsi_tot * ( SC (k) ) 

       enddo
    end subroutine scl_misc
 
                      
!    !*************************************************************************************
!
!    subroutine friction(ifrom, nm1, deltak, dt, zt, zm, var1, vart, var2)
!      implicit none
!      integer, intent(in) :: nm1, ifrom, deltak
!      real, dimension(:), intent(inout) :: var1,zt,zm
!      real, dimension(:), intent(inout), optional :: var2,vart
!
!      ! local variables
!      integer :: k
!      real zmkf,ztop,distim,c1,c2,dt
!
!      !kf = nm1 - int(deltak/2) ! orig
!      !if( kf < 10) return !ver necessidade
!      zmkf = zm(nm1 - int(deltak)) !old: float(kf )*dz
!      ztop = zm(nm1)
!
!      !distim = 60. ! orig
!      distim = min(3.*dt,60.)
!
!      c1 = 1. / (distim * (ztop - zmkf))
!      c2 = dt * c1
!
!      if(ifrom == 1) then  
!        do k = nm1,2,-1
!         if (zt(k) .le. zmkf) cycle  ! exit ???
!         vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
!        enddo
!      elseif(ifrom == 2) then
!        do k = nm1,2,-1
!         if (zt(k) .le. zmkf) cycle  ! exit ???
!         var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
!        enddo
!      endif
!      return
!    end subroutine friction 


    !**************************************************************************************
    
    subroutine mass_conservation(m1,entrainment_model,rad_p,rad_t,radius, radius_mdt,vel_p,& 
                               & zm,zt,w,wc,rho,sc, &
                               & W_treshold, hBl, dzm, dt, hs_rad_p_noDe, VISC, dt_mdt, wind_eff, &
                               & rho_e, vel_e, rbuoy,  C_epsi, C_delta, C_wind,C_delta_wind, &
                               & mass_detrainment,mass_entrainment)

      implicit none

      integer,intent(in) :: m1,entrainment_model
      real,dimension(m1),intent(inout) :: rad_p,rad_t
      real,dimension(m1),intent(in) :: radius,radius_mdt,zm,zt,w,wc,rho,sc,vel_p
      real,dimension(m1),intent(out) :: mass_detrainment,mass_entrainment
      real, intent(in) :: W_treshold, hBl, dt, dt_mdt
      real, dimension(:), intent(in) :: rho_e, vel_e, rbuoy, dzm, VISC
      real, dimension(:,:), intent(inout) :: hs_rad_p_noDe
      real, intent(in) ::  C_epsi, C_delta, C_wind,C_delta_wind
      logical, intent(in) :: wind_eff

      ! local variable
      real :: rhom1drhodz,wm1dwdz, dtlto2
      integer :: k
      real :: WBAR,epsi, delta, delta_wind, delta_Aloc,epsi_wind
      real :: DMDTM, DMDTM_sc, DMDTM_dyn
      real :: rhs_rad_p
      real :: minsv,maxsv,mass_loss,rhs_rad_p_noDe,rhs_rad_p_init,mass_gain
      real :: tau_detrainment,DMDTM_detrain,sign_w,tempkp05,tempkm05,drhoR2dz
      real :: radiuskp05,radiuskm05,sckp05,sckm05,rhokp05,rhokm05
      real :: DZ1T,DZ2T,d2rad_dz,temp,temp_raw,rad_t_Add_Mass,rad_t_Rm_Mass
      real :: entr, entr_min_limit,entr_raw
      real :: rad_p_kp05, rad_p_km05

      tau_detrainment = 10.e0
      mass_detrainment(:) = 0.e0
      mass_entrainment(:) = 0.e0
      !
      ! init
      !
      rad_t(:) = 0.e0
      !
      ! add viscosity 
      !
      do k=2,m1-1
        DZ1T   = 0.5*(ZT(K+1)-ZT(K-1))
        DZ2T = VISC (k) / (DZ1T * DZ1T)
        d2rad_dz = (rad_p(k + 1) - 2 * rad_p(k) + rad_p(k - 1) ) * DZ2T
        rad_t(k) =   rad_t(k) + d2rad_dz
      enddo
      !
      ! entrainement
      !
      do k=2,m1-1
        WBAR    = 0.5*(W(k)+W(k-1))  
!        if (k .gt. 2) then 
!            radiuskp05 = 0.5*(radius(k)+radius(k+1))
!            radiuskm05 = 0.5*(radius(k)+radius(k-1))
!            rhokp05 = 0.5*(rho(k)+rho(k+1))
!            rhokm05 = 0.5*(rho(k)+rho(k-1))
!        else
!            radiuskp05 = 0.5*(radius(k)+radius(k+1))
!            radiuskm05 = radius(k-1)
!            rhokp05 = 0.5*(rho(k)+rho(k+1))
!            rhokm05 = rho(k-1)
!        endif
        if (k .gt. 2) then 
          rad_p_km05 = 0.5 * (rad_p(k-1) + rad_p(k)  )
          rad_p_kp05 = 0.5 * (rad_p(k)   + rad_p(k+1))
          !radiuskp05 = 0.5 * (radius(k+1) + radius(k))
          !radiuskm05 = 0.5 * (radius(k-1) + radius(k))
        else
          rad_p_km05 = rad_p(k-1)  
          rad_p_kp05 = 0.5 * (rad_p(k)   + rad_p(k+1))
          !radiuskp05 = 0.5 * (radius(k+1) + radius(k))
          !radiuskm05 = radius(k-1)
        endif
        tempkp05 = w(k)   * rad_p_kp05
        tempkm05 = w(k-1) * rad_p_km05     
!        tempkp05 = w(k)   * (radiuskp05)**2 * rhokp05
!        tempkm05 = w(k-1) * (radiuskm05)**2 * rhokm05  

        call entrainment_coeff (entrainment_model, wind_eff, 'S',W_treshold, hBl, m1,    &
                                radius,radius_mdt,zm,zt,w,rho, rho_e ,vel_p, vel_e,        & 
                                dt_mdt,k, rbuoy,                                   &
                                epsi, delta, epsi_wind,delta_wind, &
                                & C_epsi, C_delta, C_wind,C_delta_wind           )
        if (wind_eff) then 
          epsi = epsi + epsi_wind
          delta = delta + delta_wind
        endif
        !
        ! compute rhs
        !
        !drhoR2dz = (tempkp05 - tempkm05) / (zm(k)-zm(k-1))
        drhoR2dz = -1 * (tempkp05 - tempkm05) * dzm(k)
        wm1dwdz     = ( w(k) - w(k-1) ) / (zm(k) - zm(k-1) ) 
        sign_w = sign(1.e0, wm1dwdz )

        !for debugging
        !-------------
        hs_rad_p_noDe(1,k) = -1.e0 * drhoR2dz 
        hs_rad_p_noDe(2,k) = rho(k) * radius(k)**2 * (epsi - delta) 
        hs_rad_p_noDe(3,k) = rad_p(k) / DT
        
        ! add rhs to the mass variable rad_p and clip to keep positive rad_p
        !-----------------------------------
        entr_raw = rad_p(k) * (epsi - delta)
        entr_min_limit =  -1. * ( rad_p(k) / DT + rad_t(k) + drhoR2dz )
        entr = max(entr_raw, entr_min_limit )
        !temp = max( rho(k) * radius(k)**2 * (epsi - delta),   & 
        !            -1 * (-1.e0 * drhoR2dz + rad_p(k) / DT)   )
       
        !write(*,'(3(3x,g15.6))')   temp_raw,    -1.e0 * drhoR2dz , rad_p(k) / DT     
        !write(*,'(4(3x,g15.6))'), temp, rho(k) , radius(k)**2, epsi
        temp = entr-entr_raw
        if (temp .ne. 0 ) then 
!          if (temp-temp_raw .ne. 0 ) then 
          if (temp .gt. 0 ) then 
            !epsi =        ( temp / (rho(k) * (0.5*(radiuskm05+radiuskp05))**2)  + delta )
            epsi = ( entr / ( 0.5* (rad_p_km05+rad_p_kp05) )  + delta )
          else
            !delta =  -1 * ( temp / (rho(k) * radius(k)**2)  - epsi )
            delta = ( entr / ( rad_p(k) )  - epsi )
          endif
        endif

        !rhs_rad_p_noDe =   -1.e0 * drhoR2dz + entr
        !rad_t(K) = rad_t(K)      +   rhs_rad_p_noDe 

        ! for mass variation monitoring
        !----------------------------------
        if (rad_t(K) .gt. 0 )  then 
            rad_t_Add_Mass =  rad_t(K) * pi * 1.e-3 *(zm(k+1)-zm(k)) * dt 
        else
            rad_t_Add_Mass = 0
        endif

        if (rad_t(K) .lt. 0 ) then 
            rad_t_Rm_Mass  =  -1*rad_t(K) * pi * 1.e-3 *(zm(k+1)-zm(k)) * dt
        else
            rad_t_Rm_Mass   = 0
        endif
        mass_loss = rho(k) * radius(k)**2 * delta * pi * 1.e-3 *(zm(k)-zm(k-1)) + rad_t_Rm_Mass
        mass_gain = rho(k) * radius(k)**2 * epsi  * pi * 1.e-3 *(zm(k)-zm(k-1)) + rad_t_Add_Mass
        mass_detrainment(k) =  mass_loss 
        mass_entrainment(k) =  mass_gain 
        !
        ! advance in time
        !
        rad_p(k) =   rad_p(k) + (rad_t(k) + drhoR2dz + entr ) * DT 
        if( abs(rad_p(k)) .lt. 1.e-6 ) rad_p(k) = 0.e0
      enddo   ! over vertical

      !call min_max_1D(m1,rad_p,minsv,maxsv,'print')

    end subroutine mass_conservation


    !******************************************************************
                               
    subroutine damp_grav_wave(ifrom,nm1,deltak,dt,zt,zm,w,t,tt,qv,qh,qi,qc,te,pe,qvenv&
                         ,vel_p,vel_t,vel_e,rad_p,rad_t)
      implicit none
      integer nm1,ifrom,deltak
      real dt
      real, dimension(:) :: w,t,tt,qv,qh,qi,qc,te,pe,qvenv,zt,zm &
                             ,vel_p,vel_t,vel_e,rad_p,rad_t

      real, dimension(nzPRM) :: dummy
      integer :: ifrom_loc                       
      dummy(:) = 0.

      !for the temperature
      if(ifrom==1) then
        call friction(ifrom,nm1,deltak,dt,zt,zm,t,tt    ,te)
        ! call friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
        ! call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
        return
      endif 

      ! for the velocity 
      if(ifrom==2) then 
        call friction(ifrom,nm1,deltak,dt,zt,zm,w,dummy ,dummy)
        !call friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
        !call friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
        !call friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
        return
      endif

!++ rp
! add damping for the radius
      if(ifrom==3) then
        ifrom_loc = 1
        call friction(ifrom_loc,nm1,deltak,dt,zt,zm,rad_p,rad_t,dummy)
        ! call friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
        ! call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
        return
      endif 

      CONTAINS

        subroutine friction(ifrom,nm1,deltak,dt,zt,zm,var1,vart,var2)
          implicit none
          real, dimension(:) :: var1,var2,vart,zt,zm
          integer k,nfpt,kf,nm1,ifrom,deltak
          real zmkf,ztop,distim,c1,c2,dt
!nfpt=50
!kf = nm1 - nfpt
!kf = nm1 - int(deltak/2) ! orig
          kf = nm1 - int(deltak)
!if( kf < 10) return !ver necessidade

          zmkf = zm(kf) !old: float(kf )*dz
          ztop = zm(nm1)

          !distim = 60. ! orig
          distim = min(3.*dt,60.)

          c1 = 1. / (distim * (ztop - zmkf))
          c2 = dt * c1

          if(ifrom == 1) then  
            do k = nm1,2,-1
             if (zt(k) .le. zmkf) cycle  ! exit ???
             vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
            enddo
          elseif(ifrom == 2) then
            do k = nm1,2,-1
             if (zt(k) .le. zmkf) cycle  ! exit ???
             var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
            enddo
          endif
          return
        end subroutine friction
!-------------------------------------------------------------------------------
     end subroutine damp_grav_wave    


    !******************************************************************************************************
    
    subroutine fallpart(m1, RHO, W, CVH, CVI, VTH, VHREL, VTI, VIREL, QH, QI, ZM, QHT, QIT)
      !   verificar se o gradiente esta correto 
      !  
      !     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
      !     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
      !     for ice, see (OT18), use F0=0.75 per argument there. rho*q
      !     values are in g/m**3, velocities in m/s
      implicit none

      integer, intent(in) :: m1
      real, intent(inout) :: VHREL, VIREL
      real, dimension(:), intent(in) :: RHO, W, QH, QI, ZM
      real, dimension(:), intent(inout) :: QHT, QIT, CVH, VTH, VTI
      real, dimension(:), intent(out) :: CVI
      
      ! Local variables
      integer :: k
      real :: vtc, dfhz,dfiz,dz1
      real, PARAMETER :: VCONST = 5.107387, F0 = 0.75  

      do k=2,m1-1
        VTC = VCONST * RHO(k) **.125   ! median volume fallspeed (KTable4)
                                
      !  hydrometeor assembly velocity calculations (K Table4)
      !  VTH(k)=-VTC*QH(k)**.125  !median volume fallspeed, water            
        VTH (k) = - 4.	    !small variation with qh
   
        VHREL = W(k) + VTH(k)  !relative to surrounding cloud
 
      !  rain ventilation coefficient for evaporation
        CVH(k) = 1.6 + 0.57E-3 * (ABS (VHREL) ) **1.5  

      !  VTI(k)=-VTC*F0*QI(k)**.125    !median volume fallspeed,ice             
        VTI(k) = - 3.                !small variation with qi

        VIREL = W(k) + VTI(k)       !relative to surrounding cloud

      !  ice ventilation coefficient for sublimation
        CVI(k) = 1.6 + 0.57E-3 * (ABS (VIREL) ) **1.5 / F0  

        IF (VHREL.GE.0.0) THEN  
          DFHZ = QH(k) * (RHO(k) *   VTH(k) - RHO(k-1) * VTH(k-1)) / RHO(k-1)
        ELSE  
          DFHZ = QH(k) * (RHO(k+1) * VTH(k+1) - RHO(k ) * VTH(k  )) / RHO(k)
        ENDIF  

        IF (VIREL.GE.0.0) THEN  
          DFIZ = QI(k) * (RHO(k  ) * VTI(k  ) - RHO(k-1) * VTI(k-1)) / RHO(k-1)
        ELSE  
          DFIZ = QI(k) * (RHO(k+1) * VTI(k+1) - RHO(k  ) * VTI(k  )) / RHO(k)
        ENDIF
   
        DZ1 = ZM(K) - ZM(K-1)
   
        !    print*,k, qht(k) , DFHZ / DZ1,qit(k) , DFIZ / DZ1
        !     print*,k,VTI(k),VTH(k)
        qht(k) = qht(k) - DFHZ / DZ1 !hydrometeors don't
        qit(k) = qit(k) - DFIZ / DZ1  !nor does ice? hail, what about
      enddo
    end subroutine fallpart


    !******************************************************************************************************

    subroutine entrainment_coeff(entrainment_model, wind_eff, type_var,W_treshold,hBl,n,  &
                                & radius,radius_mdt,zm,zt,w,rho,rho_e, vel_p, vel_e,  & 
                                & dt_mdt,k, rbuoy, &
                                & epsi, delta, epsi_wind, delta_wind, &
                                & C_epsi, C_delta, C_wind,C_delta_wind)
      implicit none

      integer,intent(in) :: entrainment_model,k,n
      logical, intent(in) :: wind_eff
      character(len=1),intent(in) :: type_var
      real,intent(in)    :: W_treshold, hBl,dt_mdt, C_epsi, C_delta, C_wind,C_delta_wind
      real,dimension(:),intent(in) :: RADIUS,radius_mdt,zm,zt,w,rho,rho_e,vel_p,vel_e, rbuoy
      real,intent(out) :: epsi, delta,epsi_wind,delta_wind

      !local variable
      real :: tempkm,tempk,pi,temp_radius
      real :: lambda, surf_fire, l_surf_fire
      integer :: i,ii,kk
      real :: delta2, delat3,delat_sc, sign_w,WBAR,rhoBAR
      real :: zkp05,zk,zkm05
      real :: tempkm05, tempkp05,dwdz,radiusBAR,radiusBAR_mdt,dudz
      real :: u_h,location_mdt
      real :: radiuskp05,radiuskm05,sckp05,sckm05,rhokp05,rhokm05,drhoR2dz
      real :: wkup,wkdown, zkup,zkdown,ukup,ukdown,VELBAR_e,VELBAR_p,rhoBAR_e
      real :: buoyBAR
      !real :: C_epsi, C_delta, C_wind,C_delta_wind
      real :: swithch_dynEnt,x

      !C_epsi = 0.55    * 10 
      !C_delta = -10.e0 * 5
      !C_wind = 0.5e0
      !C_delta_wind = 1

      !lambda = 30.e0
      !surf_fire = pi * Radius(1)**2
      !l_surf_fire = sqrt(surf_fire)

      select case (type_var)
      case ('S')
          radiusBAR = RADIUS (k)
          WBAR      = 0.5 *(w(k-1)+w(k))
          if ( (k .gt. 6) .and. (k.lt.n-6)) then 
            wkup = 0.e0; zkup = 0.e0; ukup = 0.e0
            do ii=k,k+5
                  wkup   = wkup + 0.5*(w(ii)+w(ii-1))/6.e0
                  ukup   = ukup + vel_p(ii)/6.e0
                  zkup   = zkup  + zt(ii)/6.e0
            enddo
            wkdown = 0.e0; zkdown = 0.e0; ukdown = 0.e0
            do ii=k-5,k
                  wkdown = wkdown + 0.5*(w(ii)+w(ii-1))/6.e0
                  ukdown = ukdown + vel_p(ii)/6.e0
                  zkdown = zkdown + zt(ii)/6.e0
            enddo
            dwdz = (wkup-wkdown) / (zkup-zkdown)
            dudz = (ukup-ukdown) / (zkup-zkdown)
          else
            dwdz = (w(k)-w(k-1)) / (zm(k)-zm(k-1))          ! here it is the correct one
            dudz = (vel_p(k)-vel_p(k-1)) / (zt(k)-zt(k-1))  ! derivative the half level below
          endif
          zkp05 = zm(k)
          zk    = zt(k)
          zkm05 = zm(k-1)

          rhoBAR   = rho(k)
          rhoBAR_e = rho_e(k)

          VELBAR_e = VEL_e(k)
          VELBAR_p = VEL_p(k)

      !    tempkm05  =  0.5*(radius(k)+radius(k-1)) * w(k-1) * sqrt(lambda * zm(k-1)) * 0.5*(rho(k)+rho(k-1))
      !    tempkp05  =  0.5*(radius(k)+radius(k+1)) * w(k)   * sqrt(lambda * zm(k)  ) * 0.5*(rho(k)+rho(k+1))

          buoyBAR = rbuoy(k)

      case ('W')
          radiusBAR = 0.5* (RADIUS (k) + RADIUS (k+1) )
          wBAR = w(k)
          if ( (k .gt. 6) .and. (k.lt.n-6) )then 
              wkup = 0.e0; zkup = 0.e0; ukup = 0.e0
              do ii=k,k+5
                  wkup   = wkup + w(ii)/6.e0
                  ukup   = ukup + 0.5*(vel_p(ii)+vel_p(ii+1))/6.e0
                  zkup   = zkup + zm(ii)/6.e0
              enddo
              wkdown = 0.e0; zkdown = 0.e0; ukdown = 0.e0
              do ii=k-5,k
                  wkdown = wkdown + w(ii)/6.e0
                  ukdown = ukdown + 0.5*(vel_p(ii)+vel_p(ii+1))/6.e0
                  zkdown = zkdown + zm(ii)/6.e0
              enddo
              dwdz = (wkup-wkdown) / (zkup-zkdown)
              dudz = (ukup-ukdown) / (zkup-zkdown)
          else
              dwdz = (w(k)-w(k-1)) / (zm(k)-zm(k-1))            ! derivative the half level below
              dudz = (vel_p(k)-vel_p(k-1)) / (zt(k)-zt(k-1))    ! here it is the correct one
          endif
          zkp05 = zt(k+1)
          zk    = zm(k)
          zkm05 = zt(k)

          rhoBAR   = 0.5*(rho(k)+rho(k+1))
          rhoBAR_e = 0.5*(rho_e(k)+rho_e(k+1))

          VELBAR_e = 0.5 * (VEL_e(k)+VEL_e(k+1))
          VELBAR_p = 0.5 * (VEL_p(k)+VEL_p(k+1))

      !    tempkm05 =  radius(k)   *  0.5*(w(k)+w(k-1))  * sqrt(lambda * zt(k)   ) * rho(k)
      !    tempkp05 =  radius(k+1) *  0.5*(w(k)+w(k+1))  * sqrt(lambda * zt(k+1) ) * rho(k+1)

          buoyBAR = 0.5*(rbuoy(k)+rbuoy(k-1))
      end select

      sign_w = sign(1.e0, dwdz )
      !
      ! entrainement coefficient
      !
      if ( WBAR .gt. W_treshold  ) then 
          epsi = max(0.e0, C_epsi * buoyBAR / WBAR)
      else
          epsi = 8.e-4
      endif
      !
      ! detrainement coefficient
      !
      if ( WBAR .gt. W_treshold ) then  ! we are in the plume
          delta = max(0.e0, C_delta * buoyBAR / WBAR)
      else                            ! no effect of the plume and we need to detrain the mass ! we are oustside the plume
          delta = 2.e-2     !This is a high value, it detrain from the plume
                            !everything as at this level the vertical velocity is weak                
      endif
      !
      ! dynamical entrainment
      !
      epsi_wind = 0.e0
      delta_wind = 0.e0

      !if(zk .gt. 500) then 
      if (wind_eff) then 
        !if (RADIUSBAR .gt. 1.e0 ) then 
        !    DMDTM_dyn =  (2./pi)* ( rhoBAR_e/rhoBAR * VELBAR_e - VELBAR_p ) /RADIUSBAR
        !else
        !    DMDTM_dyn = 0.e0
        !endif

        !DMDTM_dyn = DMDTM_dyn * ALPHA_dyn 

        ! detrainment due to the windshear
        if ( ( abs(WBAR - VELBAR_p)/abs(VELBAR_p) .le. 1.5 ) .and. (WBAR .gt. W_treshold)  ) then
          !if ( WBAR .gt. W_treshold ) then
          x = abs(WBAR - VELBAR_p)/abs(VELBAR_p)
          epsi_wind =  C_wind * max(0.e0,dudz) * &
                              & (atan(-(x-1.2)/.08) - (atan(-(0.3)/.08))) / &
                              & (atan(-(-1.2)/.08) - (atan(-(0.3)/.08))) 
          delta_wind = C_delta_wind * epsi_wind 
        endif
      endif
      !endif

    end subroutine entrainment_coeff
      

    
    !***************************************************************************
    
    real function  esat_hPa(t)  
      ! saturation water vapour pressure hPa
      implicit none
      real :: t
      !     esat_hPa(millibars),t(kelvin)
      real :: tc
      tc = t - t00
      esat_hPa = 6.1078 * exp((17.2693882 * tc) / (tc + tmelt))
      return
    end     
    
    !***************************************************************************
    
    SUBROUTINE LBOUND_MTT (heat_fluxW, time_loc, dat, ifDebugPrint, ifStoreDump, uDump)
      !
      ! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
      !
      ! source of equations: J.S. Turner Buoyancy Effects in Fluids
      !                      Cambridge U.P. 1973 p.172,
      !                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
      !                      TID-25075, 1969, P.28
      !
      ! fundamentally a point source below ground. at surface, this produces
      ! a velocity w(1) and temperature T(1) which vary with time. There is
      ! also a water load which will first saturate, then remainder go into
      ! QC(1).
      ! EFLUX = energy flux at ground,watt/m**2 for the last DT
      !
      implicit none
      
      real,intent(in) :: time_loc
      real,intent(in) :: heat_fluxW
      type(TPRM_data), intent(inout) :: dat
      logical, intent(in) :: ifDebugPrint, ifStoreDump
      integer, intent(in) :: uDump

      ! local variables
      real :: es,  eflux, water,  pres, c1,  c2, f, zv,  denscor, xwater ,G_redduced_ov_g
      real :: Delta_T,temp
      !            
      dat%QH (1) = dat%QH (2)   !soak up hydrometeors
      dat%QI (1) = dat%QI (2)              
      dat%QC (1) = 0.       !no cloud here
      !
      CALL BURN(EFLUX, WATER, time_loc, heat_fluxW, dat%DT, dat%FMOIST, dat%CHF2TOTALH)  
      !
      !  calculate parameters at boundary from a virtual buoyancy point source
      !
      PRES = dat%PE (1) * 1000.   !need pressure in N/m**2

      ! original initialisation

      C1 = 5. / (6. * dat%ALPHA)  !alpha is entrainment constant
      C2 = 0.9 * dat%ALPHA  
      F = EFLUX / (PRES * CP * PI)  
      F = G * R * F * dat%AREA  !buoyancy flux
      ZV = C1 * dat%RSURF  !virtual boundary height
      if(abs(C1 * ( (C2 * F) **E1)) < 1e-10)then
        dat%W(1) = 0
      else
        dat%W(1) = C1 * ( (C2 * F) **E1) / ZV**E1  !boundary velocity
      endif

      !DENSCOR = C1 * F / G / (C2 * F) **E1 / ZV**E2   !density correction
      !T (1) = TE (1) / (1. - DENSCOR)    !temperature of virtual plume at zsurf

      G_redduced_ov_g = C1 * (F / G) * (C2 * F) **(-E1) * ZV**(-E2)   !density correction
      dat%T(1) = dat%TE(1) * (1. + G_redduced_ov_g)    !temperature of virtual plume at zsurf

         !write(*,'(a,3(3x,f20.6))'), 'boundary condition', EFLUX, W(1),T(1)-Te(1)
   
      !-- para testes
      !   W(1)=-W(1)/100.
      !   T(1)=TE(1)+13.5
      !

      dat%WC(1) = dat%W(1)
      dat%VEL_P(1) = 0.
            !rad_p(2) = rsurf**2 * rho(2) * sc(1)

         !SC(1) = max(sc(1)+  1.e0 * 1.e0/heat_fluxW *heat_fluxW * ( atan(-1.e0*(time_loc-1000.e0)/20.e0)  & 
         !                                                   - atan(1000.e0/20.e0) ), 1.e-4)

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !     match dw/dz,dt/dz at the boundary. F is conserved.
      !
         !WBAR = W (1) * (1. - 1. / (6. * ZV) )  
         !ADVW = WBAR * W (1) / (3. * ZV)  
         !ADVT = WBAR * (5. / (3. * ZV) ) * (DENSCOR / (1. - DENSCOR) )  
         !ADVC = 0.  
         !ADVH = 0.  
         !ADVI = 0.  
         !ADIABAT = - WBAR * G / CP  
      
      dat%VTH(1) = - 4.  
      dat%VTI(1) = - 3.  
      dat%TXS(1) = dat%T(1) - dat%TE(1)  
      dat%VISC(1) = dat%VISCOSITY
      dat%RHO(1) = 3483.8 * dat%PE(1) / dat%T(1)   !air density at level 1, g/m**3
      XWATER = WATER / (dat%W(1) * dat%DT * dat%RHO(1) )   !firewater mixing ratio
      dat%QV(1) = XWATER + dat%QVENV(1)  !plus what's already there 

      !  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
      ES = 0.1*ESAT_hPa(dat%T(1)) !blob saturation vapor pressure, em kPa
      !  rotina do plumegen ja calcula em kPa
      !  ES       = ESAT_PR (T(1))  !blob saturation vapor pressure, em kPa

      dat%EST(1)  = ES                                  
      dat%QSAT(1) = (EPS * ES) / (dat%PE(1) - ES)   !blob saturation lwc g/g dry air
  
      IF (dat%QV(1) > dat%QSAT(1) ) THEN  
        dat%QC(1) = dat%QV(1) - dat%QSAT(1) + dat%QC(1)  !remainder goes into cloud drops
        dat%QV(1) = dat%QSAT(1)
      ENDIF  
      !
      CALL WATERBAL(dat%T(1), dat%RHO(1), dat%QC(1), dat%QH(1), dat%QI(1), dat%QV(1), &
                  & dat%QSAT(1), dat%WBAR, dat%DT, 0.0, dat%EST(1), dat%CVI(1))
      !
      dat%radius(1) = dat%rsurf
      dat%rad_p(1) = dat%rsurf**2 * dat%rho(1) !* sc(1)
      dat%sc(1) = 1.e0
      !SC(1) = eflux/heat_fluxW * 1.e0/rho(1) !SCE(1)+F/1000.*dt  ! gas/particle (g/g)   
      dat%mass_sc(1) = dat%rsurf**2 * dat%rho(1) * dat%sc(1)
      !  
      if (ifStoreDump)  write(uDump,'(A15,4(g12.5,3x))') 'EFLUX, w, sc:', time_loc, EFLUX, dat%W(1), dat%sc(1)

      temp = dat%rho(1) * (cp/1.e3) * (dat%T(1)-dat%TE(1)) * dat%w(1) !/ (radius(1)**2 * 3.14)
      
      if(ifDebugPrint) write(*,'("CHF(cpDeltaTw)= "f15.5,3x,"FRP*FRP2CHF= "f15.3,3x,"w0= "f15.3)') temp*1.e-6, &
         & eflux*1.e-6 , dat%w(1)
  
      RETURN  
    END SUBROUTINE LBOUND_MTT
      
    
    !*****************************************************************

    SUBROUTINE BURN(EFLUX, WATER, time_loc, heat_fluxW, DT, FMOIST, CHF2TOTALH)  
      !
      !- calculates the energy flux and water content at lboundary
      !
      implicit none
    
      real, parameter :: HEAT = 18.7E6 !Joules/kg - forest in Alta Floresta (MT)

      real, intent(in) :: time_loc, heat_fluxW, DT, FMOIST, CHF2TOTALH
      real, intent(out) :: EFLUX, WATER

      !
      ! The emission factor for water is 0.5. The water produced, in kg,
      ! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
      ! The fire burns for DT out of TDUR seconds, the total amount of
      ! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
      ! considered to be spread over area AREA and so the mass burned per
      ! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.

      EFLUX = max ( heat_fluxW * ( atan((time_loc-100.e0)/20.e0) - atan(-100.e0/20.e0) ), 1.e-4) ! Watts/m**2 

      !  WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST)       ! kg/m**2 
      !
      WATER = 1000. * EFLUX * (DT / HEAT) * (0.5 + FMOIST) / CHF2TOTALH ! g/m**2 
      !
      !   WATER = EFLUX*1.e-6 * .368 * (0.5 + FMOIST)
      RETURN  
    END SUBROUTINE BURN 


    !************************************************************************

    SUBROUTINE WATERBAL(T, RHO, QC, QH, QI, QV, QSAT, WBAR, DT, DQSDZ, EST, CVI)
      !
      ! Describes the evolution of water in liquid and gaseous phases
      ! Everything here is formulated for one layer
      !
      implicit none
          
      ! imported parameters
      real, intent(in) :: RHO, QSAT, WBAR, DT, DQSDZ,EST, CVI
      real, intent(inout) :: T, QC, QH, QI, QV
          
      IF (QC <= 1.0E-10) QC = 0.  !DEFEAT UNDERFLOW PROBLEM
      IF (QH <= 1.0E-10) QH = 0.  
      IF (QI <= 1.0E-10) QI = 0.  

      CALL EVAPORATE(T, RHO, QC, QH, QI, QV, QSAT, WBAR, DT, DQSDZ, EST, CVI)    !vapor to cloud,cloud to vapor  

      CALL SUBLIMATE(T, RHO, QSAT, QV)    !vapor to ice  

      CALL GLACIATE(T, RHO, QSAT, QV, QH)     !rain to ice 

      CALL MELT(T, RHO, QI, QH)         !ice to rain

      CALL CONVERT (T, RHO, QC, QH) !(auto)conversion and accretion 

      RETURN  

      CONTAINS
        
        !----------------------------------------------------------------
        
        SUBROUTINE EVAPORATE(T, RHO, QC, QH, QI, QV, QSAT, WBAR, DT, DQSDZ, EST, CVI)
          !
          !- evaporates cloud,rain and ice to saturation
          !
          implicit none
          !
          !     XNO=10.0E06
          !     HERC = 1.93*1.E-6*XN035        !evaporation constant
          !
          real, intent(in) :: RHO, QSAT, WBAR, DT, DQSDZ, EST, CVI
          real, intent(inout) :: T, QC, QH, QI, QV

          real, PARAMETER :: HERC = 5.44E-4
          
          ! Local variables
          real :: evhdt, evidt, evrate, evap, sd,	quant, dividend, divisor, devidt

          SD = QSAT - QV  !vapor deficit
          
          !IF (SD == 0.0)  RETURN      ! mas: such comparison for computed reals makes no sense
          IF (abs(SD) < 1.e-7)  RETURN  
          EVHDT = 0.  
          EVIDT = 0.  
          !evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0
          EVRATE = ABS (WBAR * DQSDZ)   !evaporation rate (Kessler 8.32)
          EVAP = EVRATE * DT   !what we can get in DT

          IF (SD <= 0.0) THEN  !     condense. SD is negative
             IF (EVAP >= ABS(SD) ) THEN    !we get it all
                QC = QC - SD  !deficit,remember?
                QV = QSAT       !set the vapor to saturation  
                T = T - SD * FRC  !heat gained through condensation per gram of dry air
                RETURN  
             ELSE  
                QC = QC + EVAP         !get what we can in DT 
                QV = QV - EVAP         !remove it from the vapor
                T = T + EVAP * FRC   !get some heat
                RETURN  
             ENDIF  
          ELSE                                !SD is positive, need some water
          !
          ! not saturated. saturate if possible. use everything in order
          ! cloud, rain, ice. SD is positive
            IF (EVAP.LE.QC) THEN
              !enough cloud to last DT  
              IF (SD.LE.EVAP) THEN          !enough time to saturate
                   QC = QC - SD       !remove cloud                                          
                   QV = QSAT          !saturate
                   T = T - SD * FRC   !cool the parcel                                          
                   RETURN  !done
              ELSE   !not enough time
                   SD = SD-EVAP               !use what there is
                   QV = QV + EVAP     !add vapor
                   T = T - EVAP * FRC !lose heat
                   QC = QC - EVAP     !lose cloud go on to rain.                                      
              ENDIF
            ELSE
              !not enough cloud to last DT
              IF (SD.LE.QC) THEN   !but there is enough to sat
                   QV = QSAT  !use it
                   QC = QC - SD  
                   T = T - SD * FRC  
                   RETURN  
                ELSE            !not enough to sat
                   SD = SD - QC
                   QV = QV + QC
                   T = T - QC * FRC         
                   QC = 0.0  !all gone
                ENDIF       !on to rain                          
            ENDIF          !finished with cloud
            !
            !  but still not saturated, so try to use some rain
            !  this is tricky, because we only have time DT to evaporate. if there
            !  is enough rain, we can evaporate it for dt. ice can also sublimate
            !  at the same time. there is a compromise here.....use rain first, then
            !  ice. saturation may not be possible in one DT time.
            !  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
            !  sd is still positive or we wouldn't be here.
            !
            IF (QH > 1.E-10)then
              !             rain evaporation in time DT
              !srf-25082005
              !  QUANT = (QC + QV - QSAT ) * RHO  !g/m**3
              QUANT = ( QSAT - QC - QV ) * RHO   !g/m**3
              EVHDT = (DT * HERC * (QUANT) * (QH * RHO) **.65) / RHO
              IF (EVHDT.LE.QH) THEN     !enough rain to last DT
                IF (SD.LE.EVHDT) THEN   !enough time to saturate	  
                   QH = QH - SD   !remove rain	  
                   QV = QSAT      !saturate	  
                   T = T - SD * FRC  !cool the parcel		  
                   !if(mintime>40) print*,'1',L,T(L)-273.15,QV(L)*1000,QH(L)*1000
                   RETURN
                ELSE                              !not enough time
                   SD = SD-EVHDT             !use what there is
                   QV = QV + EVHDT      !add vapor
                   T = T - EVHDT * FRC  !lose heat
                   QH = QH - EVHDT      !lose rain
                   !if(mintime>40.and. L<40) print*,'2',L,T(L)-273.15,QV(L)*1000,QH(L)*1000
                   !if(mintime>40.and. L<40) print*,'3',L,EVHDT,QUANT
                   !if(mintime>40.and. L<40) print*,'4',L,QC (L)*1000. , QV (L)*1000. , QSAT (L)*1000.
                ENDIF  				  !go on to ice.
              ELSE  
                !not enough rain to last DT
                IF (SD <= QH) THEN             !but there is enough to sat
                   QV = QSAT                !use it
                   QH = QH - SD  
                   T = T - SD * FRC  
                   RETURN  
                ELSE                              !not enough to sat
                   SD = SD-QH
                   QV = QV + QH
                   T = T - QH * FRC    
                   QH = 0.0                   !all gone
                ENDIF                             !on to ice
              ENDIF                                !finished with rain
            endif  ! QH(L) exists
            !
            !  now for ice
            !  equation from (OT); correction factors for units applied
            !
            IF (QI <= 1.E-10) RETURN             !no ice there
            DIVIDEND = ( (1.E6 / RHO) **0.475) * (SD / QSAT - 1) * (QI**0.525) * 1.13
            DIVISOR = 7.E5 + 4.1E6 / (10. * EST)  
            DEVIDT = - CVI * DIVIDEND / DIVISOR  !rate of change
            EVIDT = DEVIDT * DT                  !what we could get
            !
            ! logic here is identical to rain. could get fancy and make subroutine
            ! but duplication of code is easier. God bless the screen editor.
            !
            IF (EVIDT.LE.QI) THEN             !enough ice to last DT
                IF (SD.LE.EVIDT) THEN  		  !enough time to saturate
                   QI = QI - SD     !remove ice
                   QV = QSAT        !saturate
                   T = T - SD * SRC   !cool the parcel
                   RETURN           !done
                ELSE                      !not enough time
                   SD = SD-EVIDT  		  !use what there is
                   QV = QV + EVIDT  	  !add vapor
                    T =  T - EVIDT * SRC  !lose heat
                   QI = QI - EVIDT     !lose ice
                ENDIF  				   !go on,unsatisfied
             ELSE                      !not enough ice to last DT
                IF (SD.LE.QI) THEN     !but there is enough to sat
                   QV = QSAT           !use it
                   QI = QI - SD  
                    T =  T - SD * SRC  
                   RETURN  
                ELSE                   !not enough to sat
                   SD = SD-QI
                   QV = QV + QI
                   T = T - QI * SRC             
                   QI = 0.0            !all gone
                ENDIF                  !finished with ice
             ENDIF  
          ENDIF                        !finished with the SD decision
          RETURN  
        END SUBROUTINE EVAPORATE

        !---------------------------------------------------------------
        
        SUBROUTINE CONVERT(T, RHO, QC, QH) 
          !
          !- ACCRETION AND AUTOCONVERSION. Formulated for the single layer
          !
          implicit none
          
          real, intent(in) :: T, RHO
          real, intent(inout) :: QC, QH
          
          ! Local parameters
          real,      PARAMETER :: AK1 = 0.001    !conversion rate constant
          real,      PARAMETER :: AK2 = 0.0052   !collection (accretion) rate
          real,      PARAMETER :: TH  = 0.5      !Kessler threshold
          integer,   PARAMETER :: iconv = 1      !Kessler conversion
          !real, parameter :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime)
          real, parameter :: ANBASE =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
          ! na formulacao abaixo use o valor em #/cm^3  
          !real, parameter :: BDISP = 0.366       !Berry--size dispersion (maritime)
          real, parameter :: BDISP = 0.146       !Berry--size dispersion (continental)

          ! Local variables
          real ::   accrete, con, q, h, bc1,   bc2,  total

          !     selection rules
          IF (T <= TFREEZE) RETURN  !process not allowed above ice
          IF (abs(QC) < 1.e-7) RETURN  

          ACCRETE = 0.  
          CON = 0.  
          Q = RHO * QC
          H = RHO * QH

          IF (QH > 0.) ACCRETE = AK2 * Q * (H**.875)  !accretion, Kessler
          
          IF (ICONV.NE.0) THEN   !select Berry or Kessler
            !old   BC1 = 120.  
            !old   BC2 = .0266 * ANBASE * 60.  
            !old   CON = BDISP * Q * Q * Q / (BC1 * Q * BDISP + BC2) 	  
            CON = Q*Q*Q*BDISP/(60.*(5.*Q*BDISP+0.0366*ANBASE))
          ELSE  
            !   CON = AK1 * (Q - TH)   !Kessler autoconversion rate
            !   IF (CON.LT.0.0) CON = 0.0   !havent reached threshold
            CON = max(0.,AK1 * (Q - TH)) ! versao otimizada
          ENDIF  

          TOTAL = (CON + ACCRETE) * DT / RHO

          IF (TOTAL.LT.QC) THEN  
            QC = QC - TOTAL  
            QH = QH + TOTAL    !no phase change involved
          ELSE  
            QH = QH + QC    !uses all there is
            QC = 0.0  
          ENDIF  
          RETURN  
        END SUBROUTINE CONVERT

        !----------------------------------------------------------------------

        SUBROUTINE CONVERT2 (T, RHO, W, QC, QH, DZM) 
          implicit none

          real, intent(in) :: T, RHO, W, DZM
          real, intent(inout) :: QC, QH
          
          ! local variables
          LOGICAL, parameter :: AEROSOL = .true.
          real, parameter :: LAT=2.5008E6, DB=1., NB=1500. !ALPHA=0.2 
          real :: KA,KEINS,KZWEI,KDREI,VT	
          real :: CON,ACCRETE,total 
          real :: RHO_1e3
      
          RHO_1e3 =  RHO * 1.e-3 ! dens (MKS) ??

          ! autoconversion
          KA = 0.0005 
          IF( T .LT. 258.15 )THEN
          !   KEINS=0.00075
              KEINS=0.0009 
              KZWEI=0.0052
              KDREI=15.39
          ELSE
              KEINS=0.0015
              KZWEI=0.00696
              KDREI=11.58
          ENDIF
    
          !   RHO_1e3=PE/RD/TE
          VT=-KDREI* (QH/RHO_1e3)**0.125

          IF (W.GT.0.0 ) THEN
           IF (AEROSOL) THEN
             CON = 1/W  *  QC*QC*1000./( 60. *( 5. + 0.0366*NB/(QC*1000.*DB) )  )
           ELSE
             IF (QC.GT.(KA*RHO_1e3)) THEN
             !print*,'1',QC,KA*RHO_1e3
             CON = KEINS/W *(QC - KA*RHO_1e3 )
             ENDIF
           ENDIF
          ELSE
             CON = 0.0
          ENDIF

          ! accretion
          IF(W.GT.0.0) THEN
             ACCRETE = KZWEI/(W - VT) * MAX(0.,QC) *   &
                      MAX(0.001,RHO_1e3)**(-0.875)*(MAX(0.,QH))**(0.875)
          ELSE
             ACCRETE = 0.0
          ENDIF

          TOTAL = (CON + ACCRETE) * (1 / DZM)    ! DT / RHO (L)  

          IF (TOTAL.LT.QC) THEN  
             QC = QC - TOTAL  
             QH = QH + TOTAL    !no phase change involved
          ELSE  
             QH = QH + QC    !uses all there is
             QC = 0.0  
          ENDIF  
          RETURN  
        END SUBROUTINE CONVERT2
        ! ice - effect on temperature
        !      TTD = 0.0 
        !      TTE = 0.0  
        !       CALL ICE(QSATW,QSATE,T,QC,QH, &
        !               TTA,TTB,TTC,DZ,RHO_1e3,D,C,TTD,TTE)
        !       DYDX(1) = DYDX(1) + TTD  + TTE ! DT/DZ on Temp

        
        !--------------------------------------------------------------------------

        SUBROUTINE SUBLIMATE(T, RHO, QSAT, QV)
          !
          ! VAPOR TO ICE (USE EQUATION OT22)
          !
          implicit none

          real, intent(in) :: RHO, QSAT
          real, intent(inout) :: T, QV
          
          ! Local variables
          real :: dtsubh,  dividend, divisor, subl
          !
          DTSUBH = 0.  
          !
          !selection criteria for sublimation
          IF (T > TFREEZE  ) RETURN  
          IF (QV <= QSAT) RETURN  
          !
          !     from (OT); correction factors for units applied
          !
          DIVIDEND = ( (1.E6 / RHO) **0.475) * (QV / QSAT - 1) * (QI **0.525) * 1.13
          DIVISOR = 7.E5 + 4.1E6 / (10. * EST)  
          !
          DTSUBH = ABS (DIVIDEND / DIVISOR)   !sublimation rate
          SUBL = DTSUBH * DT                  !and amount possible
          !
          !     again check the possibilities
          !
          IF (SUBL < QV) THEN  
            QV = QV - SUBL             !lose vapor
            QI = QI + SUBL  	      !gain ice
            T = T + SUBL * SRC         !energy change, warms air
            RETURN  
          ELSE  
             QI = QV                    !use what there is
             T = T + QV * SRC      !warm the air
             QV = 0.0  
          ENDIF  
          RETURN  
        END SUBROUTINE SUBLIMATE

        !----------------------------------------------------------

        SUBROUTINE GLACIATE(T, RHO, QSAT, QV, QH)
          !
          ! CONVERSION OF RAIN TO ICE
          !     uses equation OT 16, simplest. correction from W not applied, but
          !     vapor pressure differences are supplied.
          !
          implicit none
          !
          real, intent(in) :: RHO, QSAT
          real, intent(inout) :: T, QV, QH

          ! Local variables
          real, PARAMETER :: GLCONST = 0.025   !glaciation time constant, 1/sec
          real :: dfrzh

          DFRZH = 0.    !rate of mass gain in ice

          !selection rules for glaciation
          IF (QH <= 0.) RETURN
          IF (QV < QSAT) RETURN
          IF (T > TFREEZE) RETURN
          !
          !      NT=TMELT-T(L)
          !      IF (NT.GT.50) NT=50
          !
          DFRZH = DT * GLCONST * QH   ! from OT(16)

          IF (DFRZH < QH) THEN  
             QI = QI + DFRZH  
             QH = QH - DFRZH  
             T = T + FRC * DFRZH  !warms air
          ELSE  
             QI = QI + QH
             T = T + FRC * QH
             QH = 0.0  
          ENDIF  
          RETURN  
        END SUBROUTINE GLACIATE

        !----------------------------------------------------------------------

        SUBROUTINE MELT(T, RHO, QI, QH)
          !
          ! MAKES WATER OUT OF ICE
          !
          implicit none
          !                                              
          real, intent(in) :: RHO
          real, intent(inout) :: T, QI, QH
          
          real, PARAMETER :: FRC = 332.27, F0 = 0.75   !ice velocity factor
          
          ! Local variables
          real DTMELT
          
          DTMELT = 0.   !conversion,ice to rain
          !
          !selection rules
          IF (QI <= 0.0  ) RETURN  
          IF (T  < TMELT) RETURN  
                                                                          !OT(23,24)
          DTMELT = DT * (2.27 / RHO) * CVI * (T - TMELT) * & 
                   ( (RHO * QI * 1.E-6)**0.525) * (F0**(- 0.42) )  !after Mason,1956
          !
          ! Anything left in ice?
          !
          IF (DTMELT < QI) THEN  
             QH = QH + DTMELT  
             QI = QI - DTMELT  
             T = T - FRC * DTMELT     !cools air
          ELSE  
             QH = QH + QI   !get all there is to get
             T = T - FRC * QI
             QI = 0.0  
          ENDIF  
          RETURN  
        END SUBROUTINE MELT
        
    END SUBROUTINE WATERBAL

    
    !***************************************************************************
    REAL FUNCTION ESAT_PR (TEM)  
      !
      ! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
      !
      implicit none
          
      ! local variabes
      real, PARAMETER :: CI1 = 6.1115, CI2 = 22.542, CI3 = 273.48
      real, PARAMETER :: CW1 = 6.1121, CW2 = 18.729, CW3 = 257.87, CW4 = 227.3

      real , intent(in) :: tem
      real :: temc
      !
      !     formulae from Buck, A.L., JAM 20,1527-1532
      !     custom takes esat wrt water always. formula for h2o only good to -40C so:
      !
      TEMC = TEM - TMELT  

      IF (TEMC <= -40.0)then
        ESAT_PR = CI1 * EXP (CI2 * TEMC / (TEMC + CI3) ) /10.  !ice, kPa
      else
        ESAT_PR = CW1 * EXP ( ( (CW2 - (TEMC / CW4) ) * TEMC) / (TEMC + CW3)) !kPa			  
      endif
      RETURN  
    END function ESAT_PR

!********************************************************
    Subroutine min_max_1d(sv, minsv, maxsv, ifPrintDebug)

      Implicit none

      ! Imported parameters
      real, dimension(:), intent(in) :: sv
      logical, intent(in) :: ifPrintDebug
      Real, intent(out) :: minsv, maxsv
      
      ! Local variables
      Integer :: i,im
 
      maxsv = -1.e10
      im = -999
      do i=1,size(sv)
         if(sv(i) .gt. maxsv) then 
            maxsv = sv(i)
            im = i
         endif
      enddo   
      if(ifPrintDebug) write(*,*) 'max = ',maxsv, 'at ', im

      minsv = 1.e10
      im = -999
      do i=1,size(sv)-2
         if(sv(i) .lt. minsv) then 
            minsv = sv(i)
            im = i
         endif
      enddo
      if(ifPrintDebug) print *, 'min = ',minsv, 'at ', im
      
      Return
    End subroutine min_max_1d    
      
      
    !***************************************************************************

    subroutine printout (dat, nrectotal, time, uDump)  
      !
      ! Prints the debug information 
      !
      implicit none

      ! imported parameters
      type(TPRM_data), intent(in) :: dat
      integer, intent(in) :: nrectotal, uDump
      real, intent(in) :: TIME
      
      ! Local variables
      integer, save :: nrec=0, nrecx=0
      integer, parameter :: interval = 1              !debug time interval,min
      integer :: ko, ii, k_initial, k_final, KK4, kl, ix, irange, k
      real :: pea, btmp,etmp,vap1,vap2,gpkc,gpkh,gpki,deficit,mass_cum, w_thresold
      real, dimension(nzPRM) :: mass, mass_fract

      Real,dimension(:),allocatable :: w_zt
      Integer*4 :: timeStr
      character(len=50) :: pltfile

      if (.false.) then     ! change if want vertical distribution
          !- vertical mass distribution
          k_initial= 0
          k_final  = 0
          w_thresold = 1.

          !- define range of the upper detrainemnt layer
          mass_cum=0.
          do ko=nzPRM-10,2,-1
              if(dat%w(ko) < w_thresold) cycle
              if(k_final==0) k_final=ko
              if(dat%w(ko)-1. > dat%w(ko-1)) then
                  k_initial=ko
                  exit
              endif
              mass_cum = mass_cum + mass(ko)
          enddo
          print*,'ki kf, mass_cum =',k_initial,k_final, mass_cum

          !- if there is a non zero depth layer, make the mass vertical distribution 
          mass_fract=0.
          if(k_final > 0 .and. k_initial > 0 ) then 
              k_initial=int((k_final+k_initial)*0.5)
              !- get the normalized mass distribution
              do ko=nzPRM-10,1,-1
                  if(dat%w(ko) < w_thresold) cycle
                  if(dat%w(ko)-1.0 > dat%w(ko-1)) exit
                  mass_fract(ko) = mass(ko) / mass_cum
              enddo
              !v-2.5
              !- parabolic vertical distribution between k_initial and k_final
              mass=0.
              KK4 = k_final-k_initial+2
              do ko=1,kk4-1
                  kl=ko+k_initial-1
                  mass(kl) = 6.* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
              enddo
              !print*,'ki kf=',int(k_initial+k_final)/2,k_final,sum(mass)*100.

              !- check if the integral is 1 (normalized)
              if(sum(mass) .ne. 1.) then
                  mass_cum = ( 1.- sum(mass) )/float(k_final-k_initial+1)
                  do ko=k_initial,k_final
                      mass(ko) = mass(ko)+ mass_cum !- values between 0 and 1.
                  enddo
                  ! print*,'new mass=',sum(mass)*100.,mass_cum
                  !pause
              endif
          endif !k_final > 0 .and. k_initial > 

          !- vertical mass distribution (horizontally integrated : kg/m)
          !mass=3.14*rad_p**2*dne*sc*100.
          !mass=mass/sum(mass)
      endif   ! if process vertical distribution
      !
      IF(dat%MINTIME == 1) nrec = 0
      !
      WRITE (uDump, fmt='(I5,A,f6.2,A,f8.2,A)')dat%MINTIME,' minutes,  dt= ', dat%DT, ' seconds  time= ',TIME,'seconds)'
      WRITE (uDump, fmt='(A,f10.2,A)') '     Ztop= ',dat%ZTOP, ' meters'   
      WRITE (uDump, '(A,A)') &
          & '   Z(km)      P(hPa)      W(m/s)      T(C)        T-TE       QV(g/kg)     SAT(g/kg)    QC(g/kg)',&
          &'    QH(g/kg)    QI(g/kg)  Buoy(1e-2 m/s2)'
      !
      DO KO = 1, 70, interval  
          PEA = dat%PE(KO) * 10.       !pressure is stored in decibars(kPa),print in mb;
          BTMP = dat%T(KO) - TMELT     !temps in Celsius
          ETMP = dat%T(KO) - dat%TE (KO)   !temperature excess
          VAP1 = dat%QV(KO)   * 1000.  !printout in g/kg for all water,
          VAP2 = dat%QSAT(KO) * 1000.  !vapor (internal storage is in g/g)
          GPKC = dat%QC(KO)   * 1000.  !cloud water
          GPKH = dat%QH(KO)   * 1000.  !raindrops
          GPKI = dat%QI(KO)   * 1000.  !ice particles 
          DEFICIT = VAP2 - VAP1     !vapor deficit
          WRITE (uDump, '(11g12.5)') dat%zt(KO)/1000., PEA, dat%W(KO), BTMP, ETMP, VAP1, &
                                   & VAP2, GPKC, GPKH, GPKI, dat%rbuoy(KO)*100. !, VTH (KO), SC(KO)
      end do

      !nrec=nrec+1
      !write (19,rec=nrec) (W (KO), KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (T (KO), KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (TE(KO), KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (QV(KO)*1000., KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) ((QC(KO)+QI(ko))*1000., KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (QH(KO)*1000., KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (rbuoy(KO), KO=1,nrectotal)
      !!   write (19,rec=nrec) (QI(KO)*1000., KO=1,nrectotal)
      !nrec=nrec+1
      !!   write (19,rec=nrec) (dwdt_entr(KO), KO=1,nrectotal)
      !!   write (19,rec=nrec) (100.*QV(KO)/QSAT(KO), KO=1,nrectotal)
      !!   write (19,rec=nrec) (THEE(KO), KO=1,nrectotal)
      !!   write (19,rec=nrec) (radius(KO), KO=1,nrectotal)
      !nrec=nrec+1
      !   write (19,rec=nrec) (QVENV(KO)*1000., KO=1,nrectotal)
      !   write (19,rec=nrec) (THEQ(KO), KO=1,nrectotal)
      !   write (19,rec=nrec) (upe(KO), KO=1,nrectotal)
      !write (19,rec=nrec) (mass(kO), KO=1,nrectotal)
      !nrec=nrec+1
      !!   write (19,rec=nrec) (tde(KO), KO=1,nrectotal)
      !write (19,rec=nrec) (vel_e(KO), KO=1,nrectotal)
      !nrec=nrec+1
      !write (19,rec=nrec) (vel_p(KO), KO=1,nrectotal) ! Pa
      !!   write (19,rec=nrec) (rbuoy(KO), KO=1,nrectotal) ! Pa
      !!   write (19,rec=nrec) (dne(KO), KO=1,nrectotal) ! Pa
      !!   write (19,rec=nrec) (vt3dh(KO), KO=1,nrectotal) ! Pa
      !!   write (19,rec=nrec) (pe(KO)*1000., KO=1,nrectotal) ! Pa

      allocate(w_zt(nrectotal)) 
      w_zt(1) = dat%w(1)
      do k =2,nrectotal
          w_zt(k) = 0.5*(dat%w(k-1) + dat%w(k))
      enddo
      ! write txt file
      !---------------------
!      write(pltfile,'(a,i6.6,a,i2.2,a)') './plumegen_',int(time),'_',int(100*(time-int(time))),'.txt'
!      open(29,file=pltfile(1:28)) 
      write(uDump,'(3a)')'#       z(km)     w (m/s)         T(K)          T-Te (K)      QV(kg/kg)      QC(kg/kg)        QH(k)',&
                      '         rbuoy      epsi               delta         epsi_wind      delta_winf         sc           ',&
                      '  rho          radius         rho*radius**2      u_h            u_he'
      do k=1,nrectotal
          write(uDump,'(f12.6,3x,17(e12.5,3x))') dat%zt(k)/1.e3, w_zt(k), dat%T(k), dat%T(k)-dat%TE(k), &
                    & dat%QV(k), dat%QC(k), dat%QH(k), dat%rbuoy(k), &
                    & dat%dwdt_entr(1:4,k), & 
                    & max(dat%sc(k),1.e-10), dat%rho(k), dat%radius(k),dat%rad_p(k), dat%vel_p(k), dat%vel_e(k)
      enddo
!      close(29)
      deallocate(w_zt)
      RETURN  
    end subroutine printout

    
    ! *********************************************************************
     
    subroutine thetae(p, t, rv, the, tdd)
      implicit none
      real, intent(in) :: p,t,rv
      real, intent(out) :: the,tdd

      real :: pit,tupo,ttd,dz,tupn,tmn
      integer :: itter

      pit=p
      tupo=t
      ttd=td(p,rv)
      tdd=ttd
      dz=cp/g*(t-ttd)
      if(dz.le.0.) goto 20
      do itter=1,50
        tupn=t-g/cp*dz
        tmn=(tupn+t)*.5*(1.+.61*rv)
        pit=p*exp(-g*dz/(r*tmn))
        if(abs(tupn-tupo).lt.0.001) goto 20
        ttd=td(pit,rv)
        tupo=tupn
        dz=dz+cp/g*(tupn-ttd)
      enddo
      stop 'stop: problems with thetae calculation - RONAN'
      20 continue
      the=tupo*(1e5/pit)**.286*exp(alvl*rv/(cp*tupo))
      return
      
    CONTAINS
      real function td(p,rs)
        implicit none
        real :: rr,rs,es,esln,p
        rr=rs+1e-8
        es=p*rr/(.622+rr)
        esln=log(es)
        td=(35.86*esln-4947.2325)/(esln-23.6837)
        return
      end 
    end subroutine thetae


END MODULE plume_rise_PRMv1
