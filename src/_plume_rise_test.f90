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
  integer :: iTmp, status, ix, iy, iz, izInverse
  logical :: eof
  real :: a, b, fABL, pressure
  character(len=100) :: strTmp
  character(len=*), parameter :: frpgrib='../Input/z_cams_c_ecmf_202104260200_gfas_an_sfc_001_frpfire.grib'
  character(len=*), parameter :: plumerisegrib='../output/z_cams_c_ecmf_202104260200_gfas_an_sfc_001_PLUMERIZE.grib'
  character(len=100), dimension(3), parameter :: metgribs = [Character(len=100) ::  &
                                                            & "../Meteo/ecglob10deg_2021042600+002grib2.ml", &
                                                            & "../Meteo/ecglob10deg_2021042600+002.lnsp", &
                                                            & "../Meteo/ecglob10deg_2021042600+002.sfc" ]

  !
  ! Create and initialise a bunch of fires
  !

  call  acquire_fires(fires, frpgrib)
  call acquire_meteo(meteo_data, fires, metgribs)
  ! stop 1


  !
  ! The main call. Dump file will be replicated to many in case of OMP run
  !
  call compute_plume_rise(fires, &                ! fires, input and output
                        & meteo_data, &           ! meteo data, input
                        & '../output/Dump.log_')
  !
  ! Print the results
  !
  call store_plume_rise(fires, frpgrib, plumerisegrib) 
  
  
end program test_plume_rise
  
