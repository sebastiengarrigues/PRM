!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module rconstants

!---------------------------------------------------------------------------
real, public, parameter ::            &
        rgas    = 287.               &
    ,   cp      = 1004.              &
    ,   cv      = 717.               &
    ,   rm      = 461.               &
    ,   p00     = 1.e5               &
    ,   t00     = 273.16             &
    ,   g       = 9.80796            &
    ,   pi      = 3.1415927          &
    ,   pi_180  = 3.1415927 / 180.   &
    ,   pi_4    = 3.1415927 / 4.     &
    ,   spcon   = 111120.            &
    ,   erad    = 6367000.           &
    ,   vonk    = 0.40               &
    ,   tkmin   = 5.e-4              &
    ,   alvl    = 2.50e6             &
    ,   alvi    = 2.834e6            &
    ,   alli    = 0.334e6            &
    ,   alvl2   = 6.25e12            &
    ,   alvi2   = 8.032e12           &
    ,   solar   = 1.3533e3           &
    ,   stefan  = 5.6696e-8          &
    ,   cww     = 4218.              &
    ,   c0      = 752.55 * 4.18684e4 &
    ,   viscos  = .15e-4             &
    ,   rowt    = 1.e3               &
    ,   dlat    = 111120.            &
    ,   omega   = 7.292e-5           &
    ,   rocp    = rgas / cp          &
    ,   p00i    = 1. / p00           &
    ,   cpor    = cp / rgas          &
    ,   rocv    = rgas / cv          &
    ,   cpi     = 1. / cp            &
    ,   cpi4    = 4. * cpi           &
    ,   cp253i  = cpi / 253.         & 
    ,   allii   = 1. / alli          &
    ,   aklv    = alvl / cp          &
    ,   akiv    = alvi / cp          &
    ,   gama    = cp / cv            &
    ,   gg      = .5 * g             &
    ,   ep      = rgas / rm          & 
    ,   p00k    = 26.870941          &  !  = p00 ** rocp  
    ,   p00ki   = 1. / p00k
!---------------------------------------------------------------------------


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
  !      XNO=10.0E06 median volume diameter raindrop (K table 4)
  !      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
  !
  real, public, parameter :: vc = 5.107387
  real, public, parameter :: r = 287.04, eps = 0.622,  tmelt = 273.3
  real, public, parameter :: heatsubl = 2.834e6, heatfus = 3.34e5, heatcond = 2.501e6
  real, public, parameter :: tfreeze = 269.3, e1 = 1./3., e2 = 5./3.
  real, public, PARAMETER :: SRC = HEATSUBL / CP, FRC = HEATFUS / CP

  REAL, PARAMETER, PUBLIC :: radians_to_degrees =  57.29577951
  REAL, PARAMETER, PUBLIC :: degrees_to_radians =  0.01745329252
  REAL, PARAMETER, PUBLIC :: earth_radius = 6378000.0
  REAL, PARAMETER, PUBLIC :: F_NAN = TRANSFER(2143289344,1.0)



end Module
