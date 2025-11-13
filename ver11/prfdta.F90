      BLOCK DATA PRFDTA
!______________________________________________________________________|
!                                                                      |
!     Copyright 2009 Spectral Sciences, Inc.                           |
!                                                                      |
!     Neither Spectral Sciences, Inc. nor its employees makes any      |
!     warranty, express or implied, regarding the uses, accuracy,      |
!     or completeness of this software.                                |
!______________________________________________________________________|

!     PRFDTA CONTAINS AEROSOL PROFILE DATA.

!     COMMONS:

!     /PRFD/ (Altitudes below are defined assuming GNDALT is 0 km)
!       I_HGT   ALTITUDE GRID FOR AEROSOL PROFILE DATA [KM].
!       HZ2K    0-3KM HAZE W/ 50, 23, 10, 5 & 2 KM VIS [550NM EXT/KM].
!       FAWI50  2-11KM  FALL/WINTER  PROFILE W/ 50KM VIS [550NM EXT/KM].
!       FAWI23  2-11KM  FALL/WINTER  PROFILE W/ 23KM VIS [550NM EXT/KM].
!       SPSU50  2-11KM SPRING/SUMMER PROFILE W/ 50KM VIS [550NM EXT/KM].
!       SPSU23  2-11KM SPRING/SUMMER PROFILE W/ 23KM VIS [550NM EXT/KM].
!       BASTFW  10-35KM  FALL/WINTER  BACKGROUND [550NM EXT/KM].
!       VUMOFW  10-35KM  FALL/WINTER  MODERATE VOLCANIC [550NM EXT/KM].
!       HIVUFW  10-35KM  FALL/WINTER  HIGH VOLCANIC [550NM EXT/KM].
!       EXVUFW  10-35KM  FALL/WINTER  EXTREME VOLCANIC [550NM EXT/KM].
!       BASTSS  10-35KM SPRING/SUMMER BA CKGROUND [550NM EXT/KM].
!       UPNATM  >30KM NORMAL UPPER ATMOSPHERIC PROFILE [550NM EXT/KM].
!       VUTONO  >30KM VOLCANIC TO NORMAL TRANSITION [550NM EXT/KM].
      INTEGER I_HGT
      REAL HZ2K,FAWI50,FAWI23,SPSU50,SPSU23,BASTFW,VUMOFW,HIVUFW,       &
     &  EXVUFW,BASTSS,VUMOSS,HIVUSS,EXVUSS,UPNATM,VUTONO
      COMMON/PRFD/I_HGT(34),HZ2K(3,5),FAWI50( 4:11),FAWI23( 4:11),      &
     &    SPSU50( 4:11),SPSU23( 4:11),BASTFW(12:27),VUMOFW(12:27),      &
     &    HIVUFW(12:27),EXVUFW(12:27),BASTSS(12:27),VUMOSS(12:27),      &
     &    HIVUSS(12:27),EXVUSS(12:27),UPNATM(27:34),VUTONO(27:34)
      SAVE /PRFD/
      DATA I_HGT/0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,  &
     &    18,19,20,21,22,23,24,25,30,35,40,45,50,70,100,99999/
      DATA HZ2K/              .0667, .0415, .0260, .1580, .0991, .0621, &
     &   .3790, .3790, .0621, .7700, .7700, .0621,1.9400,1.9400, .0621/
      DATA FAWI50/                                                      &
     &  .011400,.006430,.004850,.003540,.002310,.001410,.000980,.000787/
      DATA FAWI23/                                                      &
     &  .027200,.012000,.004850,.003540,.002310,.001410,.000980,.000787/
      DATA SPSU50/                                                      &
     &  .014600,.010200,.009310,.007710,.006230,.003370,.001820,.001140/
      DATA SPSU23/                                                      &
     &  .034600,.018500,.009310,.007710,.006230,.003370,.001820,.001140/
      DATA BASTFW/        7.14E-4, 6.64E-4, 6.23E-4, 6.45E-4,           &
     &  6.43E-4, 6.41E-4, 6.00E-4, 5.62E-4, 4.91E-4, 4.23E-4,           &
     &  3.52E-4, 2.95E-4, 2.42E-4, 1.90E-4, 1.50E-4, 3.32E-5/
      DATA VUMOFW/        1.79E-3, 2.21E-3, 2.75E-3, 2.89E-3,           &
     &  2.92E-3, 2.73E-3, 2.46E-3, 2.10E-3, 1.71E-3, 1.35E-3,           &
     &  1.09E-3, 8.60E-4, 6.60E-4, 5.15E-4, 4.09E-4, 7.60E-5/
      DATA HIVUFW/        2.31E-3, 3.25E-3, 4.52E-3, 6.40E-3,           &
     &  7.81E-3, 9.42E-3, 1.07E-2, 1.10E-2, 8.60E-3, 5.10E-3,           &
     &  2.70E-3, 1.46E-3, 8.90E-4, 5.80E-4, 4.09E-4, 7.60E-5/
      DATA EXVUFW/        2.31E-3, 3.25E-3, 4.52E-3, 6.40E-3,           &
     &  1.01E-2, 2.35E-2, 6.10E-2, 1.00E-1, 4.00E-2, 9.15E-3,           &
     &  3.13E-3, 1.46E-3, 8.90E-4, 5.80E-4, 4.09E-4, 7.60E-5/
      DATA BASTSS/        7.99E-4, 6.41E-4, 5.17E-4, 4.42E-4,           &
     &  3.95E-4, 3.82E-4, 4.25E-4, 5.20E-4, 5.81E-4, 5.89E-4,           &
     &  5.02E-4, 4.20E-4, 3.00E-4, 1.98E-4, 1.31E-4, 3.32E-5/
      DATA VUMOSS/        2.12E-3, 2.45E-3, 2.80E-3, 2.89E-3,           &
     &  2.92E-3, 2.73E-3, 2.46E-3, 2.10E-3, 1.71E-3, 1.35E-3,           &
     &  1.09E-3, 8.60E-4, 6.60E-4, 5.15E-4, 4.09E-4, 7.60E-5/
      DATA HIVUSS/        2.12E-3, 2.45E-3, 2.80E-3, 3.60E-3,           &
     &  5.23E-3, 8.11E-3, 1.20E-2, 1.52E-2, 1.53E-2, 1.17E-2,           &
     &  7.09E-3, 4.50E-3, 2.40E-3, 1.28E-3, 7.76E-4, 7.60E-5/
      DATA EXVUSS/        2.12E-3, 2.45E-3, 2.80E-3, 3.60E-3,           &
     &  5.23E-3, 8.11E-3, 1.27E-2, 2.32E-2, 4.85E-2, 1.00E-1,           &
     &  5.50E-2, 6.10E-3, 2.40E-3, 1.28E-3, 7.76E-4, 7.60E-5/
      DATA UPNATM/        3.32E-05, 4.30E-06, 1.67E-06,                 &
     &          8.00E-07, 4.20E-07, 3.20E-08, 1.90E-10, 0./
      DATA VUTONO/        7.60E-05, 8.00E-06, 1.67E-06,                 &
     &          8.00E-07, 4.20E-07, 3.20E-08, 1.90E-10, 0./
      END BLOCK DATA PRFDTA
