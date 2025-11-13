      REAL FUNCTION AERPRF(I,VIS,IHAZE,ISEASN)
!______________________________________________________________________|
!                                                                      |
!     Copyright 2009 Spectral Sciences, Inc.                           |
!                                                                      |
!     Neither Spectral Sciences, Inc. nor its employees makes any      |
!     warranty, express or implied, regarding the uses, accuracy,      |
!     or completeness of this software.                                |
!______________________________________________________________________|

!     AERPRF COMPUTES DENSITY PROFILES FOR AEROSOLS
      IMPLICIT NONE

!     ARGUMENTS:
!       I        LAYER INDEX.
!       VIS      SURFACE VISIBILITY [KM].
!       IHAZE    HAZE MODEL NUMBER.
!       ISEASN   SEASON MODEL NUMBER.
!       IVULCN   VOLCANIC AEROSOL MODEL NUMBER.
      REAL VIS
      INTEGER I,IHAZE,ISEASN

!     COMMONS:

!     /PRFD/
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
!       BASTSS  10-35KM SPRING/SUMMER BACKGROUND [550NM EXT/KM].
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

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL PRFDTA
      AERPRF=0.
      IF(IHAZE.EQ.0)RETURN


      IF(I_HGT(I).LE.2)THEN
!         INTERPOLATE BASED ON VISIBILITIES OF 50, 23, 10, 5 AND 2 KM.
          IF(VIS.GE.23.)THEN
              AERPRF=HZ2K(I,1)+(HZ2K(I,2)-HZ2K(I,1))*(1150/VIS-23)/27
          ELSEIF(VIS.GE.10.)THEN
              AERPRF=HZ2K(I,2)+(HZ2K(I,3)-HZ2K(I,2))*(230/VIS-10)/13
          ELSEIF(VIS.GE.5.)THEN
              AERPRF=HZ2K(I,3)+(HZ2K(I,4)-HZ2K(I,3))*(10/VIS-1)
          ELSE
              AERPRF=HZ2K(I,4)+(HZ2K(I,5)-HZ2K(I,4))*(10/VIS-2)/3
          ENDIF



      ELSEIF(I_HGT(I).LE.10)THEN
          IF(ISEASN.LE.1)THEN
              IF(VIS.LE.23.)THEN
                  AERPRF=SPSU23(I)
              ELSE
                  AERPRF=SPSU50(I)
                  IF(I_HGT(I).LE.4)                                     &
     &              AERPRF=AERPRF+(SPSU23(I)-AERPRF)*(1150/VIS-23)/27
              ENDIF
          ELSE
              IF(VIS.LE.23.)THEN
                  AERPRF=FAWI23(I)
              ELSE
                  AERPRF=FAWI50(I)
                  IF(I_HGT(I).LE.4)                                     &
     &              AERPRF=AERPRF+(FAWI23(I)-AERPRF)*(1150/VIS-23)/27
              ENDIF
          ENDIF


      ELSEIF(I_HGT(I).LE.30)THEN
          AERPRF=BASTSS(I)
          IF(ISEASN.LE.1)THEN
                  AERPRF=BASTSS(I)
          ELSE
                  AERPRF=BASTFW(I)
          ENDIF


      ELSE
          AERPRF=UPNATM(I)
      ENDIF
      RETURN
      END FUNCTION AERPRF
