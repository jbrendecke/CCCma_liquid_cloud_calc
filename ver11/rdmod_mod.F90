 MODULE RDMOD

 !     * Dec 10/2017  K.vonSalzen   Revised for new aerosol opt. prop.
 !    * Feb 10/2015  K.vonSalzen   New version for gcm18:
 !    *                          - Added sizes {NHS,NRS,NVS} and
 !    *                            revised calculation for IDIM2,
 !    *                            in conjunction with changes
 !    *                            to MIX3AERO and SSALTAEROP in
 !    *                            lssub.
 !    * June 6/2013. J.Li, K.Vonsalzen. Original version for gcm17.
      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)

      INTEGER, PARAMETER  :: NBS = 4, NBL = 9, NH = 7, NR = 11, &
                             NV = 5, NF1 = 11, NF2 = 13,        &
                             NHS = 11, NRS = 6, NVS = 3,        &
                             IDIM2=NHS*NRS*NVS

      REAL, DIMENSION(IDIM2,5):: SGT2, SOMGT2
      REAL, DIMENSION(IDIM2,6) :: SEXTT2
      REAL, DIMENSION(IDIM2,9) :: SABST2

END MODULE RDMOD
