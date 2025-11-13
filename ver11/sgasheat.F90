!> \file strandngh5.F90
!>\brief Compute the downward solar flux from model top to 1 hPa for optically thick atmospheres
!!
!! @author Jiangnan Li
!
 subroutine sgasheat (hrsh2o, hrso3, hrso2, hrsco2, hrsch4, tran, &
                      refl, dp, o2g, rmug, rmu3, sgwt, isun, gh,  &
                      ib, ig, pfull, il1g, il2g, ilg, lay, lev )
  !
  !     * jun 18,2021 - j.li      for new ckd
  !
  implicit none
  !
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1g
  integer, intent(in) :: il2g
  integer, intent(in) :: ilg
  integer, intent(in) :: lay
  integer, intent(in) :: lev
  real, dimension(ilg, lay, 3),  intent(inout) :: hrsh2o  !<h2o solar heating rate\f$[k/sec]\f$
  real, dimension(ilg, lay, 12), intent(inout) :: hrso3   !<o3 solar heating rate\f$[k/sec]\f$
  real, dimension(ilg, lay, 2),  intent(inout) :: hrso2   !<o2 solar heating rate\f$[k/sec]\f$
  real, dimension(ilg, lay),     intent(inout) :: hrsco2  !<co2 solar heating rate\f$[k/sec]\f$
  real, dimension(ilg, lay),     intent(inout) :: hrsch4  !<ch4 solar heating rate\f$[k/sec]\f$
  real, dimension(ilg, 2, lev),  intent(inout) :: tran    !<downward flux \f$[w/m^2]\f$
  real, dimension(ilg, 2, lev),  intent(in) :: refl       !<upward flux \f$[w/m^2]\f$
  real, dimension(ilg, lev),     intent(in) :: pfull      !<pressure at level \f$[hpa]\f$
  real, dimension(ilg, lay),     intent(in) :: dp         !<airmass path of a layer \f$[gram/cm^2]\f$
  real, dimension(ilg, lay),     intent(in) :: o2g        !<o2 mixing ratio\f$[gram/gram]\f$
  real, dimension(ilg),          intent(in) :: rmug       !<cosine of solar zenith angle\f$[0]\f$
  real, dimension(ilg),          intent(in) :: rmu3       !<factor of rmug\f$[0]\f$
  real, dimension(ilg),          intent(in) :: sgwt       !energy portion weight for each g in k space \f$[0]\f$
  !
  integer, intent(in) :: isun(ilg)  !<rearrange lat/log loop to only sun rise region\f$[0]\f$
  logical, intent(in) :: gh         !< If is true, use large gaseous optical depth group for calculations \f$[1]\f$
  !
  real :: abss
  real :: dfnet
  real :: tau
  real :: x
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: kp1

  real, parameter :: hrcoef = 9.9565841e-05
  !
  !----------------------------------------------------------------------  !
  !     Separate the heating rate for each gas for mam project             !
  !     since the interaction between each gases, the heating rate could   !
  !     not be exactly separated. for main k intervals, HRSO3 are          !
  !     obtained for each k of visible band                                !
  !                                                                        !
  !     O3(1), O2(1),  50000 - 43000 (     cm-1)                           !
  !     O3(2),         43000 - 37500                                       !
  !     O3(3),         37500 - 35714                                       !
  !     O3(4),         35714 - 33557                                       !
  !     O3(5),         33557 - 32895                                       !
  !     O3(6),         32895 - 32185                                       !
  !     O3(7),         32185 - 31746                                       !
  !     O3(8),         31746 - 30488                                       !
  !     O3(9),         30488 - 27473                                       !
  !     O3(10),        27473 - 25000                                       !
  !     O3(11),        25000 - 19000                                       !
  !     O3(12), O2(2), 19000 - 14500 -> H2O is set to 0 above 100 mb       !
  !     H2O(1), O2(2), 8400 - 14500  -> below it H2O AND O3 are combined   !
  !     H2O(2), CO2,   4200 - 8400                                         !
  !     H2O(3), CO2,   8400 - 2500                                         !
  !     HRSCO2 + HRLCO2(1) (in longwave) are 4.3 um CO2 band               !
  !     HRSCH4         8400 - 2500  only consider that in gh               !
  !----------------------------------------------------------------------  !
  !
  if (.not. gh) then
    if (ib == 1) then
      l = 12 - ig + 1
      do k = 1, lay
        kp1 = k + 1
        do i = il1g, il2g
          j = isun(i)
          dfnet        = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
          hrso3(j,k,l) =  hrcoef * max(dfnet, 0.) / dp(i,k)
        enddo
      enddo
  !
      if (ig == 1) then
  !
  !----------------------------------------------------------------------  !
  !     separate the o2 heating out above 100 mb                           !
  !----------------------------------------------------------------------  !
  !
        do i = il1g, il2g
          tran(i,2,1)       =  1.0
        enddo
  !
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            if (pfull(i,k) < 100.0) then
              tau           = (0.511e-04 - 0.881e-05 * rmu3(i)) * o2g(i,k) * dp(i,k)
              tran(i,2,kp1) =  tran(i,2,k) * exp( - tau / rmug(i))
              j = isun(i)
              dfnet         = (tran(i,2,k) - tran(i,2,kp1)) * sgwt(i)
              x             =  hrcoef * max(dfnet,0.) / dp(i,k)
              hrso3(j,k,12) =  hrso3(j,k,12) - x
              hrso2(j,k,2)  =  x
            endif
          enddo
        enddo
  !
      endif
  !
    else if (ib == 2) then
      if (ig <= 2) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet        = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
            hrso2(j,k,2) =  hrso2(j,k,2) + hrcoef * max(dfnet,0.) / dp(i,k)
          enddo
        enddo
      else
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
            hrsh2o(j,k,1) =  hrsh2o(j,k,1) + hrcoef * max(dfnet,0.) / dp(i,k)
          enddo
        enddo
      endif
  !
  !----------------------------------------------------------------------  !
  !    in band 3 and 4, the ch4 and o2 are ignored in non-gh part          !
  !----------------------------------------------------------------------  !
  !
    else if (ib == 3) then
      if (ig /= 5) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
            hrsh2o(j,k,2) =  hrsh2o(j,k,2) + hrcoef * max(dfnet,0.) / dp(i,k)
          enddo
        enddo
      else
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet       = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
            hrsco2(j,k) =  hrsco2(j,k) + hrcoef * max(dfnet,0.) / dp(i,k)
          enddo
        enddo
      endif
  !
    else if (ib == 4) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1g, il2g
          j = isun(i)
          dfnet         = (tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
          hrsh2o(j,k,3) =  hrsh2o(j,k,3) + hrcoef * max(dfnet,0.) / dp(i,k)
        enddo
      enddo
    endif
  !
  else if (gh) then
  !
  !----------------------------------------------------------------------  !
  !     calculations of heating rates for each separate gases              !
  !----------------------------------------------------------------------  !
  !
    if (ib == 1) then
      l = 3 - ig + 1
      do k = 1, lay
        kp1 = k + 1
        do i = il1g, il2g
          j = isun(i)
          dfnet           =  tran(i,2,k) - tran(i,2,kp1)
          hrso3(j,k,l)    =  hrcoef * sgwt(i) * dfnet / dp(i,k)
        enddo
      enddo
  !
      if (ig == 3) then
  !
  !----------------------------------------------------------------------  !
  !     at ib = 1, ig = 3, o2 and o3 heating rate are combined, followingc
  !     calculation for o3 heating rate only, then using (o3 + o2) - o2    !
  !     to get o3 heating, the results for each gas could not be aacuate   !
  !     below 1 mb since the strong overlpap of two gases                  !
  !----------------------------------------------------------------------  !
  !
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            tau           =  0.900007e-01 * o2g(i,k) * dp(i,k)
            abss          =  exp( - tau / rmug(i))
            tran(i,2,kp1) =  tran(i,2,k) * abss
          enddo
        enddo
  !
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrso2(j,k,1)  =  hrcoef * sgwt(i) * dfnet / dp(i,k)
            hrso3(j,k,1)  =  hrso3(j,k,1) - hrso2(j,k,1)
          enddo
        enddo

      endif
  !
    else if (ib == 2) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1g, il2g
          j = isun(i)
          dfnet           =  tran(i,2,k) - tran(i,2,kp1)
          hrso2(j,k,2)    =  hrso2(j,k,2) + hrcoef * sgwt(i) * dfnet / dp(i,k)
        enddo
      enddo
  !
    else if (ib == 3) then
      if (ig == 2) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrsh2o(j,k,2) =  hrsh2o(j,k,2) + hrcoef * sgwt(i) *  dfnet / dp(i,k)
          enddo
        enddo
      else
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrsco2(j,k)   =  hrsco2(j,k) + hrcoef * sgwt(i) * dfnet / dp(i,k)
          enddo
        enddo
      endif
  !
    else if (ib == 4) then
  !
      if (ig == 1 .or. ig == 6) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrsh2o(j,k,3) =  hrsh2o(j,k,3) + hrcoef * sgwt(i) * dfnet / dp(i,k)
          enddo
        enddo
      else if (ig == 4) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrsch4(j,k)   =  hrcoef * sgwt(i) * dfnet / dp(i,k)
          enddo
        enddo
      else
        do k = 1, lay
          kp1 = k + 1
          do i = il1g, il2g
            j = isun(i)
            dfnet         =  tran(i,2,k) - tran(i,2,kp1)
            hrsco2(j,k)   =  hrsco2(j,k) + hrcoef * sgwt(i) * dfnet / dp(i,k)
          enddo
        enddo
      endif
    endif
  endif

  return
 end subroutine sgasheat
!> \file
