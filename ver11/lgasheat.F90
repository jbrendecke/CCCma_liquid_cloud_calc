!> \file lattenu7.F90
!>\brief Attenuation of thermal radiation between the top of model and top of atmosphere
!!
!! @author Jiangnan Li
!
 subroutine lgasheat (hrlh2o, hrlco2, hrlo3, hrlch4, hrln2o, &
                      tran, refl, dp, pgw, gh, ib, ig, &
                      il1, il2, ilg, lay, lev )

  implicit none
  !
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev
  real, intent(inout), dimension(ilg,lay,9) ::  hrlh2o  !<h2o longwave heating rate\f$[k/sec]\f$
  real, intent(out),   dimension(ilg,lay,2) :: hrlco2  !<co2 longwave heating rate\f$[k/sec]\f$
  real, intent(inout), dimension(ilg,lay) :: hrlo3     !<o3 longwave heating rate\f$[k/sec]\f$
  real, intent(out),   dimension(ilg,lay) :: hrlch4    !<ch4 longwave heating rate\f$[k/sec]\f$
  real, intent(out),   dimension(ilg,lay) :: hrln2o    !<n2o longwave heating rate\f$[k/sec]\f$
  real, intent(in),    dimension(ilg,2,lev) :: tran    !<downward flux \f$[w/m^2]\f$
  real, intent(in),    dimension(ilg,2,lev) :: refl    !<upward flux \f$[w/m^2]\f$
  real, intent(in),    dimension(ilg,lay) :: dp        !<airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in) :: pgw
  logical, intent(in) :: gh  !< logical variable, if is true, use large gaseous optical depth group \f$[0]\f$
  !
  integer :: i
  integer :: k
  integer :: kp1
  real :: dfnet
  real, parameter :: hrcoef = 9.9565841e-05
  !
  if ( .not. gh) then
    if (ib /= 4 .and. ib /= 5) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet          =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlh2o(i,k,ib) =  hrlh2o(i,k,ib) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
  !
    else if (ib == 4) then
      if (ig <= 4) then      
        do k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            dfnet          =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrlh2o(i,k,ib) =  hrlh2o(i,k,ib) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      else if (ig == 5) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            dfnet       =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrln2o(i,k) =  hrln2o(i,k) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      else if (ig == 6) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            dfnet       =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrlch4(i,k) =  hrlch4(i,k) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      endif      
    !  
    else if (ib == 5) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet      =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlo3(i,k) =  hrlo3(i,k) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
    endif
  !
  else if (gh) then
  !
    if (ib == 1) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet         =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlco2(i,k,1) =  hrlco2(i,k,1) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
    endif
  !
    if (ib == 3 .or. (ib == 5 .and. ig == 2) .or. ib >= 8) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet          =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlh2o(i,k,ib) =  hrlh2o(i,k,ib) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
    endif
  !
    if (ib == 4) then
      if (ig <= 2) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            dfnet          =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrlh2o(i,k,ib) =  hrlh2o(i,k,ib) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      else if (ig == 3) then
        do k = 1, lay
        kp1 = k + 1
          do i = il1, il2
            dfnet       =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrln2o(i,k) =  hrln2o(i,k) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      else if (ig == 4) then
        do k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            dfnet       =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
            hrlch4(i,k) =  hrlch4(i,k) + hrcoef * dfnet / dp(i,k) * pgw
          enddo
        enddo
      endif
    endif
  !
    if (ib == 5 .and. ig /= 2) then
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet      =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlo3(i,k) =  hrlo3(i,k) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
    endif
  !
    if (ib == 7) then
  !
  !----------------------------------------------------------------------
  !     15 um band co2 and h2o are strongly overlapped. it is difficult
  !     to separate the cooling of them in the lower atmosphere, we put
  !     most cooling rate in h2o(7) below 100 mb (for all main k values),
  !     and assume all cooling in the upper part are from co2 (for all
  !     minor k values).
  !----------------------------------------------------------------------
  !
      do k = 1, lay
        kp1 = k + 1
        do i = il1, il2
          dfnet         =  tran(i,2,k) - tran(i,2,kp1) - refl(i,2,k) + refl(i,2,kp1)
          hrlco2(i,k,2) =  hrlco2(i,k,2) + hrcoef * dfnet / dp(i,k) * pgw
        enddo
      enddo
    endif
  !
  endif

  return
 end subroutine lgasheat
