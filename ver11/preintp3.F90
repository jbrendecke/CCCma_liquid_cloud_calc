!> \file preintp3.F90
!>\brief Determine the parameters needed to interpolate in pressure for gas optical properties
!!
!! @author Jiangnan Li
!
subroutine preintp3(inpt, dip, p, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI.     FOR NEW CKD WITH ONLY 20 STANDARD LEVELS
  !     * may 05,2006 - m.lazare. new version for gcm15e:
  !     *                         - implement rpn fix for inpt/inptm and
  !     *                           cosmetic reorganization.
  !     * original version preintp by jiangnan li.
  !----------------------------------------------------------------------
  use ckdlw5, only: ntl
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: dip !< Interpolation facter between two neighboring
                                                 !! standard input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: p !< Pressure at layer mid-point \f$[hPa]\f$
  integer, intent(inout), dimension(ilg,lay) :: inpt !< Index of standard input pressures for each model layer \f$[1]\f$
  !==================================================================
  !     determines the pressure interpretation points
  !     INPT for 20 pressure level up to 0.001 mb which used in GH case
  !==================================================================
  !
  integer :: i
  integer :: inp1
  integer :: inpdif
  integer :: j
  integer :: jends
  integer :: k
  integer :: m
  integer :: n
  real, parameter, dimension(ntl) :: standp = &   !< 20 standard pressure levels, the
                                                  !< gaseous absorption is pre-calculated
                                                  !< based these levels, then do interpolation
      [1.000E-02, 5.0000E-01,   1.2180,   1.8075,    2.6824, &
           3.9806,     5.9072,   8.7662,  13.0091,   19.3054, &
          28.6491,    42.5151,  63.0922,  93.6284,  138.9440, &
         206.1920,   305.9876, 454.0837, 673.8573, 1000.0000]
  !
  jends = 19
  do k = 1, lay
    do i = il1, il2
      inpt(i,k)   =  0
    end do
    !
    do j = 1, jends
      do i = il1, il2
        if (p(i,k) > standp(j)) then
          inpt(i,k) =  inpt(i,k) + 1
        end if
      end do
    end do
    !
    !----------------------------------------------------------------------c
    !     calculate arrays dip and dit required later for gasopt routines. c
    !     also, set values of inpt for a given level to be negative if all c
    !     longitude values are the same. this is also used in the gasopt   c
    !     routines to improve performance by eliminating the unnecessary   c
    !     indirect-addressing if inpl is negative for a given level.       c
    !     note that for inpt=0, it is assumed that levels are more or      c
    !     less horizontal in pressure, so scaling by -1 still preserves    c
    !     the value of zero and no indirect-addressing is done in the      c
    !     gasopt routines.                                                 c
    !----------------------------------------------------------------------c
    !
    inpdif =  0
    inp1   =  inpt(1,k)
    do i = il1, il2
      if (inpt(i,k) /= inp1)  inpdif = 1
      m  =  inpt(i,k)
      n  =  m + 1
      if (m > 0) then
        dip(i,k)  = (p(i,k) - standp(m)) / (standp(n) - standp(m))
      else
        dip(i,k)  =  p(i,k) / standp(1)
      end if
    end do
    !
    if (inpdif == 0) then
      do i = il1, il2
        inpt(i,k) =  inpt(i,k) + 1000
      end do
    end if
    !
  end do
  !
  return
end subroutine preintp3
!> \file
!> Determines nearest pressure interpretation points for each layer.
!! INPT for 28 pressure level up to 0.001 mb which used in GH case
!! INPTM for 18 pressure level up to 1 mb which used in non GH case
