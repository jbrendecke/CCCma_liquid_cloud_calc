!> \file tline2a.F90
!>\brief Calculate optical depth for H2O and another gas (line contribution only).
!!
!! @author Jiangnan Li
!
subroutine tline2a(taug, coef1, coef2, s1, s2, dp, dip, dt, inpt, initaug, &
                   il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI      FOR NEW CKD
  !     * feb 09,2009 - j.li.     new version for gcm15h:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     * apr 18,2008 - m.lazare. previous version tline2y for gcm15g:
  !     *                         - cosmetic change of n->lc passed
  !     *                           in and used in dimension of coef
  !     *                           array(s), so that it is not redefined
  !     *                           and passed back out to calling
  !     *                           routine, changing value in data
  !     *                           statement !
  !     *                         - cosmetic change to add threadprivate
  !     *                           for common block "trace", in support
  !     *                           of "radforce" model option.
  !     *                         - work array "s2" now local instead of
  !     *                           passed-in workspace.
  !     * may 05,2006 - m.lazare. previous version tline2x for gcm15e/f:
  !     *                         - implement rpn fix for inpt.
  !     * apr 25,2003 - j.li.     previous version tline2 for gcm15d.
  !
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  integer, intent(in)                     :: initaug !< Determinng if the previous Tau is added.
  real, intent(inout), dimension(ilg,lay) :: taug !< Gas optical thickness \f$[1]\f$
  real, intent(in), dimension(5,20) :: coef1 !< Polynomial coefficients for gas 1
                                             !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                             !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                             !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(5,20) :: coef2 !< Polynomial coefficients for gas 2
                                             !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                             !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                             !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(ilg,lay) :: s1 !< Gas 1 mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: s2 !< Gas 2 mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between adjacent standard pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1]\f$
  !
  !==================================================================
  ! Calculation of optical depth for two gases (line contribution only), the
  ! gaseous absorption coefficients in units of cm2 / gram. The absorption
  ! coefficients are calculated at the temperature T for the 18 (28 for GH)
  ! pressure levels.
  !
  ! By hydrostatic law, DP = DIFP / G = RHO * DZ,
  ! where DIFP is the layer pressure difference (in mb), G is the gravity
  ! constant, RHO is air density, and DZ is layer thickness (in cm). DP is
  ! the airmass path.
  ! S is the gaseous mixing ratio, the considered gas to air
  !                S * DP = GAS mass * DZ.
  ! we call it as the gas mass path for a model layer.
  !
  ! First use the polynomial interpolation to get results at the
  ! temperature T based on the pre-calculated results at 5 temperatures.
  ! Then, do the linear interoplation for between two pressure levels,
  ! as X1 and X2.
  ! X1, X2, Y1, Y2 are the absoption coefficients in unit cm^2/gram
  !
  !=======================================================================
  !
  integer :: i
  integer :: k
  integer :: m
  integer :: n
  real :: x1
  real :: x2
  real :: y1
  real :: y2
  !
  if (initaug == 2) then
    do k = 1, lay
      do i = il1, il2
        taug(i,k) =  0.
      end do
    end do 
  end if
  do k = 1, lay
    if (inpt(1,k) < 950) then
      do i = il1, il2
        m  =  inpt(i,k)
        n  =  m + 1
        x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + dt(i,k) * &
                    (coef1(3,n) + dt(i,k) * (coef1(4,n) + &
                     dt(i,k) * coef1(5,n))))
        !
        y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + dt(i,k) * &
                    (coef2(3,n) + dt(i,k) * (coef2(4,n) + &
                     dt(i,k) * coef2(5,n))))
        if (m > 0) then
          x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + dt(i,k) * &
                    (coef1(3,m) + dt(i,k) * (coef1(4,m) + &
                     dt(i,k) * coef1(5,m))))
          !
          y1      =  coef2(1,m) + dt(i,k) * (coef2(2,m) + dt(i,k) * &
                    (coef2(3,m) + dt(i,k) * (coef2(4,m) + &
                     dt(i,k) * coef2(5,m))))
        else
          x1      =  0.0
          y1      =  0.0
        end if
        !
        taug(i,k) =  taug(i,k) + ( (x1 + (x2 - x1) * dip(i,k)) * s1(i,k) + &
                    (y1 + (y2 - y1) * dip(i,k)) * s2(i,k) ) * dp(i,k)
      end do ! loop 1000
    else
      m  =  inpt(1,k) - 1000
      n  =  m + 1
      do i = il1, il2
        x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + dt(i,k) * &
                    (coef1(3,n) + dt(i,k) * (coef1(4,n) + &
                     dt(i,k) * coef1(5,n))))
        !
        y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + dt(i,k) * &
                    (coef2(3,n) + dt(i,k) * (coef2(4,n) + &
                     dt(i,k) * coef2(5,n))))
        if (m > 0) then
          x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + dt(i,k) * &
                    (coef1(3,m) + dt(i,k) * (coef1(4,m) + &
                     dt(i,k) * coef1(5,m))))
          !
          y1      =  coef2(1,m) + dt(i,k) * (coef2(2,m) + dt(i,k) * &
                    (coef2(3,m) + dt(i,k) * (coef2(4,m) + &
                     dt(i,k) * coef2(5,m))))
        else
          x1      =  0.0
          y1      =  0.0
        end if
        !
        taug(i,k) =  taug(i,k) + ( (x1 + (x2 - x1) * dip(i,k)) * s1(i,k) + &
                    (y1 + (y2 - y1) * dip(i,k)) * s2(i,k) ) * dp(i,k)
      end do ! loop 1500
    end if
  end do ! loop 2000
  !
  return
end subroutine tline2a
!> \file
!> Calculation of optical depth for two gases (line contribution only), the
!! gaseous absorption coefficients in units of cm2 / gram. The absorption
!! coefficients are calculated at the temperature T for the 18 (28 for GH)
!! pressure levels.
!! First use the polynomial interpolation to get results at the
!! temperature T based on the pre-calculated results at 5 temperatures.
!! Then, do the linear interoplation for between two pressure levels,
!! as X1 and X2.
