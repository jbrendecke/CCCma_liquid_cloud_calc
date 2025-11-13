!> \file tline1a.F90
!>\brief Calculate optical depth for one gas (line contribution only).
!!
!! @author Jiangnan Li
!
subroutine tline1a(taug, coef1, s, dp, dip, dt, inpt, &
                   iplus, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI.     FOR NEW CKD
  !     * feb 09,2009 - j.li.     new version for gcm15h:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     * apr 18,2008 - m.lazare. previous version tline1y for gcm15g:
  !     *                         - cosmetic change of n->lc passed
  !     *                           in and used in dimension of coef
  !     *                           array(s), so that it is not redefined
  !     *                           and passed back out to calling
  !     *                           routine, changing value in data
  !     *                           statement !
  !     *                         - cosmetic change to add threadprivate
  !     *                           for common block "trace", in support
  !     *                           of "radforce" model option.
  !     *                         - work array "s1" now local instead of
  !     *                           passed-in workspace.
  !     * may 05,2006 - m.lazare. previous version tline1x for gcm15e/f:
  !     *                         - cosmetic cleanup to eliminate duplicate
  !     *                           coding between iplus=1 and iplus=2, by
  !     *                           initializing taug if iplus=2 (this was
  !     *                           the only difference between the two
  !     *                           sections).
  !     *                         - implement rpn fix for inpt.
  !     * apr 25,2003 - j.li.     previous version tline1 for gcm15d.
  !
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: iplus
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: taug !< Gaseous optical depth \f$[1]\f$
  real, intent(in), dimension(5,20) :: coef1 !< Polynomial coefficients for temperature \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$
  real, intent(in), dimension(ilg,lay) :: s !< Gaseous mixing ratio  \f$[g/g]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[g/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation factor for pressure \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature minus 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Index of the selected standard input pressure \f$[1]\f$
  !
  !==================================================================
  ! Calculation of optical depth for single gas (line contribution only), the
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
  ! X1 and X2 are the absoption coefficients in unit cm^2/gram
  !
  !=======================================================================
  !
  integer :: i
  integer :: k
  integer :: m
  integer :: n
  real :: x1
  real :: x2
  !
  !     * initialize taug if iplus=2.
  !
  if (iplus == 2) then
    do k = 1, lay
      do i = il1, il2
        taug(i,k)     =  0.
      end do
    end do ! loop 50
  end if
  !
  do k = 1, lay
    if (inpt(1,k) < 950) then
      do i = il1, il2
        m  =  inpt(i,k)
        n  =  m + 1
        x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                    dt(i,k) * (coef1(3,n) + dt(i,k) * &
                    (coef1(4,n) + dt(i,k) * coef1(5,n))))
        if (m > 0) then
          x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                    dt(i,k) * (coef1(3,m) + dt(i,k) * &
                    (coef1(4,m) + dt(i,k) * coef1(5,m))))
        else
          x1      =  0.0
        end if
        !
        taug(i,k) =  taug(i,k) + (x1 + (x2 - x1) * dip(i,k)) * &
                    s(i,k) * dp(i,k)
      end do ! loop 1000
    else
      m  =  inpt(1,k) - 1000
      n  =  m + 1
      do i = il1, il2
        x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                    dt(i,k) * (coef1(3,n) + dt(i,k) * &
                    (coef1(4,n) + dt(i,k) * coef1(5,n))))
        if (m > 0) then
          x1      =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                    dt(i,k) * (coef1(3,m) + dt(i,k) * &
                    (coef1(4,m) + dt(i,k) * coef1(5,m))))
        else
          x1      =  0.0
        end if
        !
        taug(i,k) =  taug(i,k) + (x1 + (x2 - x1) * dip(i,k)) * &
                    s(i,k) * dp(i,k)
      end do ! loop 1500
    end if
  end do ! loop 2000
  !
  return
end subroutine tline1a
!> \file
!> Calculate the optical depth for one gas due to line contribution only for a particular
!! temperature and pressure.
!!\n
!! The absorption coefficients are computed by first using polynomial interpolation
!! to get results at the temperature \f$ T \f$ based on the pre-calculated results at
!! 5 temperatures.  We then linearly interpolate the resulting absorption coefficients
!! X1 and X2 between two pressure levels, multiply by the gas mass path, and return the
!! gas optical depth.
!!\n
!!\n
!! The gas mass path is defined as follows:
!!\n
!! Note that the layer thickness \f$ DP \f$ is given by hydrostatic law,
!!\n
!! \f{equation}{
!!  DP = DIFP / G = RHO * DZ,
!! \f}
!!\n
!! where \f$ DP \f$ is the airmass path, \f$ DIFP \f$ is the layer pressure difference (in hPa),
!! \f$ G \f$ is the gravity constant in m/s^2, \f$ RHO \f$ is air density in kg/m^3, and
!! \f$ DZ \f$ is layer thickness (in cm). If \f$ S \f$ is the gaseous mixing ratio,
!! the considered gas to air is then
!!\n
!! \f{equation}{
!!   S * DP = M_{gas} * DZ
!! \f}
!!\n
!! which we call the gas mass path for a model layer.
!!\n
!! The absorption coefficients are computed by first using polynomial interpolation
!! to get results at the temperature \f$ T \f$ based on the pre-calculated results at
!! 5 temperatures.  We then linearly interpolate between two pressure levels,
!! as X1 and X2 (where X1 and X2 are absorption coefficients in unit cm^2/g).
