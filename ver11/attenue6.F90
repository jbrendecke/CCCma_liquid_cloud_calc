!> \file attenue6.F90
!>\brief Attenuation of radiation between the top of model and top of atmosphere.
!!
!! @author Jiangnan Li
!
subroutine attenue6(atten, coef1, s1, dp, dt, rmu, isl, lc, il1, il2, ilg)
  !
  !     * may 01,2012 - m.lazare. new version for gcm16:
  !     *                         - corrects bug which (it turns out)
  !     *                           doesn't change the bit-pattern
  !     *                           because of pfull(1)~0.5 being
  !     *                           between pressure levels in
  !     *                           standp of routine preintp2, this
  !     *                           part of the routine is never
  !     *                           activated at the moment.
  !     *                           otherwise, there would have been a
  !     *                           memory fault. this refers to the
  !     *                           calculation of "n" at the start of
  !     *                           loop 2000.
  !     * feb 09,2009 - j.li.     previous version attenue4 for gcm15h/i:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     * apr 18,2008 - m.lazare. previous version attenue3 for gcm15g:
  !     *                         - cosmetic change of n->lc passed
  !     *                           in and used in dimension of coef
  !     *                           array(s), so that it is not redefined
  !     *                           and passed back out to calling
  !     *                           routine, changing value in data
  !     *                           statement !
  !     *                         - cosmetic change to add threadprivate
  !     *                           for common block "trace", in support
  !     *                           of "radforce" model option.
  !     * dec 05,2007 - m.lazare/ new version for gcm15g:
  !     *               j.li.     - work array "s1" now local instead
  !     *                           of passed in.
  !     *                         - methane contribution added through
  !     *                           ng=4.
  !     * may 05,2006 - m.lazare. previous version attenue2 for gcm15e/f:
  !     *                         - implement rpn fix for inpt.
  !     *                         previous version attenue by jiangnan li.
  !----------------------------------------------------------------------
  !     This subroutine calculates the downward flux attenuation above
  !     the model top level
  !----------------------------------------------------------------------
  !
  implicit none
  !
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: isl
  integer, intent(in) :: lc
  real, intent(inout), dimension(ilg) :: atten !< Attenuated transmission factor for solar and thermal radiation \f$[1]\f$
  real, intent(in), dimension(5, lc)  :: coef1 !< Polynomial coefficients used in equation for extinction
                                             !! \f$[cm^2/gram,cm^2/gram K^(-1),cm^2/gram
                                             !! K^(-2),cm^2/gram K^(-3),cm^2/gram K^(-4)]\f$
  real, intent(in), dimension(ilg) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg) :: s1 !< Gas mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg) :: dt !< Temperature in moon layer - 250 K \f$[K]\f$
  real, intent(in), dimension(ilg) :: rmu !< Cosine of zenith angle (\f$\mu\f$) \f$[1]\f$
  integer :: i
  real :: tau
  !==================================================================
  !     This subroutine calculates the transmission factor for downward
  !     flux iattenuation above the model top level for both solar and
  !     longwave
  !
  !     ISL = 1 for solar, ISL = 2 for infrared.
  !     NG = 1, H2O; NG = 2, O3; NG = 3, CO2; NG = 6, O2
  !==================================================================
  !
  !
  if (isl == 1) then
     do i = il1, il2
        tau= (coef1(1,1) + dt(i) * (coef1(2,1) + dt(i) * &
             (coef1(3,1) + dt(i) * (coef1(4,1) + dt(i) * &
              coef1(5,1)))) ) * s1(i) * dp(i)
        !
        atten(i) = exp( - tau / rmu(i))
     end do 
  !
  else if (isl == 2) then
     do i = il1, il2
        tau= (coef1(1,1) + dt(i) * (coef1(2,1) + dt(i) * &
             (coef1(3,1) + dt(i) * (coef1(4,1) + dt(i) * &
              coef1(5,1)))) ) * s1(i) * dp(i)
        !
        atten(i) =  tau 
     end do 
  end if
  !
  return
end subroutine attenue6
!> \file
!> Calculates the transmission factor for downward flux attenuation between the
!! the model top level and the top of the atmosphere for both solar and
!! thermal radiation.  For solar radiation Beer's law is used and the attenuation
!! is \f$exp(\frac{\tau}{\mu})\f$ and for thermal radiation the radiation is
!! diffuse and the attenuation is computed using \f$exp(\frac{\tau}{\mu_{diffuse}})\f$
!! where \f$\mu_{diffuse}f$ is the diffuse factor.
!!\n
!! It is assumed that the top of atmosphere has a pressure of 0.0005 hPa and a
!! temperature of 210 K.
