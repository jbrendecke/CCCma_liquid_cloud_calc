!> \file lattenu7.F90
!>\brief Attenuation of thermal radiation between the top of model and top of atmosphere
!!
!! @author Jiangnan Li
!
subroutine lattenu7(tau0, ib, ig, o3, q, co2, ch4, an2o, dp, &
                    dt, gh, il1, il2, ilg, lay)
  !
  !     * jun 22,2013 - j.cole.   set the radiative effect of the moon
  !     *                         layer to zero at end (atten=0.),
  !     *                         so that balt=balx.
  !     * may 01,2012 - m.lazare. previous version lattenu5 for gcm16:
  !     *                         - calls attenue5 instead of attenue4.
  !     * feb 09,2009 - j.li.     previous version lattenu4 for gcm15h/i:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls attenue4 instead of attenue3.
  !     * apr 18,2008 - m.lazare/ previous version lattenu3 for gcm15g:
  !     *               j.li.     - calls attenue3 instead of attenue2.
  !     * may 05,2006 - m.lazare. previous version lattenu2 for gcm15e/f:
  !     *                         - pass integer :: variables "init" and
  !     *                           "nit" instead of actual integer
  !     *                           values, to "attenue" routines.
  !     * original version lattenu by jiangnan li.
  !----------------------------------------------------------------------
  use ckdlw5, only: ntl, cl1co2, cl2h2o, cl3h2o, cl4h2o, cl4ch4, cl4n2o, cl5h2o, &
                    cl6h2o, cl7co2u, cl8h2o, cl9h2o, &
                    cl1co2gh, cl3h2ogh, cl4h2ogh, cl4n2ogh, cl4ch4gh, &
                    cl5o3gh, cl7co2gh, cl8h2ogh, cl9h2ogh
  implicit none
  !
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg)   :: tau0 !< lw optical depth in moon layer \f$[1]\f$
  real, intent(in), dimension(ilg)      :: o3 !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg, lay) :: q !< Water vapor \f$[gram/gram]\f$
  real, intent(in), dimension(ilg, lay) :: co2 !< CO2 \f$[gram/gram]\f$
  real, intent(in), dimension(ilg, lay) :: ch4 !<ch4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg, lay) :: an2o !<n2o mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg)      :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg, lay) :: dt !< Temperature in layer between model top and top of atmosphere - 250 K \f$[K]\f$
  logical, intent(in)  :: gh   !< logical variable, if is true, use large
                               !< gaseous optical depth group \f$[0]\f$
  !==================================================================
  !     calculation of the attenuation for the downward flux above the
  !     model top level. since the temperature at 0.005 mb is unknown we
  !     assume it is the same as that of model top level
  !     co2 effect in moon can not be fully catched by single layer radiative transfer 
  !     in such think region, we adjust co2 amount at band 7 for better downward flux
  !     at the model top level, making a better lw transfer from TOA to model top
  !==================================================================
  !
  integer :: isl
  !
  real, dimension(ilg) :: rmu  !< Cosine of zenith angle \f$[unitless]\f$
  real, dimension(ilg) :: mmr  !< Placeholder for gas concentration (mixing ratio) at model top
  real, pointer, dimension(:, :) :: coeff
  !
  isl = 2
  rmu = 0.0   ! Placeholder - rmu is not used when isl = 2
  mmr = 0.0
  nullify(coeff)

  if (.not. gh) then
    if (ib == 1) then
      coeff => cl1co2
      mmr = co2(:, 1)
  !
    else if (ib == 2) then
      coeff => cl2h2o(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 3) then
      coeff => cl3h2o(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 4) then
      if (ig <= 3) then
        coeff => cl4h2o(:, :, ig)
        mmr = q(:, 1)
      else if (ig == 4) then
        coeff => cl4n2o(:, :, ig)
        mmr = an2o(:, 1)
      else if (ig == 5) then
        coeff => cl4ch4(:, :, ig)
        mmr = ch4(:, 1)
      endif
  !
    else if (ib == 5) then
      coeff => cl5h2o(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 6) then
      coeff => cl6h2o(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 7) then
      coeff => cl7co2u(:, :, ig)
      mmr = 1.3 * co2(:, 1)
  !
    else if (ib == 8) then
      coeff => cl8h2o(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 9) then
      coeff => cl9h2o(:, :, ig)
      mmr = q(:, 1)
  !
    end if
  else if (gh) then

    if (ib == 1) then
      coeff => cl1co2gh(:, :, ig)
      mmr = co2(:, 1)
  !
    else if (ib == 3) then
      coeff => cl3h2ogh
      mmr = q(:, 1)
  !
    else if (ib == 4) then
      if (ig <= 2) then
        coeff => cl4h2ogh(:, :, ig)
        mmr = q(:, 1)
      else if (ig <= 3) then
        coeff => cl4n2ogh(:, :, ig)
        mmr = an2o(:, 1)
      else if (ig <= 4) then
        coeff => cl4ch4gh(:, :, ig)
        mmr = ch4(:, 1)
      endif
  !
    else if (ib == 5) then
      coeff => cl5o3gh(:, :, ig)
      mmr = o3
  !
    else if (ib == 7) then
      coeff => cl7co2gh(:, :, ig)
      mmr = 1.3 * co2(:, 1)
  !
    else if (ib == 8) then
      coeff => cl8h2ogh(:, :, ig)
      mmr = q(:, 1)
  !
    else if (ib == 9) then
      coeff => cl9h2ogh(:, :, ig)
      mmr = q(:, 1)
  !
    endif
  endif

  if (associated(coeff)) &
     call attenue6(tau0, coeff, mmr, dp, dt(:, 1), rmu, isl, ntl, il1, il2, ilg)
  !
  return
end subroutine lattenu7
!> \file
!> Calculation attenuation of the longwave flux between the model top level and
!! the top of the atmosphere which is assumed to at 0.005 hPa with an assumed
!! temperature that is the same as that of model top level.
!!\n
!! **Note that currently it is assumed there is no attenuation above the model
!! top.**
