!> \file gasoptl8.F90
!>\brief Compute the gas optical thickness for thermal wavelengths
!!
!! @author Jiangnan Li
!
subroutine gasoptl8(taug, gw, dp, ib, ig, o3, q, co2, ch4, an2o, &
                    f11, f12, inptr, inpt, mcont, dir, dip, dt, &
                    il1, il2, ilg, lay)
  !
  !     * JUN 18/2021 - J.LI. FOR NEW CKD
  !     * feb 02/2018 - j.li. corrects the water vapour continuum in band 6
  !     *                     which was not accurate for much warmer temps.
  !     * apr 21/2010 - j.li.     new version for gcm15i:
  !     *                         - revision of many bands for greater
  !     *                           accuracy.
  !     * feb 09,2009 - j.li.     previous version gasoptl6 for gcm15h:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls tline{1,2,3}z instead of
  !     *                           tline{1,2,3}y. also calls tlinehc4,
  !     *                           tconthl4 instead of tlinehc3,tconthl3.
  !     * apr 18,2008 - m.lazare/ previous version gasoptl5 for gcm15g:
  !     *               l.solheim/- cosmetic change to add threadprivate
  !     *               j.li.       for common block "trace", in support
  !     *                         - calls tline{1,2,3}y instead of
  !     *                           tline{1,2,3}x. also calls tlinehc3,
  !     *                           tconthl3 instead of tlinehc2,tconthl2.
  !     * jun 21,2006 - j.li.     previous version gasoptl4 for gcm15f:
  !     *                         - update the minor gas absorption.
  !     *                         - calls new tconthl2.
  !     * may 05,2006 - m.lazare. previous version gasoptl3 for gcm15e:
  !     *                         - pass integer :: variables "initaug" and
  !     *                           "mit" instead of actual integer
  !     *                           values, to "tline_" routines.
  !     * dec 07/2004 - j.li.     previous version gasoptl2 for gcm15d.
  !----------------------------------------------------------------------
 use ckdlw5, only : gwl1, cl1h2o, cl1co2, cl1n2o, &
                    gwl2, cl2h2o, cl2cs, cl2cf, &
                    gwl3, cl3h2o, cl3cs, cl3cf, &
                    gwl4, cl4h2o, cl4n2o, cl4ch4, cl4f12, cl4cs, cl4cf, &
                    gwl5, cl5h2o, cl5o3, cl5co2, cl5f11, cl5f12, cl5cs, cl5cf, &
                    gwl6, cl6h2o, cl6co2, cl6f11, cl6f12, cl6cs, cl6cf, &
                    gwl7, cl7co2u, cl7h2od, cl7co2d, cl7cs, cl7cf, cl7n2o, cl7o3, &
                    gwl8, cl8h2o, cl8cs, cl8cf, &
                    gwl9, cl9h2o, cl9cs, cl9cf

  implicit none
  real, intent(inout) :: gw
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: mcont!< Highest model level for water vapor continuum calculation \f$[1]\f$
  !
  real, intent(out), dimension(ilg,lay) :: taug !< Gas optical thickness \f$[1]\f$
  !
  real, intent(in), dimension(ilg,lay) :: dp  !< Airmass path for layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: o3  !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: q   !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: ch4 !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: an2o!< N2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: f11 !< CFC114 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: f12 !< CFC12 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dir !< Variable description\f$[units]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between two neighboring standard input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt  !< Layer temperature - 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inptr !< Number of the selected standard input H2O/CO2 ratio \f$[1]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1]\f$
  !==================================================================
  !     calculation of the optical depths due to nongray gaseous
  !     absorption for the infrared, in each layer for a given band ib
  !     and cumulative probability gw.
  !     from band1 to band4, the solar and infrared interaction is
  !     considered. the total solar energy considered in the infrared
  !     region is 11.9006 w / m2
  !     tline, etc., deal with line absorption and tcontl and tconthl
  !     deal with water vapor continuum
  !==================================================================
  !
  integer :: i
  integer :: initaug
  integer :: k
  integer :: lc
  !
  real, pointer, dimension(:, :) :: coeff1, coeff2, coeff3
  real, pointer, dimension(:, :, :) :: coeff4, coeff5
  !
  !=======================================================================
  !
  taug = 0.0

  if (ib == 1) then
    !
    !----------------------------------------------------------------------
    !     band (2500 - 2200 cm-1), nongray gaseous absorption of h2o and
    !     co2, n2o.
    !----------------------------------------------------------------------
    !
    call tline3a(taug, cl1h2o, cl1co2, cl1n2o, q, co2, an2o, dp, dip, dt, inpt, &
                 il1, il2, ilg, lay)
    !
    gw =  gwl1
    !
  else if (ib == 2) then
    !
    !----------------------------------------------------------------------
    !     band (2200 - 1900 cm-1), nongray gaseous absorption of h2o + n2oc
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1 => cl2h2o(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------
    !     simply add the n2o effect
    !----------------------------------------------------------------------
    !
    if (ig == 1) then
       do k = 1, lay
         do i = il1, il2
            taug(i,k) =  taug(i,k) + 1000. * (0.04 + 2.28 * (q(i,k) / &
                         (3. * q(i,k) + 50000. * an2o(i,k)))) * an2o(i,k) * dp(i,k)
         end do
       end do
       lc =  3
       call tcontl2(taug, cl2cs, cl2cf, q, dp, dip, dt, lc, inpt, mcont, &
                    il1, il2, ilg, lay)
    end if
    !
    gw =  gwl2(ig)
    !
  else if (ib == 3) then
    !
    !----------------------------------------------------------------------
    !     band (1900 - 1400 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1 => cl3h2o(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    if (ig <= 3) then
      lc =  4
      coeff1 => cl3cs(:, :, ig)
      coeff2 => cl3cf(:, :, ig)
      call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                   il1, il2, ilg, lay)
    end if

    gw =  gwl3(ig)
    !
  else if (ib == 4) then
    !
    !----------------------------------------------------------------------
    !     band4 (1100 - 1400 cm-1), overlapping absorption of h2o, n2o,
    !     ch4 and cfc12. direct mapping method for h2o and ch4 and n2o
    !     cfc are considered as minor gases
    !----------------------------------------------------------------------
    !
    coeff1 => cl4h2o(:, :, ig)
    coeff2 => cl4ch4(:, :, ig)
    coeff3 => cl4n2o(:, :, ig)
    call tline3a(taug, coeff1, coeff2, coeff3, q, ch4, an2o, dp, dip, dt, inpt, &
                 il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------
    !     Add the cfc12 effect
    !----------------------------------------------------------------------
    !
    if (ig <= 4) then
       initaug = 1
       coeff1 => cl4f12(:, :, ig)
       call tline1a(taug, coeff1, f12, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    end if
    !
    if (ig <= 5) then
      lc =  4
      coeff1 => cl4cs(:, :, ig)
      coeff2 => cl4cf(:, :, ig)
      call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                   il1, il2, ilg, lay)
    end if
    !
    gw =  gwl4(ig)
    !
  else if (ib == 5) then
    !
    !----------------------------------------------------------------------
    !     band5 (980 - 1100 cm-1), overlapping absorption of h2o and o3
    !     direct mapping method. co2 and cfc are simply added
    !----------------------------------------------------------------------
    !
    coeff1 => cl5h2o(:, :, ig)
    coeff2 => cl5o3(:, :, ig)
    coeff3 => cl5co2(:, :, ig)
    call tline3a(taug, coeff1, coeff2, coeff3, q, o3, co2, dp, dip, dt, &
                 inpt, il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------
    !     simply add the cfc11, cfc12 effect
    !----------------------------------------------------------------------
    !
    do k = 1, lay
      do i = il1, il2
        taug(i,k) = taug(i,k) + (cl5f11(ig) * f11(i,k) + &
                    cl5f12(ig) * f12(i,k)) * dp(i,k)
      end do
    end do
    !
    lc =  4
    coeff1 => cl5cs(:, :, ig)
    coeff2 => cl5cf(:, :, ig)
    call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                 il1, il2, ilg, lay)
    !
    gw =  gwl5(ig)
    !
  else if (ib == 6) then
    !
    !----------------------------------------------------------------------
    !     band (800 - 980 cm-1), nongray gaseous absorption of h2o.
    !     + cfc11 and cfc12
    !----------------------------------------------------------------------
    !
    coeff1 => cl6h2o(:, :, ig)
    coeff2 => cl6co2(:, :, ig)
    initaug = 2
    call tline2a(taug, coeff1, coeff2, q, co2, dp, dip, dt, inpt, initaug, &
                 il1, il2, ilg, lay)
    !
    coeff1 => cl6f11(:, :, ig)
    coeff2 => cl6f12(:, :, ig)
    initaug = 1
    call tline2a(taug, coeff1, coeff2, f11, f12, dp, dip, dt, inpt, initaug, &
                 il1, il2, ilg, lay)
    !
    lc =  6
    coeff1 => cl6cs(:, :, ig)
    coeff2 => cl6cf(:, :, ig)
    call tcontl2 (taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                  il1, il2, ilg, lay)
    !
    gw =  gwl6(ig)
    !
  else if (ib == 7) then
    !
    !----------------------------------------------------------------------
    !     band6 (540 - 800 cm-1), overlapping absorption of h2o and co2
    !     exact mapping method for h2o and co2, direct mapping for n2o
    !     o3 effect is simply added
    !----------------------------------------------------------------------
    !
    coeff1 => cl7co2u(:, :, ig)
    coeff4 => cl7h2od(:, :, :, ig)
    coeff5 => cl7co2d(:, :, :, ig)
    call tlinehc5(taug, coeff1, coeff4, coeff5, q, co2, dp, dip, &
                  dir, dt, inptr, inpt, il1, il2, ilg, lay)
    !
    coeff4 => cl7cs(:, :, :, ig)
    coeff5 => cl7cf(:, :, :, ig)
    call tconthl5(taug, coeff4, coeff5, q, dp, dip, dir, dt, inptr, inpt, mcont, &
                  il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------
    !     simply add the o3 and n2o effect
    !----------------------------------------------------------------------
    !
    initaug = 1
    coeff1 => cl7n2o(:, :, ig)
     coeff1 =  coeff1 * 0.7
    coeff2 => cl7o3(:, :, ig)
    call tline2a(taug, coeff1, coeff2,  an2o, o3, dp, dip, dt, inpt, initaug, &
                 il1, il2, ilg, lay)
    !
    gw =  gwl7(ig)
    !
  else if (ib == 8) then
    !
    !----------------------------------------------------------------------
    !     band (340 - 540 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1  => cl8h2o(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    if (ig <= 4) then
      lc =  6
      coeff1 => cl8cs(:, :, ig)
      coeff2 => cl8cf(:, :, ig)
      call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                   il1, il2, ilg, lay)
    end if
    !
    gw =  gwl8(ig)
    !
  else if (ib == 9) then
    !
    !----------------------------------------------------------------------
    !     band (0 - 340 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1  => cl9h2o(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    if (ig <= 2) then
      lc =  6
      coeff1 => cl9cs(:, :, ig)
      coeff2 => cl9cf(:, :, ig)
      call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                   il1, il2, ilg, lay)
    end if
    !
    gw =  gwl9(ig)
    !
  end if
  !
  return
end subroutine gasoptl8
!> \file
!> Calculation of the optical depths due to nongray gaseous
!! absorption for the infrared, in each layer for a given band ib
!! and cumulative probability gw.
!! from band1 to band4, the solar and infrared interaction is
!! considered. the total solar energy considered in the infrared
!! region is 11.9006 w / m2
!! tline, etc., deal with line absorption and tcontl and tconthl
!! deal with water vapor continuum.
