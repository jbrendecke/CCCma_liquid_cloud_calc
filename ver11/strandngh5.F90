!> \file strandngh5.F90
!>\brief Compute the downward solar flux from model top to 1 hPa for optically thick atmospheres
!!
!! @author Jiangnan Li
!
subroutine strandngh5 (tran, gwgh, atten, taua, tauoma, taucs, tauomc, &
                       cld, rmu, dp, o3, q, co2, ch4, an2o, o2, ib, ig,      &
                       inpt, dip, dt, cut, il1, il2, ilg, lay, lev)
  !
  !     * JUN 18,2021 - J.LI      FOR NEW CKD
  !     * jun 02,2015 - m.lazare/ new version for gcm19:
  !     *               j.cole:   - add tiled radiation calculations
  !     *                           (ie "trant")
  !     * feb 09,2009 - j.li.     previous version strandngh4 for gcm15h
  !     *                         through gcm18:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls tline{1,2}z instead of tline{1,2}y.
  !     * apr 18,2008 - m.lazare/ previous version strandngh3 for gcm15g:
  !     *               l.solheim/- cosmetic change to add threadprivate
  !     *               j.li.       for common block "trace", in support
  !     *                           of "radforce" model option.
  !     *                         - calls tline{1,2}y instead of tline{1,2}x.
  !     *                         - updating o3, adding ch4 and using
  !     *                           kurucz solar function.
  !     * may 05,2006 - m.lazare. previous version strandngh2 for gcm15e/f:
  !     *                         - pass integer :: variables "initaug" and
  !     *                           "nit" instead of actual integer
  !     *                           values, to "tline_" routines.
  !     * apr 25,2003 - j.li.     previous version strandngh for gcm15d.
  !----------------------------------------------------------------------
  use ckdsw5, only: gws1gh, cs1o3gh, cs1o2gh,  &
                    gws2gh, cs2o2gh,  &
                    gws3gh, cs3h2ogh, cs3co2gh, cs3ch4gh, &
                    gws4gh, cs4h2ogh, cs4ch4gh, cs4co2gh
  implicit none
  real, intent(in)    :: cut
  real, intent(inout) :: gwgh
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  !
  real, intent(out), dimension(ilg,2,lev) :: tran !< Downward solar flux \f$[W/m^2]\f$
  real, intent(in), dimension(ilg) :: atten !< Attenuated downward solar flux from
                                            !! top of atmosphere to model top \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,lay) :: taua !< Gas plus aerosol optical thickness \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauoma !< Gaseous plus aerosol optical depth times aerosol
                                                 !! single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taucs !< Cloud optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauomc !< Cloud optical depth times cloud single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: cld !< Cloud fraction \f$[1]\f$
  real, intent(in), dimension(ilg)     :: rmu !< 1/(cosine of solar zenith angle) \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: o3 !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: q !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: ch4 !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: an2o !< N2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: o2 !< O2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between two neighboring
                                              !! standard input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Index of standard input pressures for each model layer \f$[1]\f$
  !
  !==================================================================
  !     Calculation of the downward solar flux under the condition that the
  !     extinction coefficient of gas is very large (GH condition), the
  !     scattering effects can be neglected, thus no upward flux. Since
  !     the cloud optical depth is much smaller than the gaseous optical depth,
  !     the cloud effect is very small and be treated simply by beer's law
  !==================================================================
  !
  real :: absc
  real :: dto3
  real :: dtr1
  integer :: i
  integer :: im
  integer :: k
  integer :: kp1
  integer :: initaug
  !
  !     * internal work array.
  !
  real, dimension(ilg,lay) :: taug
  real, pointer, dimension(:, :) :: coeff1, coeff2, coeff3
  real, parameter, dimension(3) :: alpha = [0.980393, 1.043276, 1.410271]
  !
  !=======================================================================
  !
  do i = il1, il2
    tran(i,1,1)           =  atten(i)
    tran(i,2,1)           =  atten(i)
  end do
  !
  if (ib == 1) then
    !
    !----------------------------------------------------------------------
    !     band1 for uvc (35700 - 50000 cm-1), nongray gaseous absorption
    !     of o2  and o3.
    !----------------------------------------------------------------------
    !
    if (ig > 1) then
      do k = 1, lay
        do i = il1, il2
          dto3            =  dt(i,k) - 23.13
          taug(i,k)       =  1.02 * ((cs1o3gh(1,ig) +  dto3 * (cs1o3gh(2,ig) + &
                             dto3 * cs1o3gh(3,ig))) * o3(i,k) * alpha(ig) + &
                             cs1o2gh(ig-1) * o2(i,k)) * dp(i,k) + taua(i,k)
        end do
      end do
    else
      do k = 1, lay
        do i = il1, il2
          dto3            =  dt(i,k) - 23.13
          taug(i,k)       =  1.02 * (cs1o3gh(1,ig) + dto3 * (cs1o3gh(2,ig) + &
                             dto3 * cs1o3gh(3,ig))) * o3(i,k) * alpha(ig) * &
                             dp(i,k) + taua(i,k)
        end do
      end do
    end if
    gwgh =  gws1gh(ig)
    !
  else if (ib == 2) then
    !
    !----------------------------------------------------------------------
    !     band (8400 - 14500 cm-1), nongray gaseous absorption of o2
    !     and o3.
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1 => cs2o2gh (:, :, ig)
    call tline1a (taug, coeff1, o2, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay) 
    !
    gwgh =  gws2gh(ig)
    !
  else if (ib == 3) then
    !
    !----------------------------------------------------------------------
    !     band (4200 - 8400 cm-1), nongray gaseous absorption of h2o and
    !     co2.
    !----------------------------------------------------------------------
    !
    initaug = 2
    coeff1 => cs3h2ogh(:, :, ig)
    coeff2 => cs3co2gh(:, :, ig)
    call tline2a (taug, coeff1, coeff2, q, co2, dp, dip, dt, inpt, initaug, &
                  il1, il2, ilg, lay)
    !
    if (ig <=3) then
      initaug = 1
      coeff1 => cs3ch4gh(:, :, ig)
      call tline1a (taug, coeff1, ch4, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    endif
    gwgh =  gws3gh(ig)
    !
  else if (ib == 4) then
    !
    !----------------------------------------------------------------------
    !     band (2500 - 4200 cm-1), nongray gaseous absorption of h2o, co2,
    !     and ch4.
    !----------------------------------------------------------------------
    !
    coeff1 => cs4h2ogh(:, :, ig)
    coeff2 => cs4ch4gh(:, :, ig)      
    coeff3 => cs4co2gh(:, :, ig)
    call tline3a(taug, coeff1, coeff2, coeff3, q, ch4, co2, dp, dip, dt, &
                 inpt, il1, il2, ilg, lay)
    !
    gwgh =  gws4gh(ig)
    !
  end if

  do k = 1, lay
    kp1 = k + 1
    do i = il1, il2
      dtr1            =  exp(- (taug(i, k) + taua(i,k) - tauoma(i,k)) / rmu(i))
      tran(i,1,kp1)   =  tran(i,1,k) * dtr1
      !
      if (cld(i,k) < cut) then
        tran(i,2,kp1) =  tran(i,2,k) * dtr1
      else
        absc          = (1.0 - cld(i,k)) * dtr1 + cld(i,k) * &
                         exp(- (taug(i,k) + taua(i,k) + taucs(i,k) - tauomc(i,k)) / rmu(i))
        tran(i,2,kp1) =  tran(i,2,k) * absc
      end if
    end do
  end do
  !
  return
end subroutine strandngh5
!> \file
!>  Calculation of the downward solar flux under the condition that the
!! extinction coefficient of gas is very large (GH condition), the
!! scattering effects can be neglected, thus no upward flux. Since
!! the cloud optical depth is much smaller than the gaseous optical depth,
!! the cloud effect is very small and be treated simply by beer's law
