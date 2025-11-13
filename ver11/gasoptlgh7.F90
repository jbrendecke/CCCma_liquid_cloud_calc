!> \file gasoptlgh7.f90
!>\brief Compute large gas optical thicknesses at thermal wavelengths
!!
!! @author Jiangnan Li
!
subroutine gasoptlgh7(taug, gwgh, dp, ib, ig, o3, q, co2, ch4, &
                      an2o, inpt, mcont, dip, dt, il1, il2, ilg, lay)
  !
  !     * JUN 18/2021 - J.LI      FOR NEW CKD
  !     * apr 22/2010 - j.li.     new version for gcm15i:
  !     *                         - add one extra term to bandl4gh for
  !     *                           greater accuracy.
  !     * feb 09,2009 - j.li.     previous version gasoptlgh5 for gcm15h:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls tline{1,2,3}z instead of
  !     *                           tline{1,2,3}y.
  !     * apr 18,2008 - m.lazare/ previous version gasoptlgh4 for gcm15g:
  !     *               j.li.     - cosmetic change to use scalar variable
  !     *                           "initaug" (=2) in calls to tline1y
  !     *                           instead of the actual number itself.
  !     *                           similar cosmetic change to use ntl(=28)
  !     *                           scalar variable instead of the actual
  !     *                           number itself in calls to all "tline_"
  !     *                           routines.
  !     *                         - calls tline{1,2,3}y instead of
  !     *                           tline{1,2,3}x.
  !     * may 05/2006 - m.lazare. previous version gasoptlgh3 for gcm15e/f:
  !     *                         - calls new versions of:
  !     *                           tline1x,tline2x,tline3y,tcontl1,
  !     *                           tconthl1.
  !     * dec 07/2004 - j.li. previous version gasoptlgh2 for gcm15c/gcm15d:
  !     *                     bugfix to reverse cl4n2ogh and cl4ch4gh
  !     *                     in bandl4gh common block to be consistent
  !     *                     with data defined in ckdlw2.
  !     * apr 25/2003 - j.li. previous version gasoptlgh for gcm15b.
  !----------------------------------------------------------------------
  use ckdlw5, only: gwl1gh, cl1co2gh, cl1n2ogh, &
                    gwl3gh, cl3h2ogh, &
                    gwl4gh, cl4h2ogh, cl4n2ogh, cl4ch4gh, &
                    gwl5gh, cl5h2ogh, cl5o3gh, cl5co2gh, cl5csgh, cl5cfgh, &
                    gwl7gh, cl7h2ogh, cl7co2gh, &
                    gwl8gh, cl8h2ogh, &
                    gwl9gh, cl9h2ogh
  implicit none
  real, intent(inout) :: gwgh
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: mcont
  !
  real, intent(out), dimension(ilg,lay) :: taug !< Gas optical depth \f$[1]\f$
  !
  real, intent(in), dimension(ilg,lay) :: dp  !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: q   !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: o3  !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: ch4 !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: an2o!< N2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between two neighboring standard
                                              !  input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt  !< layer temperature - 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< level number of the standard input pressures \f$[0]\f$
  !==================================================================
  !     the same as gasoptl but for the k intervals close to 1 in the
  !     accumulated probability space, where the gaseous absorption
  !     coeffcient is several order larger than the k intervals not
  !     1 in the accumulated probability space. For such large absorption,
  !     the radiative transfer can be much simplified.
  !     using tlinE, etc. to deal with line absorption and tcontl to
  !     deal with water vapor continumm
  !==================================================================
  !
  integer :: lc
  !
  real, pointer, dimension(:, :) :: coeff1, coeff2, coeff3
  !
  !     * initaug is a switch used in tline1y (as "iplus") which
  !     * initializes taug to zero if its value is two. this is what
  !     * we require throughout this routine.
  !
  integer, parameter :: initaug = 2
  !=======================================================================
  taug = 0.0

  if (ib == 1) then
    !
    !----------------------------------------------------------------------
    !     band (2500 - 2200 cm-1), nongray gaseous absorption of co2.
    !----------------------------------------------------------------------
    !
    coeff1 => cl1co2gh(:, :, ig)
    coeff2 => cl1n2ogh(:, :, ig)

    call tline2a(taug, coeff1, coeff2, co2, an2o, dp, dip, dt, inpt, initaug, &
                 il1, il2, ilg, lay)
    !
    gwgh =  gwl1gh(ig)
    !
    !----------------------------------------------------------------------
    !     band (2200 - 1900 cm-1), no gh
    !----------------------------------------------------------------------
  else if (ib == 3) then
    !
    !----------------------------------------------------------------------
    !     band (1900 - 1400 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    call tline1a(taug, cl3h2ogh, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    gwgh =  gwl3gh
    !
  else if (ib == 4) then
    !
    !----------------------------------------------------------------------
    !     band3 (1100 - 1400 cm-1), overlapping absorption of h2o, n2o,
    !     and ch4. direct mapping method for h2o and ch4 and n2o
    !----------------------------------------------------------------------
    !
    coeff1 => cl4h2ogh(:, :, ig)
    coeff2 => cl4ch4gh(:, :, ig)
    coeff3 => cl4n2ogh(:, :, ig)
    call tline3a(taug, coeff1, coeff2, coeff3, q, ch4, an2o, dp, dip, dt, inpt, &
                 il1, il2, ilg, lay)
    !
    gwgh =  gwl4gh(ig)
    !
  else if (ib == 5) then
    !
    !----------------------------------------------------------------------
    !     band5 (980 - 1100 cm-1), overlapping absorption of h2o and o3
    !     direct mapping method
    !----------------------------------------------------------------------
    !
    coeff1 => cl5h2ogh(:, :, ig)
    coeff2 => cl5o3gh(:, :, ig)
    coeff3 => cl5co2gh(:, :, ig)
    call tline3a(taug, coeff1, coeff2, coeff3, q, o3, co2, dp, dip, dt, inpt, &
                 il1, il2, ilg, lay)
    !
    lc =  4
    coeff1 => cl5csgh(:, :, ig)
    coeff2 => cl5cfgh(:, :, ig)
    call tcontl2(taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                 il1, il2, ilg, lay)
    !
    gwgh =  gwl5gh(ig)
    !
    !----------------------------------------------------------------------
    !     band (800 - 980 cm-1), no gh
    !----------------------------------------------------------------------
    !
  else if (ib == 7) then
    !
    !----------------------------------------------------------------------
    !     band6 (540 - 800 cm-1), overlapping absorption of h2o and co2
    !     direct mapping method. for ig > 4, the contribution by h2o is
    !     very small.
    !----------------------------------------------------------------------
    !
    if (ig <= 3) then
      coeff1 => cl7h2ogh(:, :, ig)
      coeff2 => cl7co2gh(:, :, ig)
      call tline2a(taug, coeff1, coeff2, q, co2, dp, dip, dt, inpt, initaug, &
                   il1, il2, ilg, lay)
    !
    else
      coeff1 => cl7co2gh(:, :, ig)
      call tline1a(taug, coeff1, co2, dp, dip, dt, inpt, initaug, &
                   il1, il2, ilg, lay)
    end if
    !
    gwgh =  gwl7gh(ig)
    !
  else if (ib == 8) then
    !
    !----------------------------------------------------------------------
    !     band (340 - 540 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    coeff1  => cl8h2ogh(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    gwgh =  gwl8gh(ig)
    !
  else if (ib == 9) then
    !
    !----------------------------------------------------------------------
    !     band (0 - 340 cm-1), nongray gaseous absorption of h2o.
    !----------------------------------------------------------------------
    !
    coeff1  => cl9h2ogh(:, :, ig)
    call tline1a(taug, coeff1, q, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    gwgh =  gwl9gh(ig)
    !
  end if
  !
  return
end subroutine gasoptlgh7
!> \file
!> Compute gas optical thickness for correlated \f$k\f$-distribution quadrature points
!! close to 1 in the accumulated probability space, where the gaseous absorption
!! coeffcient is several order larger than the \f$k\f$ intervals not
!! 1 in the accumulated probability space. For such large absorption,
!! the radiative transfer can be much simplified.
!! using tlinE, etc. to deal with line absorption and tcontl to
!! deal with water vapor continuum
