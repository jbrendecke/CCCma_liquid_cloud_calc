!> \file sattenu7.F90
!>\brief Attenuation of solar radiation between the top of model and top of atmosphere
!!
!! @author Jiangnan Li
!
subroutine sattenu7 (attens, ib, ig, rmu, o3, co2, ch4, o2, dp, &
                     dt, gh, il1, il2, ilg)
  !
  !     * jun 22,2013 - j.cole.   set the radiative effect of the moon
  !     *                         layer to zero at end (attens=1.),
  !     *                         so that balt=balx (full transmission).
  !     * may 01,2012 - m.lazare. previous version sattenu5 for gcm16:
  !     *                         - calls attenue5 instead of attenue4.
  !     * feb 09,2009 - j.li.     previous version sattenu4 for gcm15h/i:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls attenue4 instead of attenue3.
  !     * apr 21,2008 - l.solheim/ previous version sattenu3 for gcm15g:
  !     *               m.lazare/  - cosmetic change to add threadprivate
  !     *               j.li.        for common block "trace", in support
  !     *                            of "radforce" model option.
  !     *                          - calls attenue3 instead of attenue2.
  !     *                          - update o3 and add ch4 effect.
  !     * may 05,2006 - m.lazare. previous version sattenu2 for gcm15e/f:
  !     *                         - pass integer :: variables "init" and
  !     *                           "nit" instead of actual integer
  !     *                           values, to "attenue" routines.
  !     * original version sattenu by jiangnan li.
  !----------------------------------------------------------------------
  use ckdsw5, only: ntl, cs1o3, cs2o2, cs3co2u, cs4ch4, &
                    cs1o3gh, cs1o2gh, cs2o2gh, cs3co2gh, cs4ch4gh, cs4co2gh
  implicit none
  !
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg) :: attens !< Attenuated transmission factor for solar \f$[1]\f$
  real, intent(in), dimension(ilg) :: rmu !< Cosine of zenith angle \f$[1]\f$
  real, intent(in), dimension(ilg) :: o3 !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg) :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg) :: ch4 !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg) :: o2 !< O2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg) :: dt !< Layer temperature - 250 K \f$[K]\f$
  logical, intent(in) :: gh !< If true use optically thick gas group \f$[1]\f$
  !==================================================================
  !     Alculation of solar attenuation above the model top level. for
  !     band1 only o3 and o2 are considered, the contribution of other
  !     gases is small. for band 3 and 4, co2 is considered for gh
  !==================================================================
  !
  real :: dto3
  integer :: i
  integer :: im
  integer :: isl
  real :: tau
  real, pointer, dimension(:, :) :: coeff
  !
  isl = 1
  if (.not. gh) then
    if (ib == 1) then
        do i = il1, il2
          dto3      =  dt(i) + 23.13
          tau       =  1.02 * ((cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + &
                       dto3 * cs1o3(3,ig))) * o3(i)) * dp(i)
          attens(i) =  exp( - tau / rmu(i))
        end do ! loop 120
    else if(ib == 2) then
      if (ig <= 2) then
        coeff => cs2o2(:, :, ig)
        call attenue6(attens, coeff, o2, dp, dt, rmu, isl, ntl, il1, il2, ilg)
      else
        do i = il1, il2
          attens(i) =  1.0
        end do !loop 200
      endif
    else if(ib == 3) then
      if (ig <= 2) then
        do i = il1, il2
          attens(i) =  1.0
        enddo
      else
        coeff => cs3co2u(:, :, ig)
        call attenue6 (attens, coeff, co2, dp, dt, rmu, isl, ntl, il1, il2, ilg)
      endif
    else if (ib == 4) then
        coeff => cs4ch4(:, :, ig)    
        call attenue6 (attens, coeff, ch4, dp, dt, rmu, isl, ntl, il1, il2, ilg)
    endif
  !
  else if (gh) then
  !
    if (ib == 1) then
      if (ig > 1) then
        do i = il1, il2
          dto3      =  dt(i) + 23.13
          tau       =  1.02 * ((cs1o3gh(1,ig) + dto3 * (cs1o3gh(2,ig) + &
                       dto3 * cs1o3gh(3,ig))) * o3(i) + cs1o2gh(ig-1) * o2(i)) * &
                       dp(i)
          attens(i) =  exp( - tau / rmu(i))
        end do ! loop 100
      else
        do i = il1, il2
          dto3      =  dt(i) + 23.13
          tau       =  1.02 * (cs1o3gh(1,ig) + dto3 * &
                      (cs1o3gh(2,ig) + dto3 * cs1o3gh(3,ig))) * o3(i) * dp(i)
          attens(i) =  exp( - tau / rmu(i))
        end do ! loop 110
      endif
  !
    else if (ib == 2) then
      coeff => cs2o2gh(:, :, ig)
      call attenue6 (attens, coeff, o2, dp, dt, rmu, isl, ntl, il1, il2, ilg)
  !
    else if (ib == 3) then
      coeff => cs3co2gh(:, :, ig)
      call attenue6 (attens, coeff, co2, dp, dt, rmu, isl, ntl, il1, il2, ilg)
  !
    else if (ib == 4) then
      if (ig == 6) then
        call attenue6 (attens, cs4co2gh, co2, dp, dt, rmu, isl, ntl, il1, il2, ilg)
      else 
        coeff => cs4ch4gh(:, :, ig)
        call attenue6 (attens, coeff, ch4, dp, dt, rmu, isl, ntl, il1, il2, ilg)
      endif
    endif
  !
  endif ! gh
  !
  return
end subroutine sattenu7
!> \file
!>  Calculation of solar attenuation above the model top level. for
!! band1 only o3 and o2 are considered, the contribution of other
!! gases is small. for band 3 and 4, co2 is considered for gh.
