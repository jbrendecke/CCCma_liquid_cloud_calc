!> \file tcontl2.F90
!>\brief Compute optical thickness due to water vapor continuum at thermal wavelengths
!!
!! @author Jiangnan Li
!
subroutine tcontl2(taug, coef1, coef2, s, dp, dip, dt, lc, inpt, &
                   mcont, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI      FOR NEW CKD
  !     * may 05,2006 - m.lazare. new version for gcm15e:
  !     *                         - implement rpn fix for inpt.
  !     * apr 25,2003 - j.li.     previous version tcontl for gcm15d.
  !
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lc
  integer, intent(in) :: mcont
  !
  real, intent(inout), dimension(ilg,lay) :: taug !< Gas optical thickness \f$[1]\f$
  real, intent(in), dimension(5,lc) :: coef1 !< Polynomial coefficients for h2o self continuum
                                             !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                             !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                             !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(5,lc) :: coef2 !< Polynomial coefficients for h2o foreign continuum
                                             !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                             !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                             !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(ilg,lay) :: s !< H2O mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between adjacent standard input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1]\f$
  !
  integer :: i
  integer :: j
  integer :: k
  integer :: m
  integer :: n
  integer :: nc
  real :: x1
  real :: x2
  real :: y1
  real :: y2
  !==================================================================
  !     infrared water vapor continuum, coef1 is the coefficient for
  !     self, coef2 is the coefficient for foreign. the continuum only
  !     applies to the layers below 138.9440 mb or even lower region
  !     depending on each band. lc is number of level for standard
  !     pressure considered in calculating the continuum.
  !     1.608 = 28.97 / 18.016, a fctr for water vapor partial pressure
  !
  !=======================================================================
  nc =  21 - lc
  !
  do k = mcont, lay
    if (inpt(1,k) < 950) then
      do i = il1, il2
        j =  inpt(i,k)
        if (j >= nc) then
          m  =  j - nc + 1
          n  =  m + 1
          x1        =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                      dt(i,k) * (coef1(3,m) + dt(i,k) * &
                      (coef1(4,m) + dt(i,k) * coef1(5,m))))
          !
          x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                      dt(i,k) * (coef1(3,n) + dt(i,k) * &
                      (coef1(4,n) + dt(i,k) * coef1(5,n))))
          !
          y1        =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                      dt(i,k) * (coef2(3,m) + dt(i,k) * &
                      (coef2(4,m) + dt(i,k) * coef2(5,m))))
          !
          y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                      dt(i,k) * (coef2(3,n) + dt(i,k) * &
                      (coef2(4,n) + dt(i,k) * coef2(5,n))))
          !
          taug(i,k) =  taug(i,k) + &
                      ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                      dip(i,k)) * 1.608 * s(i,k) + &
                      y1 + (y2 - y1) * dip(i,k) ) * &
                      s(i,k) * dp(i,k)
        end if
      end do ! loop 100
    else
      j  =  inpt(1,k) - 1000
      m  =  j - nc + 1
      n  =  m + 1
      do i = il1, il2
        if (j >= nc) then
          x1        =  coef1(1,m) + dt(i,k) * (coef1(2,m) + &
                      dt(i,k) * (coef1(3,m) + dt(i,k) * &
                      (coef1(4,m) + dt(i,k) * coef1(5,m))))
          !
          x2        =  coef1(1,n) + dt(i,k) * (coef1(2,n) + &
                      dt(i,k) * (coef1(3,n) + dt(i,k) * &
                      (coef1(4,n) + dt(i,k) * coef1(5,n))))
          !
          y1        =  coef2(1,m) + dt(i,k) * (coef2(2,m) + &
                      dt(i,k) * (coef2(3,m) + dt(i,k) * &
                      (coef2(4,m) + dt(i,k) * coef2(5,m))))
          !
          y2        =  coef2(1,n) + dt(i,k) * (coef2(2,n) + &
                      dt(i,k) * (coef2(3,n) + dt(i,k) * &
                      (coef2(4,n) + dt(i,k) * coef2(5,n))))
          !
          taug(i,k) =  taug(i,k) + &
                      ( (x1 - y1 + (x2 - x1 - y2 + y1) * &
                      dip(i,k)) * 1.608 * s(i,k) + &
                      y1 + (y2 - y1) * dip(i,k) ) * &
                      s(i,k) * dp(i,k)
        end if
      end do ! loop 150
    end if
  end do ! loop 200
  !
  return
end subroutine tcontl2
!> \file
!>  Water vapor continuum for 540-800 cm\f$^{-1}\f$. Different from tcontl,
!! variation of mass mixing ratio of H2O/CO2 is considered.
!! coef1 is the coefficient for self, coef2 is the coefficient for
!! foreign, in COEF1(5,5,4) and COEF2(5,5,4), the first index is for
!! the polynomial fit, the second is for the different H2O/CO2 rate,
!! and the third is for vertical level.
