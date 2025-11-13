!> \file tlinehc5
!>\brief Compute optical thickness due to H2O and CO2 in 540-800 cm\f$^{-1}\f$
!!
!! @author Jiangnan Li
!
subroutine tlinehc5(taug, coef2u, coef1d, coef2d, q, co2, &
                    dp, dip, dir, dt, inptr, inpt, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI.     FOR NEW CKD
  !     * feb 09,2009 - j.li.     new version for gcm15h:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     * apr 18,2008 - m.lazare. previous version tlinehc3 for gcm15g:
  !     *                         - cosmetic change to add threadprivate
  !     *                           for common block "trace", in support
  !     *                           of "radforce" model option.
  !     * may 05,2006 - m.lazare. previoius version tlinehc2 for gcm15e:
  !     *                         - implement rpn fix for inpt.
  !     * apr 25,2003 - j.li.     previous version tlinehc for gcm15d.
  !----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: taug !< Gas optical thickness \f$[1]\f$
  real, intent(in), dimension(5,8,7) :: coef1d !< Polynomial coefficients for H2O by direct mapping
                                               !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                               !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                               !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(5,13) :: coef2u !< Polynomial coefficients for CO2 by alternate sorting
                                              !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                              !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                              !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(5,8,7) :: coef2d !< Polynomial coefficients for CO2 by direct mapping
                                               !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                               !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                               !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(ilg,lay) :: q !< H2O mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: co2 !< CO2 mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between adjancent standard pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dir !< Interpretation factor for H2O/CO2 ratio \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f$
  !
  integer, intent(in), dimension(ilg,lay) :: inptr !< Index of the selected standard input H2O/CO2 ratio \f$[1]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1]\f$
  !
  !==================================================================
  ! Calculation of optical depth for h2o and co2 (in 540-800 cm, line
  ! contribution only), since this two gases are overlapped very strong,
  ! we can not use the alternate sorting method, but use the direct
  ! mapping method by sorting the two gases together. The pre-calculations
  ! are done for 5 h2o/co2 ratios of .06,  .24,  1., 4., 16, then using
  ! the interpolation to obtain results of other ratios. This method is
  ! much more expensive than alternate sorting, we only applied it to the
  ! lower atmosphere with INPT > 11, above this level we still use alternate sorting
  ! The absorption
  ! coefficients are calculated at the temperature T for the 18 (not for GH)
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
  ! X1 and X2 for h2o, Y1 and Y2 for co2, are the absoption coefficients
  ! in unit cm^2/gram
  ! The more details are in Li and Barker J. Atmos. Sci. (2005) 62 286-309
  !
  !=======================================================================
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: lp1
  integer :: m
  integer :: n
  real :: x1
  real :: x11
  real :: x12
  real :: x2
  real :: x21
  real :: x22
  real :: y1
  real :: y11
  real :: y12
  real :: y2
  real :: y21
  real :: y22
  !
  do k = 1, lay
    if (inpt(1,k) < 950) then
      do i = il1, il2
        m =  inpt(i,k)
        n =  m + 1
        if (m <= 13) then
          if (m > 0) then
            y1  =  coef2u(1,m) + dt(i,k) * (coef2u(2,m) + dt(i,k) * &
                  (coef2u(3,m) + dt(i,k) * (coef2u(4,m) + &
                   dt(i,k) * coef2u(5,m))))
          else
            y1  =  0.0
          endif
          if (m < 13) then
            y2  =  coef2u(1,n) + dt(i,k) * (coef2u(2,n) + dt(i,k) * &
                  (coef2u(3,n) + dt(i,k) * (coef2u(4,n) + &
                   dt(i,k) * coef2u(5,n))))
          else
            y2  =  coef2u(1,13) + dt(i,k) * (coef2u(2,13) + dt(i,k) * &
                  (coef2u(3,13) + dt(i,k) * (coef2u(4,13) + &
                   dt(i,k) * coef2u(5,13))))
          end if
          x1    =  0.0
          x2    =  0.0
        else
          j     =  m - 13
          n     =  j + 1
          l     =  inptr(i,k)
          if (l < 1) then
            x1  =  coef1d(1,1,j) + dt(i,k) * (coef1d(2,1,j) + dt(i,k) * &
                  (coef1d(3,1,j) + dt(i,k) * (coef1d(4,1,j) + &
                   dt(i,k) * coef1d(5,1,j))))
            x2  =  coef1d(1,1,n) + dt(i,k) * (coef1d(2,1,n) + dt(i,k) * &
                  (coef1d(3,1,n) + dt(i,k) * (coef1d(4,1,n) + &
                   dt(i,k) * coef1d(5,1,n))))
            !
            y1  =  coef2d(1,1,j) + dt(i,k) * (coef2d(2,1,j) + dt(i,k) * &
                  (coef2d(3,1,j) + dt(i,k) * (coef2d(4,1,j) + &
                   dt(i,k) * coef2d(5,1,j))))
            y2  =  coef2d(1,1,n) + dt(i,k) * (coef2d(2,1,n) + dt(i,k) * &
                  (coef2d(3,1,n) + dt(i,k) * (coef2d(4,1,n) + &
                   dt(i,k) * coef2d(5,1,n))))
            !
          else if (l < 8) then
            lp1 =  l + 1
            x11 =  coef1d(1,l,j) + dt(i,k) * (coef1d(2,l,j) + dt(i,k) * &
                  (coef1d(3,l,j) + dt(i,k) * (coef1d(4,l,j) + &
                   dt(i,k) * coef1d(5,l,j))))
            x21 =  coef1d(1,l,n) + dt(i,k) * (coef1d(2,l,n) + dt(i,k) * &
                  (coef1d(3,l,n) + dt(i,k) * (coef1d(4,l,n) + &
                   dt(i,k) * coef1d(5,l,n))))
            !
            y11 =  coef2d(1,l,j) + dt(i,k) * (coef2d(2,l,j) + dt(i,k) * &
                  (coef2d(3,l,j) + dt(i,k) * (coef2d(4,l,j) + &
                   dt(i,k) * coef2d(5,l,j))))
            y21 =  coef2d(1,l,n) + dt(i,k) * (coef2d(2,l,n) + dt(i,k) * &
                  (coef2d(3,l,n) + dt(i,k) * (coef2d(4,l,n) + &
                   dt(i,k) * coef2d(5,l,n))))
            !
            x12 =  coef1d(1,lp1,j) + dt(i,k) * (coef1d(2,lp1,j) + &
                   dt(i,k) * (coef1d(3,lp1,j) + dt(i,k) * (coef1d(4,lp1,j) + &
                   dt(i,k) * coef1d(5,lp1,j))))
            x22 =  coef1d(1,lp1,n) + dt(i,k) * (coef1d(2,lp1,n) + &
                   dt(i,k) * (coef1d(3,lp1,n) + dt(i,k) * (coef1d(4,lp1,n) + &
                   dt(i,k) * coef1d(5,lp1,n))))
            !
            y12 =  coef2d(1,lp1,j) + dt(i,k) * (coef2d(2,lp1,j) + &
                   dt(i,k) * (coef2d(3,lp1,j) + dt(i,k) * (coef2d(4,lp1,j) + &
                   dt(i,k) * coef2d(5,lp1,j))))
            y22 =  coef2d(1,lp1,n) + dt(i,k) * (coef2d(2,lp1,n) + &
                   dt(i,k) * (coef2d(3,lp1,n) + dt(i,k) * (coef2d(4,lp1,n) + &
                   dt(i,k) * coef2d(5,lp1,n))))
            !
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1d(1,8,j) + dt(i,k) * (coef1d(2,8,j) + dt(i,k) * &
                  (coef1d(3,8,j) + dt(i,k) * (coef1d(4,8,j) + &
                   dt(i,k) * coef1d(5,8,j))))
            x2  =  coef1d(1,8,n) + dt(i,k) * (coef1d(2,8,n) + dt(i,k) * &
                  (coef1d(3,8,n) + dt(i,k) * (coef1d(4,8,n) + &
                   dt(i,k) * coef1d(5,8,n))))
            y1  =  coef2d(1,8,j) + dt(i,k) * (coef2d(2,8,j) + dt(i,k) * &
                  (coef2d(3,8,j) + dt(i,k) * (coef2d(4,8,j) + &
                   dt(i,k) * coef2d(5,8,j))))
            y2  =  coef2d(1,8,n) + dt(i,k) * (coef2d(2,8,n) + dt(i,k) * &
                  (coef2d(3,8,n) + dt(i,k) * (coef2d(4,8,n) + &
                   dt(i,k) * coef2d(5,8,n))))
          end if
        end if
        !
        taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * q(i,k) + &
                    (y1 + (y2 - y1) * dip(i,k)) * co2(i,k) ) * dp(i,k)
      end do ! loop 100

    else
      m =  inpt(1,k) - 1000
      n =  m + 1
      do i = il1, il2
        if (m <= 13) then
          if (m > 0) then
            y1  =  coef2u(1,m) + dt(i,k) * (coef2u(2,m) + dt(i,k) * &
                  (coef2u(3,m) + dt(i,k) * (coef2u(4,m) + &
                   dt(i,k) * coef2u(5,m))))
          else
            y1  =  0.0
          endif
          if (m < 13) then
            y2  =  coef2u(1,n) + dt(i,k) * (coef2u(2,n) + dt(i,k) * &
                  (coef2u(3,n) + dt(i,k) * (coef2u(4,n) + &
                   dt(i,k) * coef2u(5,n))))
          else
            y2  =  coef2u(1,13) + dt(i,k) * (coef2u(2,13) + dt(i,k) * &
                  (coef2u(3,13) + dt(i,k) * (coef2u(4,13) + &
                   dt(i,k) * coef2u(5,13))))
          end if
          x1    =  0.0
          x2    =  0.0
        else
          j     =  m - 13
          n     =  j + 1
          l     =  inptr(i,k)
          if (l < 1) then
            x1  =  coef1d(1,1,j) + dt(i,k) * (coef1d(2,1,j) + dt(i,k) * &
                  (coef1d(3,1,j) + dt(i,k) * (coef1d(4,1,j) + &
                   dt(i,k) * coef1d(5,1,j))))
            x2  =  coef1d(1,1,n) + dt(i,k) * (coef1d(2,1,n) + dt(i,k) * &
                  (coef1d(3,1,n) + dt(i,k) * (coef1d(4,1,n) + &
                   dt(i,k) * coef1d(5,1,n))))
            !
            y1  =  coef2d(1,1,j) + dt(i,k) * (coef2d(2,1,j) + dt(i,k) * &
                  (coef2d(3,1,j) + dt(i,k) * (coef2d(4,1,j) + &
                   dt(i,k) * coef2d(5,1,j))))
            y2  =  coef2d(1,1,n) + dt(i,k) * (coef2d(2,1,n) + dt(i,k) * &
                  (coef2d(3,1,n) + dt(i,k) * (coef2d(4,1,n) + &
                   dt(i,k) * coef2d(5,1,n))))
            !
          else if (l < 8) then
            lp1 =  l + 1
            x11 =  coef1d(1,l,j) + dt(i,k) * (coef1d(2,l,j) + dt(i,k) * &
                  (coef1d(3,l,j) + dt(i,k) * (coef1d(4,l,j) + &
                   dt(i,k) * coef1d(5,l,j))))
            x21 =  coef1d(1,l,n) + dt(i,k) * (coef1d(2,l,n) + dt(i,k) * &
                  (coef1d(3,l,n) + dt(i,k) * (coef1d(4,l,n) + &
                   dt(i,k) * coef1d(5,l,n))))
            !
            y11 =  coef2d(1,l,j) + dt(i,k) * (coef2d(2,l,j) + dt(i,k) * &
                  (coef2d(3,l,j) + dt(i,k) * (coef2d(4,l,j) + &
                   dt(i,k) * coef2d(5,l,j))))
            y21 =  coef2d(1,l,n) + dt(i,k) * (coef2d(2,l,n) + dt(i,k) * &
                  (coef2d(3,l,n) + dt(i,k) * (coef2d(4,l,n) + &
                   dt(i,k) * coef2d(5,l,n))))
            !
            x12 =  coef1d(1,lp1,j) + dt(i,k) * (coef1d(2,lp1,j) + &
                   dt(i,k) * (coef1d(3,lp1,j) + dt(i,k) * (coef1d(4,lp1,j) + &
                   dt(i,k) * coef1d(5,lp1,j))))
            x22 =  coef1d(1,lp1,n) + dt(i,k) * (coef1d(2,lp1,n) + &
                   dt(i,k) * (coef1d(3,lp1,n) + dt(i,k) * (coef1d(4,lp1,n) + &
                   dt(i,k) * coef1d(5,lp1,n))))
            !
            y12 =  coef2d(1,lp1,j) + dt(i,k) * (coef2d(2,lp1,j) + &
                   dt(i,k) * (coef2d(3,lp1,j) + dt(i,k) * (coef2d(4,lp1,j) + &
                   dt(i,k) * coef2d(5,lp1,j))))
            y22 =  coef2d(1,lp1,n) + dt(i,k) * (coef2d(2,lp1,n) + &
                   dt(i,k) * (coef2d(3,lp1,n) + dt(i,k) * (coef2d(4,lp1,n) + &
                   dt(i,k) * coef2d(5,lp1,n))))
            !
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1d(1,8,j) + dt(i,k) * (coef1d(2,8,j) + dt(i,k) * &
                  (coef1d(3,8,j) + dt(i,k) * (coef1d(4,8,j) + &
                   dt(i,k) * coef1d(5,8,j))))
            x2  =  coef1d(1,8,n) + dt(i,k) * (coef1d(2,8,n) + dt(i,k) * &
                  (coef1d(3,8,n) + dt(i,k) * (coef1d(4,8,n) + &
                   dt(i,k) * coef1d(5,8,n))))
            y1  =  coef2d(1,8,j) + dt(i,k) * (coef2d(2,8,j) + dt(i,k) * &
                  (coef2d(3,8,j) + dt(i,k) * (coef2d(4,8,j) + &
                   dt(i,k) * coef2d(5,8,j))))
            y2  =  coef2d(1,8,n) + dt(i,k) * (coef2d(2,8,n) + dt(i,k) * &
                  (coef2d(3,8,n) + dt(i,k) * (coef2d(4,8,n) + &
                   dt(i,k) * coef2d(5,8,n))))
          end if
        end if
        !
        taug(i,k) = ( (x1 + (x2 - x1) * dip(i,k)) * q(i,k) + &
                    (y1 + (y2 - y1) * dip(i,k)) * co2(i,k) ) *  dp(i,k)
      end do
    end if
  end do 

  return
end subroutine tlinehc5
!> \file
!> Calculation of optical depth for h2o and co2 (in 540-800 cm, line
!! contribution only), since this two gases are overlapped very strong,
!! we can not use the alternate sorting method, but use the direct
!! mapping method by sorting the two gases together. The pre-calculations
!! are done for 5 h2o/co2 ratios of .06,  .24,  1., 4., 16, then using
!! the interpolation to obtain results of other ratios. This method is
!! much more expensive than alternate sorting, we only applied it to the
!! lower atmosphere with INPT > 11, above this level we still use alternate sorting
!! The absorption coefficients are calculated at the temperature T for the 18 (not for GH)
!! pressure levels.
!!\n
!! First use the polynomial interpolation to get results at the
!! temperature T based on the pre-calculated results at 5 temperatures.
!! Then, do the linear interoplation for between two pressure levels,
!! X1 and X2 for h2o, Y1 and Y2 for co2, are the absoption coefficients
!! in unit cm^2/gram.
!!\n
!! The more details are in \cite Li2005a.
