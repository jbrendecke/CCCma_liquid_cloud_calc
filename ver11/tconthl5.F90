!> \file
!>\brief Compute optical thickness due to water vapor continuum in 540-800 cm\f$^{-1}\f$
!!
!! @author Jiangnan Li
!
subroutine tconthl5 (taug, coef1, coef2, s, dp, dip, dir, dt, &
                     inptr, inpt, mcont, il1, il2, ilg, lay)
  !
  !     * feb 9,2009  - j.li.      new version for gcm15h:
  !     *                          - remove unused trace common block.
  !     * apr 21/2008 - l.solheim. previous version tconthl3 for gcm15g:
  !     *                          - cosmetic change to add threadprivate
  !     *                            for common block "trace", in support
  !     *                            of "radforce" model option.
  !     * jun 21/2006 - j.li.     previous version tconthl2 for gcm15f:
  !     *                         - new scheme to properly account for the
  !     *                           h2o/c)2 ratio.
  !     * jun 19/2006 - m.lazare. previous version tconthl1 for gcm15e:
  !     *                         - implement rpn fix for inpt.
  !     *                         - use variable instead of constant
  !     *                           in intrinsics such as "max",
  !     *                           so that can compile in 32-bit mode
  !     *                           with real  .
  !     * apr 25,2003 - j.li.     previous version tconthl for gcm15d.
  !----------------------------------------------------------------------
  !     water vapor continuum for 540-800 cm-1. different from tcontl,
  !     variation of mass mixing ratio for h2o and co2 is consider.
  !
  !     taug:  gaseous optical depth
  !     s:     input h2o mixing ratio for each layer
  !     dp:    air mass path for a model layer (exlained in raddriv).
  !     dip:   interpretation factor for pressure between two
  !            neighboring standard input data pressure levels
  !     dir:  interpretation factor for mass ratio of h2o / co2 between
  !           two neighboring standard input ratios
  !     dt:    layer temperature - 250 k
  !     inpr: number of the ratio level for the standard 5 ratios
  !     inpt:  number of the level for the standard input data pressures
  !     mcont: the highest level for water vapor continuum calculation
  !----------------------------------------------------------------------
  !
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: mcont
  !
  real, intent(inout), dimension(ilg,lay) :: taug !< Gas optical thickness \f$[1]\f$
  real, intent(in), dimension(5,8,4) :: coef1 !< Polynomial coefficients for h2o self continuum
                                              !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                              !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                              !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(5,8,4) :: coef2 !< Polynomial coefficients for h2o foreign continuum
                                              !! \f$[cm^2/gram,cm^2/gram K^{-1},cm^2/gram
                                              !! K^{-2},cm^2/gram K^{-3},cm^2/gram K^{-4}]\f$
                                              !! \f$[0,K^{-1},K^{-2},K^{-3},K^{-4}]\f$

  real, intent(in), dimension(ilg,lay) :: s !< H2O mixing ratio  \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between adjacent standard pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dir !< Interpretation factor for H2O/CO2 ratio \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f$
  integer, intent(in), dimension(ilg,lay) :: inptr !< Index of the slected standard input H2O/CO2 ratio \f$[1]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1  !
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
  !==================================================================
  !     Water vapor continuum for 540-800 cm-1. different from tcontl,
  !     variation of mass mixing ratio of h2o/co2 is considered.
  !     coef1 is the coefficient for self, coef2 is the coefficient for
  !     foreign, in COEF1(5,5,4) and COEF2(5,5,4), the first index is for
  !     the polynomial fit, the second is for the different h2o/co2 rate,
  !     and the third is for vertical level.
  !==================================================================
  !-----------------------------------------------------------------------
  do k = mcont, lay
    if (inpt(1,k) < 950) then
      do i = il1, il2
        m =  inpt(i,k) - 16
        n =  m + 1
        if (m >= 1) then
          l   =  inptr(i,k)
          if (l < 1) then
            x1  =  coef1(1,1,m) + dt(i,k) * (coef1(2,1,m) + dt(i,k) * (coef1(3,1,m) + &
                   dt(i,k) * (coef1(4,1,m) + dt(i,k) * coef1(5,1,m))))
            x2  =  coef1(1,1,n) + dt(i,k) * (coef1(2,1,n) + dt(i,k) * (coef1(3,1,n) + &
                   dt(i,k) * (coef1(4,1,n) + dt(i,k) * coef1(5,1,n))))
            !
            y1  =  coef2(1,1,m) + dt(i,k) * (coef2(2,1,m) + dt(i,k) * (coef2(3,1,m) + &
                   dt(i,k) * (coef2(4,1,m) + dt(i,k) * coef2(5,1,m))))
            y2  =  coef2(1,1,n) + dt(i,k) * (coef2(2,1,n) + dt(i,k) * (coef2(3,1,n) + &
                   dt(i,k) * (coef2(4,1,n) + dt(i,k) * coef2(5,1,n))))
            !
          else if (l < 8) then
            lp1 =  l + 1
            x11 =  coef1(1,l,m) + dt(i,k) * (coef1(2,l,m) + dt(i,k) * (coef1(3,l,m) + &
                   dt(i,k) * (coef1(4,l,m) + dt(i,k) * coef1(5,l,m))))
            x21 =  coef1(1,l,n) + dt(i,k) * (coef1(2,l,n) + dt(i,k) * (coef1(3,l,n) + &
                   dt(i,k) * (coef1(4,l,n) + dt(i,k) * coef1(5,l,n))))
            !
            y11 =  coef2(1,l,m) + dt(i,k) * (coef2(2,l,m) + dt(i,k) * (coef2(3,l,m) + &
                   dt(i,k) * (coef2(4,l,m) + dt(i,k) * coef2(5,l,m))))
            y21 =  coef2(1,l,n) + dt(i,k) * (coef2(2,l,n) + dt(i,k) * (coef2(3,l,n) + &
                   dt(i,k) * (coef2(4,l,n) + dt(i,k) * coef2(5,l,n))))
            !
            x12 =  coef1(1,lp1,m) + dt(i,k) * (coef1(2,lp1,m) + dt(i,k) * (coef1(3,lp1,m) + &
                   dt(i,k) * (coef1(4,lp1,m) + dt(i,k) * coef1(5,lp1,m))))
            x22 =  coef1(1,lp1,n) + dt(i,k) * (coef1(2,lp1,n) + dt(i,k) * (coef1(3,lp1,n) + &
                   dt(i,k) * (coef1(4,lp1,n) + dt(i,k) * coef1(5,lp1,n))))
            !
            y12 =  coef2(1,lp1,m) + dt(i,k) * (coef2(2,lp1,m) + dt(i,k) * (coef2(3,lp1,m) + &
                   dt(i,k) * (coef2(4,lp1,m) + dt(i,k) * coef2(5,lp1,m))))
            y22 =  coef2(1,lp1,n) + dt(i,k) * (coef2(2,lp1,n) + dt(i,k) * (coef2(3,lp1,n) + &
                   dt(i,k) * (coef2(4,lp1,n) + dt(i,k) * coef2(5,lp1,n))))
            !
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1(1,8,m) + dt(i,k) * (coef1(2,8,m) + dt(i,k) * (coef1(3,8,m) + &
                   dt(i,k) * (coef1(4,8,m) + dt(i,k) * coef1(5,8,m))))
            x2  =  coef1(1,8,n) + dt(i,k) * (coef1(2,8,n) + dt(i,k) * (coef1(3,8,n) + &
                   dt(i,k) * (coef1(4,8,n) + dt(i,k) * coef1(5,8,n))))
            y1  =  coef2(1,8,m) + dt(i,k) * (coef2(2,8,m) + dt(i,k) * (coef2(3,8,m) + &
                   dt(i,k) * (coef2(4,8,m) + dt(i,k) * coef2(5,8,m))))
            y2  =  coef2(1,8,n) + dt(i,k) * (coef2(2,8,n) + dt(i,k) * (coef2(3,8,n) + &
                   dt(i,k) * (coef2(4,8,n) + dt(i,k) * coef2(5,8,n))))
          end if
          !
          taug(i,k) =  taug(i,k) + ((x1 - y1 + (x2 - x1 - y2 + y1) * dip(i,k)) * 1.608 * s(i,k) + &
                       y1 + (y2 - y1) * dip(i,k)) * s(i,k) * dp(i,k)
        end if
      end do ! loop 100
    else
      j =  inpt(1,k) - 1000
      m =  j - 16
      n =  m + 1
      do i = il1, il2
        if (m >= 1) then
          l   =  inptr(i,k)
          if (l < 1) then
            x1  =  coef1(1,1,m) + dt(i,k) * (coef1(2,1,m) + dt(i,k) * (coef1(3,1,m) + &
                   dt(i,k) * (coef1(4,1,m) + dt(i,k) * coef1(5,1,m))))
            x2  =  coef1(1,1,n) + dt(i,k) * (coef1(2,1,n) + dt(i,k) * (coef1(3,1,n) + &
                   dt(i,k) * (coef1(4,1,n) + dt(i,k) * coef1(5,1,n))))
            !
            y1  =  coef2(1,1,m) + dt(i,k) * (coef2(2,1,m) + dt(i,k) * (coef2(3,1,m) + &
                   dt(i,k) * (coef2(4,1,m) + dt(i,k) * coef2(5,1,m))))
            y2  =  coef2(1,1,n) + dt(i,k) * (coef2(2,1,n) + dt(i,k) * (coef2(3,1,n) + &
                   dt(i,k) * (coef2(4,1,n) + dt(i,k) * coef2(5,1,n))))
            !
          else if (l < 8) then
            lp1 =  l + 1
            x11 =  coef1(1,l,m) + dt(i,k) * (coef1(2,l,m) + dt(i,k) * (coef1(3,l,m) + &
                   dt(i,k) * (coef1(4,l,m) + dt(i,k) * coef1(5,l,m))))
            x21 =  coef1(1,l,n) + dt(i,k) * (coef1(2,l,n) + dt(i,k) * (coef1(3,l,n) + &
                   dt(i,k) * (coef1(4,l,n) + dt(i,k) * coef1(5,l,n))))
            !
            y11 =  coef2(1,l,m) + dt(i,k) * (coef2(2,l,m) + dt(i,k) * (coef2(3,l,m) + &
                   dt(i,k) * (coef2(4,l,m) + dt(i,k) * coef2(5,l,m))))
            y21 =  coef2(1,l,n) + dt(i,k) * (coef2(2,l,n) + dt(i,k) * (coef2(3,l,n) + &
                   dt(i,k) * (coef2(4,l,n) + dt(i,k) * coef2(5,l,n))))
            !
            x12 =  coef1(1,lp1,m) + dt(i,k) * (coef1(2,lp1,m) + dt(i,k) * (coef1(3,lp1,m) + &
                   dt(i,k) * (coef1(4,lp1,m) + dt(i,k) * coef1(5,lp1,m))))
            x22 =  coef1(1,lp1,n) + dt(i,k) * (coef1(2,lp1,n) + dt(i,k) * (coef1(3,lp1,n) + &
                   dt(i,k) * (coef1(4,lp1,n) + dt(i,k) * coef1(5,lp1,n))))
            !
            y12 =  coef2(1,lp1,m) + dt(i,k) * (coef2(2,lp1,m) + dt(i,k) * (coef2(3,lp1,m) + &
                   dt(i,k) * (coef2(4,lp1,m) + dt(i,k) * coef2(5,lp1,m))))
            y22 =  coef2(1,lp1,n) + dt(i,k) * (coef2(2,lp1,n) + dt(i,k) * (coef2(3,lp1,n) + &
                   dt(i,k) * (coef2(4,lp1,n) + dt(i,k) * coef2(5,lp1,n))))
            !
            x1  =  x11 + (x12 - x11) * dir(i,k)
            x2  =  x21 + (x22 - x21) * dir(i,k)
            y1  =  y11 + (y12 - y11) * dir(i,k)
            y2  =  y21 + (y22 - y21) * dir(i,k)
          else
            x1  =  coef1(1,8,m) + dt(i,k) * (coef1(2,8,m) + dt(i,k) * (coef1(3,8,m) + &
                   dt(i,k) * (coef1(4,8,m) + dt(i,k) * coef1(5,8,m))))
            x2  =  coef1(1,8,n) + dt(i,k) * (coef1(2,8,n) + dt(i,k) * (coef1(3,8,n) + &
                   dt(i,k) * (coef1(4,8,n) + dt(i,k) * coef1(5,8,n))))
            y1  =  coef2(1,8,m) + dt(i,k) * (coef2(2,8,m) + dt(i,k) * (coef2(3,8,m) + &
                   dt(i,k) * (coef2(4,8,m) + dt(i,k) * coef2(5,8,m))))
            y2  =  coef2(1,8,n) + dt(i,k) * (coef2(2,8,n) + dt(i,k) * (coef2(3,8,n) + &
                   dt(i,k) * (coef2(4,8,n) + dt(i,k) * coef2(5,8,n))))
          end if
          !
          taug(i,k) =  taug(i,k) + ((x1 - y1 + (x2 - x1 - y2 + y1) * dip(i,k)) * 1.608 * s(i,k) + &
                       y1 + (y2 - y1) * dip(i,k)) * s(i,k) * dp(i,k)
        end if
      end do
    end if
  end do
  !
  return
end subroutine tconthl5
!> \file
!>  Water vapor continuum for 540-800 cm\f$^{-1}\f$. Different from tcontl,
!! variation of mass mixing ratio of H2O/CO2 is considered.
!! coef1 is the coefficient for self, coef2 is the coefficient for
!! foreign, in COEF1(5,5,4) and COEF2(5,5,4), the first index is for
!! the polynomial fit, the second is for the different H2O/CO2 rate,
!! and the third is for vertical level.
