!> \file raylei.F90
!>\brief Compute Rayleigh scatter optical thickness for wavelength intervals 2, 3 and 4
!!
!! @author Jiangnan Li
!
subroutine raylei (taur, ib, dp, il1, il2, ilg, lay)
  !
  !     * apr 25,2003 - j.li.
  !----------------------------------------------------------------------c
  !     rayleigh scattering for bands2-bands4, near infrared region      c
  !                                                                      c
  !     taur: rayleigh optical depth                                     c
  !     dp:   air mass path for a model layer (exlained in raddriv).     c
  !----------------------------------------------------------------------c
  implicit none
  integer, intent(in) :: ib
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: taur !< Raylegh scattering optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  !
  !==================================================================
  integer :: i
  integer :: ibm1
  integer :: k
  real, dimension(3), parameter :: ri = &   !<extinction coefficient \f$[cm^2/gram]\f$
                                        [1.649498e-05, 1.7997e-06, 0.13586e-07]
  !=======================================================================
  !
  ibm1 = ib - 1
  do k = 1, lay
    do i = il1, il2
      taur(i,k) = 1.08 * ri(ibm1) * dp(i,k)
    end do
  end do ! loop 100
  !
  return
end subroutine raylei
!> \file
!> Compute the Rayleigh optical properties for solar wavelength intervals 2, 3 and 4,
!! i.e., the intervals in the near infrared.
