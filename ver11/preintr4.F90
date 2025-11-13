!> \file preintr4.F90
!>\brief Determine the parameters needed to interpolate in H2O/CO2 ratio for gas optical properties
!!
!! @author Jiangnan Li
subroutine preintr4 (inpr, dir, q, co2, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI.     FOR NEW CKD WITH 6 VALUES OF H2O/CO2
  !     * feb 09,2009 - j.li.     new version for gcm15h:
  !     *                         - co2 passed directly, thus no need
  !     *                           for "trace" common block.
  !     * apr 21,2008 - l.solheim. previous version printr2 for gcm15g:
  !     *                          - cosmetic change to add threadprivate
  !     *                            for common block "trace", in support
  !     *                            of "radforce" model option.
  !     * apr 25,2003 - j.li.  previous version preintr for gcm15e/gcm15f.
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: dir !< Interpretation factor for H2O/CO2 ratio \f$[1]\f$
  real, intent(in), dimension(ilg,lay)    :: q   !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)    :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  integer, intent(inout), dimension(ilg,lay) :: inpr !< Index of standard input H2O/CO2 ratio for each model layer \f$[1]\f$
  !==================================================================
  !     this subroutine determines the interpretion points for the ratio of h2o and co2.
  !     STAND shows 8 standard values of the ratios of h2o/co2, the
  !     absorption coeffients of the comined h2o and co2 are
  !     pre-calculated for the 8 ratios, general results are the
  !     interpolations based on the 8 values.
  !==================================================================
  integer :: i
  integer :: j
  integer :: jendr
  integer :: k
  integer :: l
  integer :: lp1
  real, dimension(ilg,lay) :: rhc !< Ratio of H2O / CO2 \f$[1]\f$
  !
  !8 pre-set values of h20/co2 ratio,results are interpolated based on the 8 values
  real, dimension(8), parameter :: standr = [0.02, 0.2, .39, 0.75, 1.2, 1.8, 4., 16.]
  !=======================================================================
  jendr  =  8
  do k = 1, lay
    !
    do i = il1, il2
      inpr(i,k)   =  0
      rhc(i,k)    =  q(i,k) / max(co2(i,k), 1.0e-10)
    end do
    !
    do j = 1, jendr
      do i = il1, il2
        if (rhc(i,k) > standr(j)) then
          inpr(i,k) =  inpr(i,k) + 1
        end if
      end do
    end do
    !
    do i = il1, il2
      l   =  inpr(i,k)
      lp1 =  l + 1
      if (l >= 1 .and. l < 8) then
        dir(i,k)  = (rhc(i,k) - standr(l)) / (standr(lp1) - standr(l))
      else
        !
        !----------------------------------------------------------------------c
        !     dir is not used with values of {0,6} in tlinehc, but we          c
        !     initialize here to avoid problems with nan when used             c
        !     in multitasking mode.                                            c
        !----------------------------------------------------------------------c
        !
        dir(i,k)  =  0.0
      end if
    end do
    !
  end do
  !
  return
end subroutine preintr4
!> \file
!> Determines the interpretion points for the ratio of H2O and CO2.
!! STAND shows 5 standard values of the ratios of H2O/CO2, the
!! absorption coeffients of the comined H2O and CO2 are
!! pre-calculated for the 5 ratios, general results are the
!! interpolations based on the 5 values.
