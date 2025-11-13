!> \file
!> \brief Brief description of the routine purpose.
!!
!! @author Routine author name(s)
!
subroutine initspecsampl2(nsample_sw,nsample_lw,nbs,nbl, &
                          maxng, ivers, iradforce, max_sam)
  !
  ! apr 30/2012  - j.cole. new version for gcm16:
  !                        add extra check, which requires passing
  !                        in "IRADFORCE" and "MAX_SAM".
  ! june 5, 2006 - j.cole. previous version initspecsampl.
  !
  ! update subroutine so that the variable vers defines the level of noise
  ! to use: 1 -> high    (clds)
  !         2 -> medium  (spec)
  !         3 -> low     (ref)
  !
  ! the computational time increases dramatically going from 1 to 3.  version
  ! 1 should (hopefully) be nearly as fast as the current code while 3 should
  ! be just able to complete one month in the wallclock time limit.
  !
  ! apr. 14, 2005 - j. cole
  !
  ! this code generates the spectral sampling for each ig correlated
  ! k-distribution integration point within each ib wavelength interval.
  ! the sampling was determined by tests by petri raisanen.
  !
  ! for each interval ib the spectral sampling is generated first for
  ! the kgs points and then the kgsgh points.  i.e., for the first band
  ! there are 7 values, 6 for kgs and 1 for kgsgh.  since for the kgsgh
  ! points clouds are ignored due to large gas optical thicknesses, the
  ! number of spectral samples is always 1.
  !
  ! there is more advanced sampling to be done based on location, i.e.,
  ! over land, over ocean, over snow, minimized error in atm. heating,
  ! surface forcing, etc..

  implicit none

  integer, intent(in) :: nbs   !< Number of wavelength intervals for solar radiative transfer \f$[unitless]\f$
  integer, intent(in) :: nbl   !< Number of wavelength intervals for thermal radiative transfer \f$[unitless]\f$
  integer, intent(in) :: maxng !< Variable description\f$[units]\f$
  integer, intent(in) :: ivers !< Variable description\f$[units]\f$
  integer, intent(in) :: iradforce !< Variable description\f$[units]\f$
  integer, intent(in) :: max_sam !< Variable description\f$[units]\f$

  integer, intent(out) , dimension(nbs,maxng,2) :: nsample_sw !< NSAMPLE_*(:,:,1) is for GH=.FALSE. and NSAMPLE_*(:,:,2) is for GH = .TRUE.\f$[units]\f$
  integer, intent(out) , dimension(nbl,maxng,2) :: nsample_lw !< NSAMPLE_*(:,:,1) is for GH=.FALSE. and NSAMPLE_*(:,:,2) is for GH = .TRUE.\f$[units]\f$
  !==================================================================
  ! physical (adjustable) parameters
  !
  ! define and document here any adjustable parameters.
  ! this should be variable described using the doxygen format above as
  ! well as a description of its minimum/default/maximum.
  !
  ! here is an example,
  !
  ! real :: beta !< This is the adjustable factor for computing liquid cloud effective radius \f$[]\f$
  !           !! It is compute differently when using bulk or PAM aerosols.
  !           !! For bulk aerosols its minimum/default/maximum is (1.0/1.3/1.5).
  !==================================================================

  integer :: ib
  integer :: ig
  integer :: igh

  ! set all to 1 spectral sample initially
  do igh = 1, 2
    do ig = 1, maxng
      do ib = 1, nbs
        nsample_sw(ib,ig,igh) = 1
      end do
      do ib = 1, nbl
        nsample_lw(ib,ig,igh) = 1
      end do ! ib
    end do ! ig
  end do  ! igh

  if (ivers == 1) then
    ! do nothing since we will use 1 sample

  else if (ivers == 2) then

    ! set the spectral samples so that a total of 45 calculations
    ! are made in the shortwave that minimize error for land and atm.

    nsample_sw(1,1:3,1) = 3
    nsample_sw(1,4:6,1) = 2

    nsample_sw(2,1:4,1) = 3
    nsample_sw(2,5,1)   = 2

    nsample_sw(3,1,1)   = 3
    nsample_sw(3,2:7,1) = 2

    ! set the spectral samples so that a total of 71 calculations
    ! are made in the longwave that minimize error for the atm.

    nsample_lw(3,1,1)   = 2
    nsample_lw(4,1,1)   = 3
    nsample_lw(4,2:4,1) = 2
    nsample_lw(5,1,1)   = 3
    nsample_lw(5,2,1)   = 2
    nsample_lw(6,1:2,1) = 3
    nsample_lw(7,1:4,1) = 3
    nsample_lw(8,1:4,1) = 3
    nsample_lw(8,5,1)   = 2

  else
     write( * , * ) 'UNKNOWN VALUE FOR IVERS ',ivers
     write( * , * ) 'MUST USE 1, 2, OR 3 !! '
  end if
  !
  !     * extra checks.
  !
  if (iradforce /= 0) then
    if (sum(nsample_lw) > max_sam .or. &
        sum(nsample_lw) > max_sam) then
  !    call xit('INITSPECSAMPL2/RADDRIV8', - 999)
    end if
  end if

  return
end
!> \file
!> This is an example of adding text at the end of the routine.
!! Your detailed description of the routine can be put here, including scientific
!! numerics, and another other important information.
!! Doxygen should be able to translate LaTeX and it is possible to include
!! references using "\cite", for example, \cite vonSalzen2013.
!! Equations can be included as well, as inline equations \f$ F=ma \f$,
!! or in the equation environment \n
!! (NOTE that HTML will not number the equation but it will in LaTeX),
!! \f{equation}{
!!  F_1=ma_1
!! \f}
