!> \file raylev4.F90
!>\brief Compute Rayleigh scatter optical thickness for wavelength interval 1
!!
!! @author Jiangnan Li
!
subroutine raylev3 (taur, ig, dp, rmu3, il1, il2, ilg, lay)
  !
  !     * Jun 18,2021 - J.LI.  FOR NEW CKD, INCREASE SOLR SUBBAND IN UV
  !     * dec 05,2007 - j.li.  new version for gcm15g:
  !     *                      - revised data for ri0,ri2.
  !     * apr 25,2003 - j.li.  previous version raylev up through gcm15f.
  !----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lay) :: taur !< Raylegh scattering optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of a layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg) :: rmu3 !< Factor related to solar zenth angle \f$[1]\f$
  !==================================================================
  !     rayleigh scattering for each sub-band in bands1, visible region
  !     taur is the optical depth rayleigh scattering for a given layer
  !     for uvc (35700 - 50000 cm^-1), since the optical depth of o3 and
  !     o2 are very large, rayleigh scattering effect is neglected, it
  !     is shown even for 10% o3 amount of the standard atmo, the
  !     rayleigh scattering for uvc still can be neglected.
  !     for par and uva, since their spectral ranges are very wide, small
  !     errors could occur for large zenith angle, slightly adjustment
  !     is needed, this does mean the rayleigh optical depth is related
  !     solar zenith angle for multiple scattering process in swtran.
  !==================================================================
  !
  integer :: i
  integer :: k
  real, dimension(9), parameter :: ri0 = & 
       [0.67642E-04, 0.19855E-03, 0.42255E-03, 0.64011E-03, 0.86476E-03, &
        0.10037E-02, 0.11081E-02, 0.12575E-02, 0.14229E-02]
  real, dimension(6), parameter :: ri2 = & 
       [0.20000E-05, 0.12000E-04, 0.95000E-05, 0.20000E-04, 0.80000E-05, -.70000E-04]
  !=======================================================================
  if (ig <= 6) then
    do k = 1, lay
      do i = il1, il2
        taur(i,k) = (ri0(ig) - ri2(ig) * rmu3(i)) * dp(i,k)
      end do
    end do
  else
    do k = 1, lay
      do i = il1, il2
        taur(i,k) =  ri0(ig)  * dp(i,k)
      end do
    end do       
  endif
  !
  return
end subroutine raylev3
!> \file
