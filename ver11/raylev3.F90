!> \file raylev4.F90
!>\brief Compute Rayleigh scatter optical thickness for wavelength interval 1
!!
!! @author Jiangnan Li
!
subroutine raylev3 (taur, ig, dp, rmu, il1, il2, ilg, lay)
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
  real, intent(in), dimension(ilg) :: rmu !< Cosine of solar zenth angle \f$[1]\f$
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
  real, dimension(ilg,9) :: factormu
  real, dimension(9), parameter :: ri0 = & 
       [0.48040E-04, 0.82006E-04, 0.20968E-03, 0.56998E-03, 0.89894E-03, &
        0.98420E-03, 0.11055E-02, 0.13062E-02, 0.14980E-02]
  !=======================================================================
    do i = il1, il2
      factormu(i,1) = 9.97922e-01 -2.28904e-04/rmu(i) + 1.79688e-06*rmu(i)
      factormu(i,2) = 9.91919e-01 -5.78653e-04/rmu(i) + 1.03984e-05*rmu(i)
      factormu(i,3) = 9.32745e-01 -4.82633e-03/rmu(i) + 1.31798e-03*rmu(i)
      factormu(i,4) = 9.43895e-01 -7.29559e-03/rmu(i) + 4.41385e-03*rmu(i)
      factormu(i,5) = 9.99007e-01 -2.16808e-04/rmu(i) + 6.90153e-06*rmu(i)
      factormu(i,6) = 0.999285
      factormu(i,7) = 0.994034
      factormu(i,8) = 0.992955
      factormu(i,9) = 0.992025
    enddo
    do k = 1, lay
      do i = il1, il2
        taur(i,k) = ri0(ig) * factormu(i,ig) * dp(i,k)
      end do
    end do
  !
  return
end subroutine raylev3
!> \file
