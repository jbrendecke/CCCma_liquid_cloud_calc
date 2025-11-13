!> \file
!>\brief Longwave radiative transfer calculation with no subgrid-scale horizontal variability.
!!
!! @author Jiangnan Li and Jason Cole
!
subroutine lwtran_mcica3(fut, fdt, &
                         slwf, tauc, omc, gc, fl, taual, &
                         taug, bf, bst, urbf, dbf, em0t, &
                         cld, nct, cut, &
                         itile, il1, il2, ilg, lay, lev, ntile)
  !
  !     * jun 02,2015 - m.lazare/ new version for gcm19:
  !     *               j.cole:   - add tiled radiation calculations
  !     *                           (ie "fut","fdt"), under control of
  !     *                           "itilrad".
  !     * feb 11,2009 - j.cole.  previous version lwtran_mcica2 for gcm15h
  !     *                        through gcm18:
  !     * feb 20/2007 - j. cole. previous version lwtran_mcica for gcm15g:
  !     *                        - update to match lwtran2, automatic arrays.
  !     *                        - used petri's hard work and modified to
  !     *                          suit my version. no longer has subcolumns,
  !     *                          works on only one subcolumn per gcm
  !     *                          column.
  !     * may 31/2004 - p.raisanen:
  !
  !     this code version has been updated for mcica radiation calculations
  !
  !       - clear-sky fluxes are calculated using mean profiles for the gcm
  !         column, but cloudy-sky fluxes are calculated separately
  !         for nsub subcolumns and then averaged. it is assumed that
  !         there are equally many subcolumns for each gcm column il1...il2.
  !
  !       - treatment of cloud overlap and horizontal variability eliminated.
  !         cloud fraction for subcolumns (cldsub) is still included as
  !         input, but the code assumes it is either 0 or 1.
  !
  !
  !     * feb 04,2004 - j.li, m.lazare.
  !
  implicit none
  real, intent(in) :: cut
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: ntile  !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,ntile,2,lev) :: fut !< Tiled upward infrared flux \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: fdt !< Tiled downward infrared flux \f$[W/m^2]\f$
  real, intent(in), dimension(ilg) :: slwf !< Input solar flux at model top level \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,lay) :: tauc !< Cloud optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: omc !< Cloud single sacttering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: gc !< Cloud symmetry factor \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: fl !< Cloud optical depth scaling factor \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taual !< Aerosol optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taug !< Gaseous optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lev) :: bf !< Blackbody intensity integrated over each
  real, intent(in), dimension(ilg,lay) :: urbf !< Diffuse factor times the difference of
  real, intent(in), dimension(ilg,lay) :: dbf !< Difference of bf for two neighbor levels
  real, intent(in), dimension(ilg,lay) :: cld !< Cloud fraction \f$[1]\f$
  real, intent(in), dimension(ilg,ntile) :: bst !< Tiled blackbody intensity integrated over
  real, intent(in), dimension(ilg,ntile) :: em0t !< Tiled surface emisivity \f$[1]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< Surface tile number \f$[1]\f$
  integer, intent(in), dimension(ilg) :: nct !< Highest cloud level \f$[1]\f$
  !
  real :: cow
  real :: dtr2
  real :: embk
  real :: epsd
  real :: epsu
  integer :: i
  integer :: k
  integer :: km1
  integer :: kp1
  integer :: m
  integer :: maxc
  real :: rtaul1
  real :: sf
  real :: ssalb
  real :: tau2
  real :: taudtr
  real :: taul1
  real :: taul2
  real :: ubeta
  real :: w
  real :: wgrcow
  real :: x
  real :: zeta
  !
  !     * internal work arrays.
  !
  real, dimension(ilg,lay) :: emisw       !<factor for emissivity\f$[0]\f$
  real, dimension(ilg,lay) :: scatbk      !<back scattering factor\f$[0]\f$
  real, dimension(ilg,lay) :: scatfw      !<forward scattering factor\f$[0]\f$
  real, dimension(ilg,ntile) :: embst     !<emisivity * BST\f$[W/M^2/SR]\f$
  real, dimension(ilg,ntile) :: abse0t    !<factor for tilted emissivity\f$[0]\f$
  real, dimension(ilg,2,lev) :: fd        !<downward infrared flux\f$[w/m^2]\f$
  real, dimension(ilg,2,lay) :: scatsm !<internal scattering factor\f$[0]\f$
  real, dimension(ilg,2,lay) :: taum   !<sacled cloud optical depth\f$[0]\f$
  real, dimension(ilg,2,lay) :: xd     !<the emission part in the downward flux transmission\f$[0]\f$
  real, dimension(ilg,2,lay) :: xu     !<the emission part in the upward flux transmission\f$[0]\f$
  real, dimension(ilg,2,lay) :: dtr    !<direct transmission\f$[0]\f$
  !
  real, parameter :: ru = 1.6487213
  !=======================================================================
  !
  !
  !=======================================================================
  !     Calculation of longwave radiative transfer using absorption
  !     approximation, but including corrections for scattering
  !     (unperturbed + backward scattering effect +  forward scattering
  !      effect + internal scattering effect  (Li JAS 2002)
  !     This code is for McICA, the cloud fraction is either 1 or 0 in the
  !     calculation, therefore the cloud overlap effect is not considered.
  !     Also, since the subgrid cloud effect has been considered in McICA
  !     the scaling for optical depth due to horizontal variability is
  !     not considered.
  !=======================================================================
  !
  !-------------------
  ! initializations
  !-------------------
  !
  do k = 1, lev
    do i = il1,il2
      fd(i,:,k) = 0.
    end do
  end do
  do k = 1,lev
    do m = 1, ntile
      do i = il1, il2
        fdt(i,m,:,k) = 0.0
        fut(i,m,:,k) = 0.0
      end do ! i
    end do ! m
  end do ! k

  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        embst(i,m)  = em0t(i,m) * bst(i,m)
        abse0t(i,m) = 1.0 - em0t(i,m)
      end if
    end do ! i
  end do ! m
  !
  !----------------------------------------------------------------------
  !     calculate the downward fluxes first without scattering
  !     combine the optical properties for the infrared,
  !     1, aerosol + gas; 2, cloud + aerosol + gas.
  !     fd (fu) is down (upward) flux.
  !     gaussian integration and diffusivity factor, ru (li jas 2000)
  !     above maxc, exponential source function is used
  !     below maxc, linear source function is used
  !----------------------------------------------------------------------
  !
  !----------------------------------------------------------------------
  ! computation of clear-sky downward fluxes using exponential source
  ! function
  !----------------------------------------------------------------------

  maxc = lev

  do i = il1, il2
    fd(i,1,1) =  slwf(i)
    maxc      =  min (nct(i), maxc)
  end do
  !
  do k = 2, lev
    km1 = k - 1
    do i = il1, il2
      taul1                =  taual(i,km1) + taug(i,km1)
      rtaul1               =  taul1 * ru
      dtr(i,1,km1)         =  exp ( - rtaul1)
      ubeta                =  urbf(i,km1) / (taul1 + 1.e-20)
      epsd                 =  ubeta + 1.0
      epsu                 =  ubeta - 1.0
      !
      if (abs(epsd) > 0.001) then
        xd(i,1,km1)        = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
      else
        xd(i,1,km1)        =  rtaul1 * bf(i,km1) * dtr(i,1,km1)
      end if
      !
      if (abs(epsu) > 0.001) then
        xu(i,1,km1)        = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
      else
        xu(i,1,km1)        =  rtaul1 * bf(i,k) * dtr(i,1,km1)
      end if
      !
      fd(i,1,k)            =  fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
    end do
  end do
  do k = 1,lev
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          fdt(i,m,1,k)     =  fd(i,1,k)
        end if
      end do ! i
    end do ! m
  end do ! k

  !----------------------------------------------------------------------
  ! computation of clear-sky upward fluxes
  !----------------------------------------------------------------------

  do m = 1, ntile
    do i = il1,il2
      if (itile(i,m) > 0) then
        fut(i,m,1,lev)     =  embst(i,m) + abse0t(i,m) * fdt(i,m,1,lev)
      end if
    end do ! i
  end do ! m

  do k = lev - 1, 1, - 1
    kp1 = k + 1
    do m = 1, ntile
      do i = il1,il2
        if (itile(i,m) > 0) then
          fut(i,m,1,k)     =  fut(i,m,1,kp1) * dtr(i,1,k) + xu(i,1,k)
        end if
      end do ! i
    end do ! m
  end do ! k

  !----------------------------------------------------------------------
  ! all-sky downward fluxes above the highest cloud top
  !----------------------------------------------------------------------

  do i = il1,il2
    fd(i,2,1)              =  slwf(i)
  end do

  do k = 2, lev
    do i = il1,il2
      if (k <= nct(i)) then
        fd(i,2,k)          =  fd(i,1,k)
      end if
    end do
  end do

  !----------------------------------------------------------------------
  !     add the layers downward from the highest cloud layer to the
  !     surface. determine the xu for the upward path.
  !     using exponential source function also for all sky flux in cloud
  !     free layers.
  !----------------------------------------------------------------------
  !
  do k = maxc + 1, lev
    km1 = k - 1
    do i = il1, il2
      if (k > nct(i)) then
        if (cld(i,km1) < cut) then
          fd(i,2,k)        =  fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
        else
          taul2            =  tauc(i,km1) + taug(i,km1)
          ssalb            =  omc(i,km1) / taul2
          sf               =  ssalb * fl(i,km1)
          w                = (ssalb - sf) / (1.0 - sf)
          cow              =  1.0 - w
          taum(i,1,km1)    = (cow * taul2 * (1.0 - sf) + taual(i,km1)) * ru
          zeta             =  dbf(i,km1) / taum(i,1,km1)
          tau2             =  taum(i,1,km1) + taum(i,1,km1)
          !
          dtr(i,2,km1)     =  exp( - taum(i,1,km1))
          dtr2             =  exp( - tau2)
          emisw(i,km1)     =  zeta * (1.0 - dtr(i,2,km1))
          embk             =  dtr2 - 1.0
          !
          xd(i,2,km1)      =  bf(i,k) - bf(i,km1) * dtr(i,2,km1) - emisw(i,km1)
          xu(i,2,km1)      =  bf(i,km1) - bf(i,k) * dtr(i,2,km1) + emisw(i,km1)
          !
          wgrcow           =  w * gc(i,km1) / cow
          taudtr           =  taum(i,1,km1) * dtr(i,2,km1)
          scatfw(i,km1)    =  wgrcow * taudtr
          scatbk(i,km1)    =  0.5 * wgrcow * (dtr2 - 1.0)
          !
          x                =  wgrcow * (2.0 * emisw(i,km1) + &
                                (0.5 * embk - taudtr) * zeta)
          scatsm(i,1,km1)  =  - scatbk(i,km1) * bf(i,k) - &
                                 scatfw(i,km1) * bf(i,km1) - x
          scatsm(i,2,km1)  =  - scatbk(i,km1) * bf(i,km1) - &
                                 scatfw(i,km1) * bf(i,k) + x
          fd(i,2,k)        =  fd(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
        end if
      end if
    end do
  end do
  do k = 1,lev
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          fdt(i,m,2,k)     =  fd(i,2,k)
        end if
      end do ! i
    end do ! m
  end do ! k

  !----------------------------------------------------------------------
  ! computation of all-sky upward fluxes
  !----------------------------------------------------------------------
  ! initialization for surface
  !----------------------------------------------------------------------
  k = lev - 1
  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        fut(i,m,2,lev)     =  embst(i,m) + abse0t(i,m) * fdt(i,m,2,lev)
      end if
    end do ! i
  end do ! m

  !----------------------------------------------------------------------
  !     add the layers upward from the first layer to maxc
  !     scattering effect for upward path is included
  !----------------------------------------------------------------------
  !
  do k = lev - 1, maxc, - 1
    kp1 = k + 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          if (k >= nct(i)) then
            if (cld(i,k) < cut) then
              fut(i,m,2,k) =  fut(i,m,2,kp1) * dtr(i,1,k) + xu(i,1,k)
            else
              fut(i,m,2,k) =  fut(i,m,2,kp1) * &
                             (dtr(i,2,k) + scatfw(i,k)) + &
                              fdt(i,m,2,k) * scatbk(i,k) + &
                              scatsm(i,2,k) + xu(i,2,k)
            end if
          end if
        end if
      end do ! i
    end do ! m
  end do ! k
  !
  !----------------------------------------------------------------------
  !     add the layers upward above the highest cloud to the toa, no
  !     scattering
  !----------------------------------------------------------------------
  !
  do k = lev - 1, 1, - 1
    kp1 = k + 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          if (kp1 <= nct(i)) then
            fut(i,m,2,k)   =  fut(i,m,2,kp1) * dtr(i,1,k) + xu(i,1,k)
          end if
        end if
      end do ! i
    end do ! m
  end do ! k
  !
  !----------------------------------------------------------------------
  !     scattering effect for downward path in from maxc to the surface
  !----------------------------------------------------------------------
  !
  do k = maxc + 1, lev
    km1 = k - 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          if (km1 >= nct(i)) then
            if (cld(i,km1) < cut) then
              fdt(i,m,2,k) =  fdt(i,m,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
            else
              fdt(i,m,2,k) =  fdt(i,m,2,km1) * &
                             (dtr(i,2,km1) + scatfw(i,km1)) + &
                              fut(i,m,2,k) * scatbk(i,km1) + &
                              scatsm(i,1,km1) + xd(i,2,km1)
            end if
          end if
        end if
      end do ! i
    end do ! m
  end do ! k
  !

  return
end subroutine lwtran_mcica3
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
