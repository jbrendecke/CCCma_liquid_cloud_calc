!> \file lwtragh4.F90
!>\brief Longwave radiative transfer calculations for optically thick atmosphere
!!
!! @author Jiangnan Li
!
subroutine lwtragh4(fut, fdt, &
                    slwf, tauci, omci, taual, taug, bf, &
                    urbf, cld, em0t, bst, itile, cut, &
                    il1, il2, ilg, lay, lev, &
                    ntile)
  !
  !     * jun 02,2015 - m.lazare/ new version for gcm19:
  !     *               j.cole:   - add tiled radiation calculations
  !     *                           (ie "fut","fdt")
  !     * feb 11,2009 - j.cole.  previous version lwtragh4 for gcm15h
  !     *                        through gcm18:
  !     *                         - correct bug in specification of abse0.
  !     *                         - initialize clear-sky fx in top layer.
  !     *                         - change any work arrays having a middle
  !     *                           dimension of "4" to "2", as this is
  !     *                           what is used.
  !     * dec 05,2007 - m.lazare. previous version lwtragh3 for gcm15g:
  !     *                         - support added for emissivity<1.
  !     * nov 22/2006 - m.lazare. previous version lwtragh2 for gcm15f:
  !     *                         - bound of 1.e-10 used for 1-cld
  !     *                           instead of subtracting an epsilon.
  !     *                         - work arays for this routine now
  !     *                           are automatic arrays and their
  !     *                           space is not passed in.
  !     * aug 29,2003 - j.li.       previous version lwtragh up to gcm15e.
  implicit none
  real, intent(in) :: cut
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: ntile  !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  real, intent(in), dimension(ilg) :: slwf !< Input solar flux at model top level \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,lay) :: tauci !< Cloud optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: omci !< Cloud single sacttering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taual !< Aerosol optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taug !< Gaseous optical depth \f$[1]\f$
  real, intent(in), dimension(ilg,lev) :: bf !< Blackbody intensity integrated over each band at each level \f$[W/m^2/sr]\f$
  real, intent(in), dimension(ilg,lay) :: urbf !< Diffuse factor times the difference of log(BF) for two neighbor levels used in
  real, intent(in), dimension(ilg,lay) :: cld !< Cloud fraction\f$[1]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: fut !< Tiled upward infrared flux \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: fdt !< Tiled downward infrared flux\f$[W/m^2]\f$
  real, intent(in), dimension(ilg,ntile) :: em0t !< Tiled surface emisivity \f$[1]\f$
  real, intent(in), dimension(ilg,ntile) :: bst !< Tiled blackbody intensity integrated over each band at each level \f$[W/m^2/sr]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< Surface tile number \f$[1]\f$
  !
  !==================================================================
  !     in the g space with interval close 1 (very large optical depth)  !
  !     or in the case with cloud absorption is very small or the weight !
  !     of flux and cooling rate are very small. the cloud radiative     !
  !     process can be highly simplified. the absorption approximation   !
  !     method is used and cloud random and maximum overlap is           !
  !     considered, but cloud scattering and inhomogeneity are ignored.  !
  !     the exponential source planck function is used which is more     !
  !     accurate in the region above 200 mb in comparison with linear    !
  !     source function                                                  !
  !                                                                      !
  !     FYT and FX are used to deal with the cloud overlap, but it is not
  !     used when McICA is applied.
  !==================================================================
  !
  real :: abse0
  real :: cow
  real :: crtaul2
  real :: ctaul2
  real :: embs
  real :: epsd
  real :: epsu
  integer :: i
  integer :: k
  integer :: km1
  integer :: km2
  integer :: kp1
  integer :: m
  real :: rtaul1
  real :: taul1
  real :: taul2
  real :: ubeta
  !
  real, dimension(ilg,2,lev) :: fd
  real, dimension(ilg,2,lay) :: xu  !< the emission part in the upward flux transmission (li, 2002 jas p3302)\f$[0]\f$
  real, dimension(ilg,2,lay) :: xd  !< the emission part in the downward flux transmission
  real, dimension(ilg,2,lay) :: dtr !< direct transmission
  real, dimension(ilg,2,lev) :: fx  !< the same as fy but for the downward flux
  real, dimension(ilg,ntile,2,lev) :: fyt !< Variable description\f$[units]\f$
  !
  real, parameter :: ru = 1.6487213
  !
  !----------------------------------------------------------------------c
  !     initialization for first layer. calculate the downward flux in   c
  !     the second layer                                                 c
  !     combine the optical properties for the infrared,                 c
  !     1, aerosol + gas; 2, cloud + aerosol + gas.                      c
  !     fd (fu) is down (upward) flux                                    c
  !     the overlap between solar and infrared in 4 - 10 um is           c
  !     considered, slwf is the incoming solar flux                      c
  !     singularity for xd and xu has been considered as li jas 2002     c
  !----------------------------------------------------------------------c
  !
  do i = il1, il2
    fd(i,1,1)         =  slwf(i)
    fd(i,2,1)         =  slwf(i)
    fx(i,1,1)         =  slwf(i)
    fx(i,2,1)         =  slwf(i)
    !
    taul1             =  taual(i,1) + taug(i,1)
    rtaul1            =  taul1 * ru
    dtr(i,1,1)        =  exp ( - rtaul1)
    ubeta             =  urbf(i,1) / (taul1 + 1.e-20)
    epsd              =  ubeta + 1.0
    epsu              =  ubeta - 1.0
    !
    if (abs(epsd) > 0.001) then
      xd(i,1,1)       = (bf(i,2) - bf(i,1) * dtr(i,1,1)) / epsd
    else
      xd(i,1,1)       =  rtaul1 * bf(i,1) * dtr(i,1,1)
    end if
    if (abs(epsu) > 0.001) then
      xu(i,1,1)       = (bf(i,2) * dtr(i,1,1) - bf(i,1)) / epsu
    else
      xu(i,1,1)       =  rtaul1 * bf(i,2) * dtr(i,1,1)
    end if
    !
    fd(i,1,2)         =  fd(i,1,1) * dtr(i,1,1) + xd(i,1,1)
    !
    if (cld(i,1) < cut) then
      fx(i,1,2)       =  fd(i,1,2)
      fx(i,2,2)       =  fd(i,1,2)
      fd(i,2,2)       =  fd(i,1,2)
    else
      taul2           =  tauci(i,1) + taul1
      cow             =  1.0 - omci(i,1) / taul2
      ctaul2          =  cow * taul2
      crtaul2         =  ctaul2 * ru
      dtr(i,2,1)      =  exp ( - crtaul2)
      ubeta           =  urbf(i,1) / (ctaul2)
      epsd            =  ubeta + 1.0
      epsu            =  ubeta - 1.0
      !
      if (abs(epsd) > 0.001) then
        xd(i,2,1)     = (bf(i,2) - bf(i,1) * dtr(i,2,1)) / epsd
      else
        xd(i,2,1)     =  crtaul2 * bf(i,1) * dtr(i,2,1)
      end if
      if (abs(epsu) > 0.001) then
        xu(i,2,1)     = (bf(i,2) * dtr(i,2,1) - bf(i,1)) / epsu
      else
        xu(i,2,1)     =  crtaul2 * bf(i,2) * dtr(i,2,1)
      end if
      !
      fx(i,1,2)       =  fx(i,1,1) * dtr(i,1,1) + xd(i,1,1)
      fx(i,2,2)       =  fx(i,2,1) * dtr(i,2,1) + xd(i,2,1)
      fd(i,2,2)       =  fx(i,1,2) + &
                        cld(i,1) * (fx(i,2,2) - fx(i,1,2))
    end if
  end do ! loop 100
  !
  do k = 3, lev
    km1 = k - 1
    km2 = km1 - 1
    do i = il1, il2
      taul1           =  taual(i,km1) + taug(i,km1)
      rtaul1          =  taul1 * ru
      dtr(i,1,km1)    =  exp ( - rtaul1)
      ubeta           =  urbf(i,km1) / (taul1 + 1.e-20)
      epsd            =  ubeta + 1.0
      epsu            =  ubeta - 1.0
      !
      if (abs(epsd) > 0.001) then
        xd(i,1,km1)   = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
      else
        xd(i,1,km1)   =  rtaul1 * bf(i,km1) * dtr(i,1,km1)
      end if
      if (abs(epsu) > 0.001) then
        xu(i,1,km1)   = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
      else
        xu(i,1,km1)   =  rtaul1 * bf(i,k) * dtr(i,1,km1)
      end if
      !
      fd(i,1,k)       =  fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
      !
      if (cld(i,km1) < cut) then
        fd(i,2,k)     =  fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
        fx(i,1,k)     =  fd(i,2,k)
        fx(i,2,k)     =  fd(i,2,k)
      else
        taul2         =  tauci(i,km1) + taul1
        cow           =  1.0 - omci(i,km1) / taul2
        ctaul2        =  cow * taul2
        crtaul2       =  ctaul2 * ru
        dtr(i,2,km1)  =  exp ( - crtaul2)
        ubeta         =  urbf(i,km1) / (ctaul2)
        epsd          =  ubeta + 1.0
        epsu          =  ubeta - 1.0
        !
        if (abs(epsd) > 0.001) then
          xd(i,2,km1) = (bf(i,k) - bf(i,km1) * dtr(i,2,km1)) / epsd
        else
          xd(i,2,km1) =  crtaul2 * bf(i,km1) * dtr(i,2,km1)
        end if
        if (abs(epsu) > 0.001) then
          xu(i,2,km1) = (bf(i,k) * dtr(i,2,km1) - bf(i,km1)) / epsu
        else
          xu(i,2,km1) =  crtaul2 * bf(i,k) * dtr(i,2,km1)
        end if
        !
        if (cld(i,km1) <= cld(i,km2)) then
          fx(i,1,k)   = ( fx(i,2,km1) + (1.0 - cld(i,km2)) / &
                        (max(1.0 - cld(i,km1),1.e-10)) * &
                        (fx(i,1,km1) - fx(i,2,km1)) ) * &
                        dtr(i,1,km1) + xd(i,1,km1)
          fx(i,2,k)   =  fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
        else if (cld(i,km1) > cld(i,km2)) then
          fx(i,1,k)   =  fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
          fx(i,2,k)   = (fx(i,1,km1) + cld(i,km2) / cld(i,km1) * &
                        (fx(i,2,km1) - fx(i,1,km1))) * &
                        dtr(i,2,km1) + xd(i,2,km1)
        end if
        !
        fd(i,2,k)     =  fx(i,1,k) + cld(i,km1) * (fx(i,2,k) - &
                        fx(i,1,k))
      end if
    end do ! loop 200
  end do ! loop 250
  !
  do k = 1, lev
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          fdt(i,m,1,k) = fd(i,1,k)
          fdt(i,m,2,k) = fd(i,2,k)
        end if
      end do ! i
    end do ! m
  end do ! k

  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        embs             =  em0t(i,m) * bst(i,m)
        abse0            =  1.0 - em0t(i,m)
        fut(i,m,1,lev)   =  embs + abse0 * fdt(i,m,1,lev)
        fyt(i,m,1,lev)   =  embs + abse0 * fx(i,1,lev)
        fyt(i,m,2,lev)   =  embs + abse0 * fx(i,2,lev)
        fut(i,m,2,lev)   =  fyt(i,m,1,lev) + &
                           cld(i,lay) * &
                           (fyt(i,m,2,lev) - fyt(i,m,1,lev))
        !
        fut(i,m,1,lay)   =  fut(i,m,1,lev) * dtr(i,1,lay) + &
                           xu(i,1,lay)
        fyt(i,m,1,lay)   =  fyt(i,m,1,lev) * dtr(i,1,lay) + &
                           xu(i,1,lay)
        !
        if (cld(i,lay) < cut) then
          fyt(i,m,2,lay) =  fyt(i,m,2,lev) * dtr(i,1,lay) + &
                           xu(i,1,lay)
          fut(i,m,2,lay) =  fyt(i,m,1,lay)
        else
          fyt(i,m,2,lay) =  fyt(i,m,2,lev) * dtr(i,2,lay) + &
                           xu(i,2,lay)
          fut(i,m,2,lay) =  fyt(i,m,1,lay) + &
                           cld(i,lay) * &
                           (fyt(i,m,2,lay) - fyt(i,m,1,lay))
        end if
      end if
    end do ! i
  end do ! m
  !
  do k = lev - 2, 1, - 1
    kp1 = k + 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          fut(i,m,1,k) = fut(i,m,1,kp1) * dtr(i,1,k) + xu(i,1,k)
          !
          if (cld(i,k) < cut) then
            fut(i,m,2,k) = fut(i,m,2,kp1) * dtr(i,1,k) + &
                           xu(i,1,k)
            fyt(i,m,1,k) = fut(i,m,2,k)
            fyt(i,m,2,k) = fut(i,m,2,k)
          else
            if (cld(i,k) < cld(i,kp1)) then
              fyt(i,m,1,k) = (fyt(i,m,2,kp1) + (1.0 - cld(i,kp1)) &
                             / (1.0 - cld(i,k)) * (fyt(i,m,1,kp1) &
                             -  fyt(i,m,2,kp1)) ) * dtr(i,1,k) &
                             +  xu(i,1,k)
              fyt(i,m,2,k) = fyt(i,m,2,kp1) * dtr(i,2,k) + &
                             xu(i,2,k)
            else
              fyt(i,m,1,k) = fyt(i,m,1,kp1) * dtr(i,1,k) + &
                             xu(i,1,k)
              fyt(i,m,2,k) = (fyt(i,m,1,kp1) + cld(i,kp1) &
                             / cld(i,k) &
                             * (fyt(i,m,2,kp1) - fyt(i,m,1,kp1))) &
                             * dtr(i,2,k) + xu(i,2,k)
            end if
            !
            fut(i,m,2,k) = fyt(i,m,1,k) + &
                           cld(i,k) * (fyt(i,m,2,k) - fyt(i,m,1,k))
          end if
        end if
      end do ! i
    end do ! m
  end do ! k
  !
  return
end subroutine lwtragh4
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
