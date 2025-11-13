!> \file
!>\brief Solar radiative transfer calculation with no subgrid-scale horizontal variability.
!!
!! @author Jiangnan Li and Jason Cole
!
subroutine swtran_mcica2(refl, tran, itile, &
                         ucumdtr, tran0, taua, taur, taug, &
                         tauoma, tauomga, f1, f2, tauc, &
                         tauomc, tauomgc, cld, rmu, c1, c2, &
                         albsurt, csalbt, nct, cut, &
                         il1, il2, ilg, lay, lev, ntile)
  !
  !     * sept 3/2018 - j. cole  promoted more variables to real*8
  !     *                        hopefully this will make code more robust.
  !
  !     * feb 10/2014 - j.li.    revised version for gcm18:
  !     *                        - compute the direct beam using the
  !     *                          unscaled cloud optical properties.
  !     * feb 20/2007 - j. cole: updated code to follow changes made between
  !                              swtran and swtran2
  !
  !     * apr 11/2005 - j. cole: used petri's hard work and modified to suit my version.
  !                              no longer has subcolumns, works on only one column per
  !                              gcm column.
  !     * may 31/2004 - p.raisanen:
  !
  !     this code version has been updated for mcica radiation calculations
  !
  !       - clear-sky reflectances and transmittances are calculated
  !         using mean profiles for the gcm column, but cloudy-sky
  !         reflectances and transmittances are calculated separately
  !         for nsub subcolumns and then averaged. it is assumed that
  !         there are equally many subcolumns for each gcm column i1...i2.
  !
  !       - treatment of cloud overlap and horizontal variability eliminated.
  !         cloud fraction for subcolumns (cldsub) is still included as
  !         input, but the code assumes it is either 0 or 1.
  !
  !       - security checks added to avoid values of x2 and x4 equal to 1.
  !
  !     * apr 25,2003 - j.li.
  !
  implicit none
  real, intent(in) :: cut
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: ntile  !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  real, intent(inout), dimension(ilg,2,lev) :: ucumdtr !< Direct transmission for multi-layers using unscale cloud properties \f$[1]\f$
  real, intent(inout), dimension(ilg) :: tran0 !< Downward solar flux at highest computation level \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: refl !< Layer reflectivity for all sky and clear sky (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: tran !< Layer transmission for all sky and clear sky (tiled) \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,ntile) :: albsurt !< All sky surface albedo (tiled) \f$[1]\f$
  real, intent(in), dimension(ilg,ntile) :: csalbt !< Clear sky surface albedo (tiled) \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taua !< Aerosol optical thickness \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taur !< Rayleigh optical thickness \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: taug !< Gaseous optical thickness \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauoma !< Aerosol optical thickness times aerosol single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauomga !< TAUOMA times asymmetry factor \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: f1 !< TAUOMGA times asymmetry factor \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: f2 !< TAUOMGC times single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauc !< Cloud optical thickness \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauomc !< TAUC times cloud single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: tauomgc !< TAUOMC times single scattering albedo \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: cld !< Cloud fraction \f$[1]\f$
  real, intent(in), dimension(ilg) :: rmu !< Diffuse reflection from model bottom level \f$[W/m^2]\f$
  real, intent(in), dimension(ilg) :: c1 !< Factor based on solar zenith \f$[1]\f$
  real, intent(in), dimension(ilg) :: c2 !< Factor based on solar zenith \f$[1]\f$
  integer, intent(in), dimension(ilg) :: nct !< Index of highest cloud level \f$[1]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< Index of surface tile \f$[1]\f$
  !
  ! Local variables
  !
  real :: dmm
  real :: dpp
  real :: fmm
  real :: fpp
  integer :: i
  integer :: k
  integer :: km1
  integer :: l
  integer :: lp1
  integer :: m
  real :: umm
  real :: upp
  !
  real, dimension(ilg,2,lev) :: cumdtr !<direct transmission for mult-layers\f$[0]\f$
  real*8 , dimension(ilg,2,lay) :: rdf   !<layer diffuse reflection\f$[w/m^2]\f$
  real*8 , dimension(ilg,2,lay) :: tdf   !<layer diffuse transmission\f$[w/m^2]\f$
  real*8 , dimension(ilg,2,lay) :: rdr   !<layer direct reflection\f$[w/m^2]\f$
  real*8 , dimension(ilg,2,lay) :: tdr   !<layer diffuse transmission\f$[w/m^2]\f$
  real*8 , dimension(ilg,2,lay) :: dtr   !<direct transmission\f$[0]\f$
  real*8 , dimension(ilg,2,lay) :: udtr  !<direct transmission using unscale
                                         ! cloud properties\f$[0]\f$
  real*8 , dimension(ilg,ntile,2,lev) :: rmdf !<diffuse reflection from model top leve
  real*8 , dimension(ilg,ntile,2,lev) :: tmdr !<direct transmission from model top level
  real*8 , dimension(ilg,ntile,2,lev) :: rmur !<direct reflection from model bottom level
  real*8 , dimension(ilg,ntile,2,lev) :: rmuf !<diffuse reflection from model bottom level
  !==================================================================
  !     delta-eddington approximation and adding process for clear and
  !     all sky, the adding method by coakley et al (1983). this code
  !     can deal with solar radiative transfer through atmosphere with
  !     cloud fraction either 1 or 0.
  !==================================================================

  !
  !     * scalars promoted to 64-bit in loop 200 to avoid pole singularity
  !     * from two-stream calculation in clear-sky.
  !
  real*8 :: extopt
  real*8 :: omars
  real*8 :: ssalb
  real*8 :: sf1
  real*8 :: sf
  real*8 :: tau1
  real*8 :: om1
  real*8 :: cow
  real*8 :: ssgas1
  real*8 :: cowg
  real*8 :: alamd
  real*8 :: u
  real*8 :: u2
  real*8 :: uu
  real*8 :: efun
  real*8 :: efun2
  real*8 :: x1
  real*8 :: rn
  real*8 :: x2
  real*8 :: y
  real*8 :: yx
  real*8 :: dm
  real*8 :: gscw
  real*8 :: appgm
  real*8 :: apmgm
  real*8 :: omarcs
  real*8 :: sf2
  real*8 :: tau2
  real*8 :: om2
  real*8 :: ssgas2
  real*8 :: x3
  real*8 :: x4
  real*8 :: x

  real*8 :: rdf8
  real*8 :: tdf8
  real*8 :: rdr8
  real*8 :: tdr8
  real*8 :: dtr8
  real*8 :: udtr8
  real*8 :: arg_alamd
  real*8 :: arg_dtr
  real*8 :: arg_udtr
  real*8 :: arg_efun
  real*8 , dimension(ilg) :: rmu8
  real*8 , dimension(ilg) :: c18
  real*8 , dimension(ilg) :: c28
  real*8 , dimension(ilg,lay) :: taua8
  real*8 , dimension(ilg,lay) :: taur8
  real*8 , dimension(ilg,lay) :: taug8
  real*8 , dimension(ilg,lay) :: tauoma8
  real*8 , dimension(ilg,lay) :: tauomga8
  real*8 , dimension(ilg,lay) :: f18
  real*8 :: f28 (ilg,lay)
  real*8 :: tauc8(ilg,lay)
  real*8 :: tauomc8(ilg,lay)
  real*8 :: tauomgc8 (ilg,lay)

  real*8, parameter :: r_0_5 = 0.5_8
  real*8, parameter :: r_one = 1.0_8
  real*8, parameter :: r_1_5 = 1.5_8
  real*8, parameter :: r_three = 3.0_8
  real*8, parameter :: r_zero = 0.0_8
  real*8, parameter :: r_mone = - 1.0_8
  real*8, parameter :: eps_10 = 1.0d-10
  real*8, parameter :: eps_20 = 1.0d-20
  real*8, parameter :: lowbnd = 1.0d-15

  !-------------------
  ! initializations
  !-------------------

  do k = 1, lev
    do i = il1,il2
      cumdtr(i,1:2,k) = r_zero
      ucumdtr(i,1:2,k) = r_zero
    end do ! loop 20
  end do ! loop 10

  do k = 1, lev
    do m = 1, ntile
      do i = il1,il2
        refl(i,m,1:2,k) = r_zero
        tran(i,m,1:2,k) = r_zero
      end do ! i
    end do ! m
  end do ! k

  do i = il1, il2
    rmu8(i) = real(rmu(i),8)
    c18(i)  = real(c1(i),8)
    c28(i)  = real(c2(i),8)
  end do

  do k = 1, lay
    do i = il1, il2
      taua8(i,k)    = real(taua(i,k),8)
      taur8(i,k)    = real(taur(i,k),8)
      taug8(i,k)    = real(taug(i,k),8)
      tauoma8(i,k)  = real(tauoma(i,k),8)
      tauomga8(i,k) = real(tauomga(i,k),8)
      tauc8(i,k)    = real(tauc(i,k),8)
      tauomc8(i,k)  = real(tauomc(i,k),8)
      tauomgc8(i,k) = real(tauomgc(i,k),8)
      f18(i,k)      = real(f1(i,k),8)
      f28(i,k)      = real(f2(i,k),8)
    end do ! i
  end do ! k
  !
  !----------------------------------------------------------------------
  !     combine the optical properties for solar,
  !     1, aerosol + rayleigh + gas; 2, cloud + aerosol + rayleigh + gas
  !     calculate the direct and diffuse reflection and transmission in
  !     the scattering layers using the delta-eddington method.
  !----------------------------------------------------------------------
  !
  !----------------------------------------------------------------------
  ! first, computation of clear-sky fluxes                   1
  !----------------------------------------------------------------------

  do k = 1, lay
    do i = il1, il2
      extopt            =  taua8(i,k) +  taur8(i,k) +  taug8(i,k)
      omars             =  tauoma8(i,k) + taur8(i,k)
      ssalb             =  omars / (extopt + eps_20)
      sf1               =  f18(i,k) / omars
      sf                =  ssalb * sf1
      tau1              =  extopt * (r_one - sf)
      om1               = (ssalb - sf) / (r_one - sf)
      cow               =  r_one - om1 + eps_10
      ssgas1            = (tauomga8(i,k) / omars - sf1) / (r_one - sf1)
      cowg              =  r_one - om1 * ssgas1
      arg_alamd         =  r_three * cow * cowg
      alamd             =  sqrt(arg_alamd)
      !
      arg_dtr           =  r_mone * tau1 / rmu8(i)
      arg_udtr          =  r_mone * extopt / rmu8(i)
      dtr8              =  exp(arg_dtr)
      udtr8             =  exp(arg_udtr)
      u                 =  r_1_5 * cowg / alamd
      u2                =  u + u
      uu                =  u * u
      arg_efun          =  r_mone * alamd * tau1
      efun              =  exp(arg_efun)
      efun2             =  efun * efun
      x1                = (uu - u2 + r_one) * efun2
      rn                =  r_one / (uu + u2 + r_one - x1)
      x2                =  alamd * rmu8(i)
      y                 =  r_one - x2 * x2
      yx                =  sign( max( abs(y), lowbnd), y)
      dm                =  om1 / yx
      gscw              =  ssgas1 * cow
      appgm             = (c18(i) + r_0_5 + gscw * (c18(i) + c28(i))) * dm
      apmgm             = (c18(i) - r_0_5 + gscw * (c18(i) - c28(i))) * dm
      rdf8              = (uu - r_one) * (r_one  - efun2) * rn
      tdf8              = (u2 + u2) * efun * rn
      rdr8              =  appgm * rdf8 +  apmgm * (tdf8 * dtr8 - r_one)
      tdr8              =  appgm * tdf8 + (apmgm * rdf8 - appgm + r_one) * dtr8

      rdf(i,1,k)        = real(rdf8,4)
      tdf(i,1,k)        = real(tdf8,4)
      rdr(i,1,k)        = real(rdr8,4)
      tdr(i,1,k)        = real(tdr8,4)
      dtr(i,1,k)        = real(dtr8,4)
      udtr(i,1,k)       = real(udtr8,4)
    end do
  end do ! loop 200

  !
  do i = il1, il2
    !
    !----------------------------------------------------------------------
    !     initialization for the first level.
    !----------------------------------------------------------------------
    !
    cumdtr(i,1,1)       =  tran0(i)
    ucumdtr(i,1,1)      =  tran0(i)
  end do ! loop 300

  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        !
        !----------------------------------------------------------------------
        !     initialization for the first level.
        !----------------------------------------------------------------------
        !
        tmdr(i,m,1,1)   =  tran0(i)
        rmdf(i,m,1,1)   =  1.0 - tran0(i)
        !
        !----------------------------------------------------------------------
        !     initialization for the ground layer.
        !----------------------------------------------------------------------
        !
        rmur(i,m,1,lev) =  csalbt(i,m)
        rmuf(i,m,1,lev) =  csalbt(i,m)
      end if
    end do ! i
  end do ! m
  !
  !----------------------------------------------------------------------
  !     add the layers downward from the second layer to the surface.
  !----------------------------------------------------------------------
  !
  do k = 2, lev
    km1 = k - 1
    do i = il1, il2
      cumdtr(i,1,k)     =  cumdtr(i,1,km1) *  dtr(i,1,km1)
      ucumdtr(i,1,k)    =  ucumdtr(i,1,km1) * udtr(i,1,km1)
    end do ! i
  end do ! k

  do k = 2, lev
    km1 = k - 1
    l = lev - k + 1
    lp1 = l + 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          dmm           =  tdf(i,1,km1) / (1.0 - rdf(i,1,km1) * rmdf(i,m,1,km1))
          fmm           =  rmdf(i,m,1,km1) * dmm
          tmdr(i,m,1,k) =  cumdtr(i,1,km1) * (tdr(i,1,km1) + rdr(i,1,km1) * fmm) + &
                          (tmdr(i,m,1,km1) - cumdtr(i,1,km1)) * dmm
          rmdf(i,m,1,k) =  rdf(i,1,km1) + tdf(i,1,km1) * fmm

          !
          !----------------------------------------------------------------------
          !     add the layers upward from one layer above surface to the level 1.
          !----------------------------------------------------------------------

          umm           =  tdf(i,1,l) / (1.0 - rdf(i,1,l) * rmuf(i,m,1,lp1))
          fmm           =  rmuf(i,m,1,lp1) * umm
          rmur(i,m,1,l) =  rdr(i,1,l) + dtr(i,1,l) * rmur(i,m,1,lp1) * umm + &
                          (tdr(i,1,l) - dtr(i,1,l)) * fmm
          rmuf(i,m,1,l) =  rdf(i,1,l) + tdf(i,1,l) * fmm
        end if
      end do ! loop 400
    end do ! loop 450
  end do ! loop 451
  !
  !----------------------------------------------------------------------
  !     add downward to calculate the resultant reflectances and
  !     transmittance at flux levels.
  !----------------------------------------------------------------------
  !
  do k = 1, lev
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          dmm           =  1.0 / (1.0 - rmuf(i,m,1,k) * rmdf(i,m,1,k))
          x             =  cumdtr(i,1,k) * rmur(i,m,1,k)
          y             =  tmdr(i,m,1,k) - cumdtr(i,1,k)
          tran(i,m,1,k) =  cumdtr(i,1,k) + (x * rmdf(i,m,1,k) + y) * dmm
          refl(i,m,1,k) = (x + y * rmuf(i,m,1,k)) * dmm
        end if
      end do ! loop 500
    end do ! loop 551
  end do ! loop 550

  !----------------------------------------------------------------------
  ! second, computation of all-sky fluxes
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! all-sky fluxes.
  !----------------------------------------------------------------------

  do k = 1, lay
    do i = il1, il2
      if (cld(i,k) < cut) then
        rdf(i,2,k)      =  rdf(i,1,k)
        tdf(i,2,k)      =  tdf(i,1,k)
        rdr(i,2,k)      =  rdr(i,1,k)
        tdr(i,2,k)      =  tdr(i,1,k)
        dtr(i,2,k)      =  dtr(i,1,k)
        udtr(i,2,k)     =  udtr(i,1,k)
      else
        extopt          =  tauc8(i,k) + taua8(i,k) + taur8(i,k) + taug8(i,k)
        omarcs          =  tauomc8(i,k) + taur8(i,k)
        ssalb           =  omarcs / extopt
        sf2             =  f28(i,k) / omarcs
        sf              =  ssalb * sf2
        tau2            =  extopt * (r_one - sf)
        om2             = (ssalb - sf) / (r_one - sf)
        cow             =  r_one - om2
        ssgas2          = (tauomgc8(i,k)/ omarcs - sf2) / (r_one - sf2)
        cowg            =  r_one - om2 * ssgas2
        arg_alamd       =  r_three * cow * cowg
        alamd           =  sqrt(arg_alamd)
        !
        arg_dtr         =  r_mone * tau2 / rmu8(i)
        arg_udtr        =  r_mone * extopt / rmu8(i)
        dtr8            =  exp(arg_dtr)
        udtr8           =  exp(arg_udtr)
        u               =  r_1_5 * cowg / alamd
        u2              =  u + u
        uu              =  u * u
        arg_efun        =  r_mone * alamd * tau2
        efun            =  exp(arg_efun)
        efun2           =  efun * efun
        x3              = (uu - u2 + r_one) * efun2
        rn              =  r_one / (uu + u2 + r_one - x3)
        x4              =  alamd * rmu8(i)
        y               =  r_one - x4 * x4
        yx              =  sign( max( abs(y), lowbnd), y)
        dm              =  om2 / yx
        gscw            =  ssgas2 * cow
        appgm           = (c18(i) + r_0_5 + gscw * (c18(i) + c28(i))) * dm
        apmgm           = (c1(i) - r_0_5 +  gscw * (c18(i) - c28(i))) * dm

        rdf8            = (uu - r_one) * (r_one - efun2) * rn
        tdf8            = (u2 + u2) * efun * rn
        rdr8            =  appgm * rdf8 + apmgm * (tdf8 * dtr8 - r_one)
        tdr8            =  appgm * tdf8 + (apmgm * rdf8 - appgm + r_one) * dtr8

        rdf(i,2,k)      = real(rdf8,4)
        tdf(i,2,k)      = real(tdf8,4)
        rdr(i,2,k)      = real(rdr8,4)
        tdr(i,2,k)      = real(tdr8,4)
        dtr(i,2,k)      = real(dtr8,4)
        udtr(i,2,k)     = real(udtr8,4)

      end if
    end do
  end do ! loop 600
  !
  do i = il1, il2
    !
    !----------------------------------------------------------------------
    !     initialization for the first level (level 1).
    !----------------------------------------------------------------------
    !
    cumdtr(i,2,1)       =  tran0(i)
    ucumdtr(i,2,1)      =  tran0(i)

  end do ! loop 700

  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        !
        !----------------------------------------------------------------------
        !     initialization for the first level (level 1).
        !----------------------------------------------------------------------
        !

        tmdr(i,m,2,1)   =  tran0(i)
        rmdf(i,m,2,1)   =  1.0 - tran0(i)
        !
        !----------------------------------------------------------------------
        !     initialization for the ground layer.
        !----------------------------------------------------------------------
        !
        rmur(i,m,2,lev) =  albsurt(i,m)
        rmuf(i,m,2,lev) =  albsurt(i,m)
      end if
    end do ! i
  end do ! m
  !----------------------------------------------------------------------
  !     add the layers downward from the second layer to the surface.
  !----------------------------------------------------------------------
  !
  do k = 2, lev
    km1 = k - 1
    do i = il1, il2
      if (nct(i) == lev) then
        cumdtr(i,2,k)    =  cumdtr(i,1,k)
        ucumdtr(i,2,k)   =  ucumdtr(i,1,k)
      else
        if (k <= nct(i)) then
          cumdtr(i,2,k)  =  cumdtr(i,1,k)
          ucumdtr(i,2,k) =  ucumdtr(i,1,k)
        else
          cumdtr(i,2,k)  =  cumdtr(i,2,km1) * dtr(i,2,km1)
          ucumdtr(i,2,k) =  ucumdtr(i,2,km1) * udtr(i,2,km1)
        end if
      end if
    end do ! i
  end do ! k

  do k = 2, lev
    km1 = k - 1
    l = lev - k + 1
    lp1 = l + 1
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          if (nct(i) <= lay) then
            if (k <= nct(i)) then
              tmdr(i,m,2,k) =  tmdr(i,m,1,k)
              rmdf(i,m,2,k) =  rmdf(i,m,1,k)
            else
              dpp           =  tdf(i,2,km1) / (1.0 - rmdf(i,m,2,km1) * rdf(i,2,km1))
              fpp           =  rmdf(i,m,2,km1) * dpp
              tmdr(i,m,2,k) =  cumdtr(i,2,km1) * (tdr(i,2,km1) + rdr(i,2,km1) * fpp) + &
                              (tmdr(i,m,2,km1) - cumdtr(i,2,km1)) * dpp
              rmdf(i,m,2,k) =  rdf(i,2,km1) + tdf(i,2,km1) * fpp
            end if
            !
            !----------------------------------------------------------------------
            !     add the layers upward from one layer above surface to the level 1.
            !----------------------------------------------------------------------
            !
            upp             =  tdf(i,2,l) / (1.0 - rmuf(i,m,2,lp1) * rdf(i,2,l))
            fpp             =  rmuf(i,m,2,lp1) * upp
            rmur(i,m,2,l)   =  rdr(i,2,l) + dtr(i,2,l) * rmur(i,m,2,lp1) * upp + &
                              (tdr(i,2,l) - dtr(i,2,l)) * fpp
            rmuf(i,m,2,l)   =  rdf(i,2,l) + tdf(i,2,l) * fpp
          else
            tmdr(i,m,2,k)   =  1.0
            rmdf(i,m,2,k)   =  0.0
            rmur(i,m,2,l)   =  0.0
            rmuf(i,m,2,l)   =  0.0
          end if
        end if
      end do ! loop 800
    end do ! loop 850
  end do ! loop 851

  !
  !----------------------------------------------------------------------
  !     add downward to calculate the resultant reflectances and
  !     transmittance at flux levels.
  !----------------------------------------------------------------------
  !

  do k = 1, lev
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          if (nct(i) <= lay) then
            dpp             =  1.0 / (1.0 - rmuf(i,m,2,k) * rmdf(i,m,2,k))
            x               =  cumdtr(i,2,k) * rmur(i,m,2,k)
            y               =  tmdr(i,m,2,k) - cumdtr(i,2,k)
            tran(i,m,2,k)   =  cumdtr(i,2,k) + (x * rmdf(i,m,2,k) + y) * dpp
            refl(i,m,2,k)   = (x + y * rmuf(i,m,2,k)) * dpp
          else
            tran(i,m,2,k)   =  tran(i,m,1,k)
            refl(i,m,2,k)   =  refl(i,m,1,k)
          end if
        end if
      end do ! loop 900
    end do ! loop 951
  end do ! loop 950

  return
end subroutine swtran_mcica2
!> \file
!> Compute solar radiative transfer using delta-Eddington approximation and adding process for clear and
!! all sky, the adding method by \cite Coakley1983. This code can deal with solar radiative transfer
!! through atmosphere with cloud fraction either 1 or 0.
