!> \file sw_cld_props4.F90
!>\brief Compute cloud optical properties at solar wavelengths
!!
!! @author Jiangnan Li and Jason Cole
!
subroutine sw_cld_props4(cldg, taucsg, tauomc, tauomgc, tauomgc_4str, f2, wrkag, wrkbg, & ! output
                         dz, cld, clw, cic, radeqvw, radeqvi, taua_g, tauoma_g,         & ! input
                         tauomga_g, tauomga_g_4str, f1_g, anu, cldt, rmug, eta, isun,   & ! input
                         mcica, cut, lengath, il1, il2, ilg, lay, ib, nstrm)              ! control
  !
  !     * nov 2021    - j.li.     add four-stream cloud optical property
  !     * feb 10/2015 - j.li.     new version for gcm18+:
  !     *                         - revised ice cloud optical properties
  !     *                           from yang et al (2012).
  !     * jun 14/2013 - j. li.    previous version sw_cld_props3 for gcm17+:
  !                               - add semi-direct effect.
  !     * feb 10/2009 - m.lazare. previous version sw_cld_props2 for gcm15h
  !     *                         through gcm16).
  !     *                         - calculate optical properties regardless
  !     *                           of cloud layer water content, it no
  !     *                           cutoff to zero below certain threshold.
  !     * oct 22/2007 - j.cole/   previous version sw_cld_props for gcm15g:
  !     *               m.lazare. compute cloud s/w optical properties
  !     *                         required for radiation. it is composed
  !     *                         of three parts:
  !     *                          1. assigns "cld" to "cldg".
  !     *                          2. calculation of {taucs,omcs,gcs}
  !     *                             migrated from clouds15 routines in
  !     *                             physics (separated into s/w and l/w
  !     *                             separate components).
  !     *                          3. calculation of resulting optical
  !     *                             properties arrays required for
  !     *                             radiation.
  !     *                          4. calculates wrka and wrkb as input
  !     *                             to ensuing surface albedo
  !     *                             calculation.
  !     *                         note that this means {taucs,omcs,gcs}
  !     *                         merely become internal work arrays.
  !
  use datcldop4, only :  aws, bws, cws, ais, bis, cis, som, sg
  implicit none
  !
  !     * output data.
  !
  real, intent(out) , dimension(ilg,lay)   :: cldg!< cloud fraction \f$[1]\f$
  real, intent(out) , dimension(ilg,lay)   :: taucsg!< cloud optical depth \f$[1]\f$
  real, intent(out) , dimension(ilg,lay)   :: tauomc!< Cloud optical depth times single scattering albedo \f$[1]\f$
  real, intent(out) , dimension(ilg,lay)   :: tauomgc!< TAUOMC times asymmetry factor \f$[1]\f$
  real, intent(out) , dimension(ilg,lay,2) :: tauomgc_4str!< TAUOMC times asymmetry factor \f$[1]\f$
  real, intent(out) , dimension(ilg,lay)   :: f2!< TAUOMGC times asymmetry factor \f$[1]\f$

  real, intent(out) , dimension(ilg)       :: wrkag!< Vertically integrated aerosol optical thickness \f$[1]\f$
  real, intent(out) , dimension(ilg)       :: wrkbg!< Vertically integrated aerosol absorption thickness \f$[1]\f$
  !
  !     * input data.
  !
  real, intent(in) , dimension(ilg,lay)   :: dz   !< Layer geometeric thickness\f$[m]\f$
  real, intent(in) , dimension(ilg,lay)   :: cld !< Total vertically projected cloud fraction \f$[1]\f$
  real, intent(in) , dimension(ilg,lay)   :: clw !< Liquid cloud water content \f$[gram/m^3]\f$
  real, intent(in) , dimension(ilg,lay)   :: cic !< Ice cloud water content \f$[gram/m^3]\f$
  real, intent(in) , dimension(ilg,lay)   :: radeqvw !< Liquid cloud effective radius\f$[um]\f$
  real, intent(in) , dimension(ilg,lay)   :: radeqvi !< Ice cloud effective radius\f$[um]\f$
  real, intent(in) , dimension(ilg,lay)   :: taua_g !< Aerosol optical thickness \f$[1]\f
  real, intent(in) , dimension(ilg,lay)   :: tauoma_g !< Aerosol optical thickness times single scattering albedo \f$[1]\f$
  real, intent(in) , dimension(ilg,lay)   :: tauomga_g !< TAUOMA times asymmetry factor \f$[1]\f$
  real, intent(in) , dimension(ilg,lay,2) :: tauomga_g_4str !< TAUOMA times asymmetry factor \f$[1]\f$
  real, intent(in) , dimension(ilg,lay)   :: f1_g !< TAUOMA times asymmetry factor \f$[1]\f$
  real, intent(in) , dimension(ilg,lay)   :: anu !< Relative horizontal variability assuming \f$\Gamma\f$-distribution \f$[1]\f$
  real, intent(in) , dimension(ilg,lay)   :: eta !< Volume fraction of black carbon in liquid cloud droplet \f$[m^3/m^3]\f$

  real, intent(in) , dimension(ilg)       :: cldt !< Total vertically projected cloud fraction \f$[1]\f$
  real, intent(in) , dimension(ilg)       :: rmug !< 1/(cosine of solar zenith) \f$[1]\f$

  integer, intent(in) , dimension(ilg)    :: isun !< Mapping of sunlit points along longitude circle \f$[1]\f$

  real, intent(in) :: cut !< Threshold minimum cloud fraction \f$[1]\f$
  
  integer, intent(in) :: mcica   !< Switch to use McICA or determinstic radiative transfer (0=deterministic, 1=McICA) \f$[unitless]\f$
  integer, intent(in) :: lengath !< Number of sunlit atmospheric columns <\f$[1]\f$
  integer, intent(in) :: il1   !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2   !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg   !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay   !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: ib    !< Solar band number \f$[1]\f$
  integer, intent(in) :: nstrm !< Switch to 2 or 4 stream radiative transfer \f$[unitless]\f$
  !==================================================================
  !     Cloud solar optical properties. It contains three parts:
  !     1. parameterization for solar optical properties based on effective
  !     radius and cloud liquid (ice) water contents, the effective
  !     variance is fixed, which shown in Dobbie et al (1999 JGR) and
  !     Yang et al (2012 JGR)
  !     2. BC semi-direct effect, the parameterization is based on ETA,
  !     the details are shown in Li et al. (2011 JGR)
  !     3. We can choose McICA or 4 column cloud method for radiative
  !     transfer, if the 4 column cloud method is used, the clould optical
  !     depth needs to be scaled down. For the 4 column cloud method, the
  !     physics is shown in Li et al (2005 QJRMS).
  !     4. Parameterization for ice cloud:
  !     Ice cloud optical property based on yang et al.(2012) jas
  !     and a series of papers by group of u a&m and u wisconsin
  !
  !     dg is the effective diameter
  !     taucs,omcs,gcs(taucl,omcl,gcl): optical depth,single
  !     scattering albedo,asymmetry factor for solar(infrared).
  !     radeqvw: effective radius(in micrometer) for water cloud
  !     radeqvi: effective radius(in micrometer) for ice cloud
  !     dg: geometry length for ice cloud
  !     wcdw(wcdi): liquid water(ice) content(in gram / m^3)
  !     wclw(wcli): liquid water(ice) path length(in gram / m^2)
  !     ccld: cloud fraction
  !==================================================================
  !
  !     * local data.
  !
  real, dimension(ilg) :: vi_taucs
  !
  real  :: wcdw
  real  :: wcdi
  real  :: wclw
  real  :: wcli
  real  :: rew1
  real  :: rew2
  real  :: rew3
  real  :: dg
  real  :: dg2
  real  :: dg3
  real  :: tausw
  real  :: omsw
  real  :: gsw, gsw2, gsw3, gsw4
  real  :: tausi
  real  :: omsi
  real  :: gsi, gsi2, gsi3, gsi4
  real  :: y1
  real  :: y2
  real  :: y12
  real  :: omsm
  real  :: gsm, gsm2, gsm3, gsm4
  real  :: seta
  real  :: x1
  real  :: x2
  real  :: dlt
  real  :: c20
  real  :: taucsi
  real  :: omcs
  real  :: gcs, gcs2, gcs3, gcs4
  real  :: x
  real  :: x_cas
  real  :: vi_taucsm
  real  :: rmugfac
  real  :: anufac
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: lp1
  !
  real, parameter, dimension(6) :: eta_tab = [1.e-08, 1.e-07, 1.e-06, 1.e-05, 1.e-04, 1.e-03]
  !
  !----------------------------------------------------------------------
  !     cloud radiative properties for radiation.
  !     taucsi, omcs, gcs (taucl, omcl, gcl): optical depth, single
  !     scattering albedo, asymmetry factor for solar (infrared).
  !     radeqvw: effective radius(in micrometer) for water cloud
  !     radeqvi: effective radius(in micrometer) for ice cloud
  !     dg: geometry length for ice cloud
  !     wcdw (wcdi): liquid water (ice) content (in gram / m^3)
  !     wclw (wcli): liquid water (ice) path length (in gram / m^2)
  !     cld/cldg: cloud fraction
  !     parameterization for water cloud:
  !     dobbie, etc. 1999, jgr, 104, 2067-2079
  !     lindner, t. h. and j. li., 2000, j. clim., 13, 1797-1805.
  !     parameterization for ice cloud:
  !     ice cloud optical property based on yang et al. (2012) jas
  !     and a serirs papers by group of u a &m and u wisconsin
  !     dg is the effective diameter
  !----------------------------------------------------------------------
  !
  do i = 1, lengath
    vi_taucs(i) = 0.0
    wrkag(i)    = 0.0
    wrkbg(i)    = 0.0
  end do
  !
  do k = 1, lay
    do i = 1, lengath
      j = isun(i)
      cldg(i,k) = cld(j,k)
      if (cld(j,k) <= cut) then
        x_cas        =  0.
        taucsi       =  0.
        vi_taucs(i)  =  0.
        taucsg(i,k)  =  0.
        tauomc(i,k)  =  0.
        tauomgc(i,k) =  0.
        tauomgc_4str(i,k,:) =  0.
        f2(i,k)      =  0.
      else
        wcdw  =  clw(j,k)
        wcdi  =  cic(j,k)
        wclw  =  wcdw * dz(j,k)
        wcli  =  wcdi * dz(j,k)

        if (nstrm >= 2) then
           rew1 = min(max(radeqvw(j,k), 5.0), 20.0)
        else ! nstrm == 0
           rew1 = radeqvw(j,k)
        end if
        rew2  =  rew1 * rew1
        rew3  =  rew2 * rew1
        dg    =  2.0 * min(max(radeqvi(j,k), 5.0), 60.0)
        dg2   =  dg  * dg
        dg3   =  dg2 * dg
        !
        !----------------------------------------------------------------------c
        !     add semi-direct effect inside cloud                              c
        !----------------------------------------------------------------------c
        !
        if (eta(j,k) > 1.e-03 .or. eta(j,k) < 1.e-08) then
          omsm =  0.0
          gsm  =  0.0
          gsm2 =  0.0
          gsm3 =  0.0
          gsm4 =  0.0
        else
          seta =  log10(eta(j,k)) + 9.
          l    =  int(seta)
          lp1  =  l + 1
          x1   =  som(1,ib,l) + som(2,ib,l) * rew1 + som(3,ib,l) * rew2
          y1   =  sg(1,ib,l) + sg(2,ib,l) * rew1 + sg(3,ib,l) * rew2
          if (lp1 <= 6) then
            x2 =  som(1,ib,lp1) + som(2,ib,lp1) * rew1 + som(3,ib,lp1) * rew2
            y2 =  sg(1,ib,lp1) + sg(2,ib,lp1) * rew1 + sg(3,ib,lp1) * rew2
          else
            x2 =  0.
            y2 =  0.
          end if
          dlt  = (eta(j,k) - eta_tab(l)) / (eta_tab(lp1) - eta_tab(l))
          omsm =  x1 + dlt * (x2 - x1)
          gsm  =  y1 + dlt * (y2 - y1)
          gsm2 =  gsm * gsm
          gsm3 =  gsm2 * gsm
          gsm4 =  gsm3 * gsm
        end if
        !
        tausw  =  wclw * (aws(1,ib) + aws(2,ib) / rew1 + &
                  aws(3,ib) / rew2 + aws(4,ib) / rew3)
        omsw   =  1.0 - (bws(1,ib) + bws(2,ib) * rew1 + &
                  bws(3,ib) * rew2 + bws(4,ib) * rew3) + omsm
        omsw   =  min(omsw, 0.999999)
        gsw    =  cws(1,ib,1) + cws(2,ib,1) * rew1 + &
                  cws(3,ib,1) * rew2 + cws(4,ib,1) * rew3 + gsm
        gsw2   =  cws(1,ib,2) + cws(2,ib,2) * rew1 + &
                  cws(3,ib,2) * rew2 + cws(4,ib,2) * rew3 + gsm2
        gsw3   =  cws(1,ib,3) + cws(2,ib,3) * rew1 + &
                  cws(3,ib,3) * rew2 + cws(4,ib,3) * rew3 + gsm3
        gsw4   =  cws(1,ib,4) + cws(2,ib,4) * rew1 + &
                  cws(3,ib,4) * rew2 + cws(4,ib,4) * rew3 + gsm4
        !
        tausi  =  wcli * ( ais(1,ib) + ais(2,ib) / dg + ais(3,ib) / dg2)
        omsi   =  1.0 - (bis(1,ib) + bis(2,ib) * dg + &
                  bis(3,ib) * dg2 + bis(4,ib) * dg3)
        omsi   =  min(omsi, 0.999999)
        gsi    =  cis(1,ib,1) + cis(2,ib,1) * dg + cis(3,ib,1) * dg2 + &
                  cis(4,ib,1) * dg3
        gsi2   =  cis(1,ib,2) + cis(2,ib,2) * dg + cis(3,ib,2) * dg2 + cis(4,ib,2) * dg3
        gsi3   =  cis(1,ib,3) + cis(2,ib,3) * dg + cis(3,ib,3) * dg2 + cis(4,ib,3) * dg3
        gsi4   =  cis(1,ib,4) + cis(2,ib,4) * dg + cis(3,ib,4) * dg2 + cis(4,ib,4) * dg3
        !
        taucsi =  tausw + tausi
        if (taucsi > 0.) then
          y1   =  omsw * tausw
          y2   =  omsi * tausi
          y12  =  y1 + y2
          omcs =  y12 / taucsi
          gcs  = (y1 * gsw + y2 * gsi) / y12
          gcs2 = (y1 * gsw2 + y2 * gsi2) / y12
          gcs3 = (y1 * gsw3 + y2 * gsi3) / y12
          gcs4 = (y1 * gsw4 + y2 * gsi4) / y12
        else
          omcs =  0.0
          gcs  =  0.0
          gcs2 =  0.0
          gcs3 =  0.0
          gcs4 =  0.0
        end if
      !
      !----------------------------------------------------------------------c
      !     scaling the cloud optical properties                             c
      !----------------------------------------------------------------------c
      !
        vi_taucsm     =  vi_taucs(i)
        vi_taucs(i)   =  vi_taucs(i) + taucsi
        if (mcica == 0) then
          rmugfac     =  (2.0 - rmug(i)) ** 0.40
          anufac      =  1.0 / (1.0 + 5.68 * anu(j,k) ** 1.4)
          x           =  taucsi + 9.2 * sqrt(vi_taucsm)
          taucsg(i,k) =  taucsi / (1.0 + 0.185 * x * rmugfac * anufac)
        else
          taucsg(i,k) =  taucsi
        end if
        if (cldt(j) > 0.) then
          x_cas       =  max(min(cldg(i,k)/cldt(j),1.),0.)
        else
          x_cas       = 0.0
        end if
        c20           =  taucsg(i,k) * omcs
        tauomc(i,k)   =  tauoma_g(i,k) + c20
        tauomgc(i,k)  =  tauomga_g(i,k) + c20 * gcs
        tauomgc_4str(i,k,1)  =  tauomga_g_4str(i,k,1) + c20 * gcs2
        tauomgc_4str(i,k,2)  =  tauomga_g_4str(i,k,2) + c20 * gcs3
        if (nstrm <= 2) then
           f2(i,k)    =  f1_g(i,k) + c20 * gcs2
        else !(nstrm == 4)
           f2(i,k)    =  f1_g(i,k) + c20 * gcs4
        end if
      end if
      !
      wrkag(i)        =  wrkag(i) + taua_g(i,k)
      wrkbg(i)        =  wrkbg(i) + x_cas * taucsi + taua_g(i,k)
    end do
  end do
  !

  return
end subroutine sw_cld_props4
!> \file
!> !! Computation of cloud solar optical properties. It contains three parts:
!! 1. Parameterization for solar optical properties based on effective
!! radius and cloud liquid (ice) water contents, the effective
!! variance is fixed, which shown in \cite Dobbie1999 for liquid clouds and
!! \cite Yang2013 for ice clouds.
!!\n
!!\n
!! 2. Black carbon semi-direct effect, the parameterization is based on ETA,
!! the details of which are detailed in \cite Li2011.
!!\n
!!\n
!! 3. We can choose McICA or 4 column cloud method for radiative
!! transfer, if the 4 column cloud method is used, the cloud optical
!! depth needs to be scaled down. For the 4 column cloud method, the
!! physics is shown in \cite Li2005.
