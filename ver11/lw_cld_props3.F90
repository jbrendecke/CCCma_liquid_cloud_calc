!> \file lw_cld_props3.F90
!>\brief Compute cloud optical properties at thermal wavelengths
!!
!! @author Jiangnan Li and Jason Cole
!
subroutine lw_cld_props3(tauci, omci, gci, psi, f2,      &   ! output
                         dz, cld, clw, cic, radeqvw, radeqvi, &   ! input
                         cut, il1, il2, ilg, lay, ib, nstrm)      ! control
  !
  !     * feb 10/2015 - j.li.     new version for gcm18+:
  !     *                         - revised ice cloud optical properties
  !     *                           from yang et al (2012).
  !     * feb 10/2009 - m.lazare. previous version lw_cld_props2 for gcm15h:
  !     *                         - calculate optical properties regardless
  !     *                           of cloud layer water content, it no
  !     *                           cutoff to zero below certain threshold.
  !     * oct 22/2007 - j.cole/   previous version lw_cld_props for gcm15g:
  !     *               m.lazare. compute cloud l/w optical properties
  !     *                         required for radiation. it is composed
  !     *                         of three parts:
  !     *                          1. assigns "cld" to "cldg".
  !     *                          2. calculation of {taucl,omcl,gcl}
  !     *                             migrated from clouds15 routines in
  !     *                             physics (separated into s/w and l/w
  !     *                             separate components).
  !     *                          3. calculation of resulting optical
  !     *                             properties arrays required for
  !     *                             radiation.
  !     *                         note that this means {taucl,omcl,gcl}
  !     *                         merely become scalars.
  !
  use datcldop4, only : awl_2str, bwl_2str, cwl, cwl_2str, awl_4str, bwl_4str, &
                        ail, bil, cil

  implicit none
  !
  !     * output data.
  !
  real, intent(out) , dimension(ilg,lay) :: tauci !< Cloud optical depth \f$[1]\f$
  real, intent(out) , dimension(ilg,lay) :: omci !< Cloud single scattering albedo \f$[1]\f$
  real, intent(out) , dimension(ilg,lay) :: gci !< Cloud asymmetry factor \f$[1]\f$
  real, intent(out) , dimension(ilg,lay,6) :: psi !< !<6 psi function for 4stream \f$[0]\f$
  real, intent(out) , dimension(ilg,lay) :: f2 !< Delta-Edington scaling factor, GCI\f$^2\f$ \f$[1]\f$
  !
  !     * input data.
  !
  real, intent(in) , dimension(ilg,lay) :: dz   !< Layer geometeric thickness \f$[m]\f$
  real, intent(in) , dimension(ilg,lay) :: cld !< Input cloud fraction \f$[1]\f
  real, intent(in) , dimension(ilg,lay) :: clw !< Cloud liquid water content \f$[gram/m^2]\f$
  real, intent(in) , dimension(ilg,lay) :: cic !< Cloud ice water content \f$[gram/m^2]\f$
  real, intent(in) , dimension(ilg,lay) :: radeqvw !< Liquid cloud effective radius \f$[um]\f$
  real, intent(in) , dimension(ilg,lay) :: radeqvi !< Ice cloud effective radius \f$[um]\f$

  real, intent(in) :: cut !< Minimum threshold for input cloud to be considered for radaitive transfer \f$[1]\f$

  integer, intent(in) :: il1   !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2   !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg   !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay   !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: ib    !< Longwave band number \f$[1]\f$
  integer, intent(in) :: nstrm !< Switch to 2 or 4 stream radiative transfer \f$[unitless]\f$
  !
  !=======================================================================
  !     Cloud longwave optical properties. It contains two parts:
  !     1. parameterization for longwave optical properties based on effective
  !     radius and cloud liquid (ice) water contents, the effective
  !     variance is fixed, which shown in Lindner and Li. (2000, J. CLIM.) and
  !     Yang et al (2012 JGR)
  !     2. parameter factor for longwave radiative perturbation method
  !      which is shown in Li (2000 Jas).
  !=======================================================================
  !
  !     * local data (scalar)
  !     uu3 = 3 * u * u = 1.1036383, u = 1 / sqrt(e)
  real, parameter :: uu3 = 1.1036383  !< Factor for perturbation method in longwave radiative transfer solution \f$[1]\f$
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
  real  :: taulw
  real  :: omlw
  real  :: glw, glw2, glw3, glw4
  real  :: tauli
  real  :: omli
  real  :: gli, gli2, gli3, gli4
  real  :: taucl
  real  :: omcl
  real  :: gcl, gcl2, gcl3, gcl4
  real  :: xinv, x1, x2, x3
  real  :: y1
  real  :: y2
  real  :: y12
  real :: absli
  integer :: i
  integer :: k
  !
  real, dimension(:, :), pointer :: awl !< Pointers to cloud water/ice LW optical properties data
  real, dimension(:, :), pointer :: bwl
  !
  !----------------------------------------------------------------------
  !     cloud radiative properties for radiation.
  !     taucs, omcs, gcs (taucl, omcl, gcl): optical depth, single
  !     scattering albedo, asymmetry factor for solar (infrared).
  !     radeqvw: effective radius(in micrometer) for water cloud
  !     radeqvi: effective radius(in micrometer) for ice cloud
  !     dg: geometry length for ice cloud
  !     wcdw (wcdi): liquid water (ice) content (in gram / m^3)
  !     wclw (wcli): liquid water (ice) path length (in gram / m^2)
  !     ccld: cloud fraction
  !     parameterization for water cloud:
  !     dobbie, etc. 1999, jgr, 104, 2067-2079
  !     lindner, t. h. and j. li., 2000, j. clim., 13, 1797-1805.
  !     parameterization for ice cloud:
  !     ice cloud optical property based on yang et al. (2012) jas
  !     and a serirs papers by group of u a &m and u wisconsin
  !     dg is the effective diameter

  !----------------------------------------------------------------------
  !
  if (nstrm == 0) then
     awl => awl_2str
     bwl => bwl_2str
  else if (nstrm == 4 .or. nstrm == 2) then
     awl => awl_4str
     bwl => bwl_4str
  end if

  psi(:, :, :) = 0.0
  gci(:, :) = 0.0
  do k = 1, lay
    do i = il1, il2
      if (cld(i,k) <= cut) then
        tauci(i,k) =  0.0
        omci(i,k)  =  0.0
        gci(i,k)   =  0.0
        f2(i,k)    =  0.0
      else
        wcdw       =  clw(i,k)
        wcdi       =  cic(i,k)
        wclw       =  wcdw * dz(i,k)
        wcli       =  wcdi * dz(i,k)
        if (nstrm == 0) then
           rew1    =  radeqvw(i,k)
        else
           rew1    =  min(max(radeqvw(i,k), 5.0), 20.0)
        end if
        rew2       =  rew1 * rew1
        rew3       =  rew2 * rew1
        dg         =  2.0 * min (max (radeqvi(i,k), 5.0), 60.0)
        dg2        =  dg  * dg
        dg3        =  dg2 * dg
        !
        taulw      =  wclw * (awl(1,ib) + awl(2,ib) * rew1 + awl(3,ib) / rew1 + &
                      awl(4,ib) / rew2  + awl(5,ib) / rew3)
        omlw       =  1.0 - (bwl(1,ib) + bwl(2,ib) / rew1 + bwl(3,ib) * rew1 + &
                      bwl(4,ib) * rew2)
        if (nstrm == 0) then
           glw     =  cwl_2str(1,ib) + cwl_2str(2,ib) / rew1 + cwl_2str(3,ib) * rew1 + &
                      cwl_2str(4,ib) * rew2
        else if (nstrm == 4 .or. nstrm == 2) then
           glw     =  cwl(1,ib,1) + cwl(2,ib,1) / rew1 + cwl(3,ib,1) * rew1 + cwl(4,ib,1) * rew2
        end if
        glw2       =  cwl(1,ib,2) + cwl(2,ib,2) / rew1 + cwl(3,ib,2) * rew1 + cwl(4,ib,2) * rew2
        glw3       =  cwl(1,ib,3) + cwl(2,ib,3) / rew1 + cwl(3,ib,3) * rew1 + cwl(4,ib,3) * rew2
        glw4       =  cwl(1,ib,4) + cwl(2,ib,4) / rew1 + cwl(3,ib,4) * rew1 + cwl(4,ib,4) * rew2
        !
        !----------------------------------------------------------------------
        !     since in yang, it is a param for absorptance depth
        !----------------------------------------------------------------------
        !
        absli      =  ail(1,ib) + ail(2,ib) / dg + ail(3,ib) / dg2
        omli       =  1.0 - (bil(1,ib) + bil(2,ib) * dg + bil(3,ib) * dg2 + &
                      bil(4,ib) * dg3)
        tauli      =  wcli * (absli / (1.0 - omli))
        gli        =  cil(1,ib,1) + cil(2,ib,1) * dg + cil(3,ib,1) * dg2 + cil(4,ib,1) * dg3
        gli2       =  cil(1,ib,2) + cil(2,ib,2) * dg + cil(3,ib,2) * dg2 + cil(4,ib,2) * dg3
        gli3       =  cil(1,ib,3) + cil(2,ib,3) * dg + cil(3,ib,3) * dg2 + cil(4,ib,3) * dg3
        gli4       =  cil(1,ib,4) + cil(2,ib,4) * dg + cil(3,ib,4) * dg2 + cil(4,ib,4) * dg3
        !
        taucl      =  taulw + tauli
        if (taucl > 0.) then
          y1       =  omlw * taulw
          y2       =  omli * tauli
          y12      =  y1 + y2
          omcl     =  y12 / taucl
          gcl      = (glw * y1 + gli * y2) / y12
          gcl2     = (glw2 * y1 + gli2 * y2) / y12
          gcl3     = (glw3 * y1 + gli3 * y2) / y12
          gcl4     = (glw4 * y1 + gli4 * y2) / y12
        else
          omcl     =  0.
          gcl      =  0.
          gcl2     =  0.
          gcl3     =  0.
          gcl4     =  0.
        end if
        !
        tauci(i,k) =  taucl
        omci(i,k)  =  omcl * tauci(i,k)
        if (nstrm <= 2) then
           if (nstrm == 0) then
              f2(i,k) =  gcl * gcl
           else
              f2(i,k) =  gcl2
           end if
           gci(i,k)   = (gcl - f2(i,k)) / (1.0 - f2(i,k))
           gci(i,k)   =  - 0.5 * (1.0 - uu3 * gci(i,k))
        else if (nstrm == 4) then
!
!...  aerosol scattering in longwave is not considered, uu3 see li (2002)
!     psi1 = psi(11), psi2=psi(12)=psi(21), psi3=psi(22), psi4 = psi(1-1),
!     psi5=psi(1-2)=psi(2-1), psi6=psi(2-2)
!
           f2(i,k)    =  gcl4
           xinv       =  1.0 / (1.0 - f2(i,k))
           x1         =  3.0 * (gcl - f2(i,k)) * xinv
           x2         =  5.0 * (gcl2 - f2(i,k)) * xinv
           x3         =  7.0 * (gcl3 - f2(i,k)) * xinv
!
!  in psi(i,k,1) and psi(i,k,3), factor of 2 because the weight = 0.5
!
           psi(i,k,1) =  2.0 - 0.5 * (1.0 + 0.0446582 * x1 + 0.1875000 * x2 + 0.0860800 * x3)
           psi(i,k,2) = -0.5 * (1.0 + 0.1666667 * x1 - 0.1875000 * x2 - 0.0127315 * x3)
           psi(i,k,3) =  2.0 - 0.5 * (1.0 + 0.6220084 * x1 + 0.1875000 * x2 + 0.0018830 * x3)
           psi(i,k,4) = -0.5 * (1.0 - 0.0446582 * x1 + 0.1875000 * x2 - 0.0860800 * x3)
           psi(i,k,5) = -0.5 * (1.0 - 0.1666667 * x1 - 0.1875000 * x2 + 0.0127315 * x3)
           psi(i,k,6) = -0.5 * (1.0 - 0.6220084 * x1 + 0.1875000 * x2 - 0.0018830 * x3)
        end if
      end if
    end do
  end do

  return
end subroutine lw_cld_props3
