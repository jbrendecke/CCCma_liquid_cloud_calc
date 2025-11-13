!> \file
!>\brief 4-stream spherical harmonic expansion for sw radiative transfer.
!!
!! @author Jiangnan Li
!
 subroutine swtran4st(reflt, trant, cumdtr, attens, taua, taur, taug,  &
                      tauoma, tauomga, tauomga_4str, f1, f2, tauc,     &
                      tauomc, tauomgc, tauomgc_4str, cld, rmu, albsur, &
                      csalb, nct, cut, il1, il2, ilg, lay, lev, ntile)
!
!----------------------------------------------------------------------
!     4-stream spherical harmonic expansion for sw radiative transfer
!     doubling/adding based on zhang & li (jas 2013), single layer
!     solution based on li & ramaswamy (jas 1996)
!----------------------------------------------------------------------
!
!     refl:      reflectivity (1) clear sky; (2) all sky
!     tran:      transmitivity
!     cumdtr:    direct transmission for mult-layers (no scaling)
!     scumdtr:   direct transmission for mult-layers (scaling)
!     attens:     direct flux transmission at level 1
!     taua:      aerosol optical depth
!     taur:      rayleigh optical depth
!     taug:      gaseous optical depth
!     tauoma:    aerosol optical depth times aerosol single scattering albedo
!     tauomga:   tauoma times aerosol asymmetry factor
!     f1:        square of aerosol asymmetry factor
!     f2:        square of cloud asymmetry factor
!     tauc:      cloud optical depth
!     tauomc:    cloud optical depth times cloud single
!                scattering albedo
!     tauomgc:   tauomc times cloud asymmetry factor
!     cld:       cloud fraction (assumed to be 0 or 1)
!     nct:       the highest cloud top level, if nct=lev, no cloud in this column
!     rmu:       cos of solar zenith angle
!     albsur:    surface albedo
!     csalb:     clear sky surface albedo
!----------------------------------------------------------------------
!
  implicit none
!
  integer, intent(in) :: il1, il2, ilg, lay, lev
  integer, intent(in) :: ntile !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  real,    intent(in) :: cut   !<cloud fraction criterion \f$[0]\f$
!
  real, intent(out), dimension(ilg,2,lev) :: cumdtr !< Direct transmission for multi-layers using unscale cloud properties \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: reflt !< Tiled layer reflectivity for all sky and clear sky \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,2,lev) :: trant !< Tiled layer transmission for all sky and clear sky \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,lay)   :: taua    !< Aerosol optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: taur    !< Rayleigh optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: taug    !< Gaseous optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauc    !< Cloud optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauoma  !< taua times aerosol single scattering albedo\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauomc  !< tauc times cloud single scattering albedo\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauomga !< tauoma times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay,2) :: tauomga_4str !< tauoma times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauomgc !< tauomc times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay,2) :: tauomgc_4str !< tauomc times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: f1      !< tauomga times the square of first momemt\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: f2      !< tauomgc times the square of first momemt\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: cld     !< cloud fraction\f$[0]\f$
  real, intent(in), dimension(ilg)       :: attens  !< direct flux transmission at level 1 \f$[Wm^(-2)]\f$
  real, intent(in), dimension(ilg)       :: rmu     !< 1/(cosine of solar zenith angle\f$[0]\f$
  real, intent(in), dimension(ilg,ntile) :: albsur  !< all sky surface albedo \f$[0]\f$
  real, intent(in), dimension(ilg,ntile) :: csalb   !< clear sky surface albedo \f$[0]\f$
  integer, intent(in), dimension(ilg)    :: nct     !< the highest cloud top level in the column \f$[0]\f$
!
  real*8, dimension(ilg,2,lev)             :: scumdtr !<direct transmission (scaling)\f$[Wm^(-2)]\f$
  real*8, dimension(ilg,2,lev)             :: dtr     !<layer direct transmission 1 clear sky, 2 all sky \f$[0]\f$
!
  real*8, dimension(ilg,lay,2,2)           :: rdfc, rdfa, tdfc, tdfa
  real*8, dimension(ilg,lay,2)             :: rdrc, rdra, tdrc, tdra
  real*8, dimension(ilg,ntile,lev,2,2)     :: rmdfc, rmdfa, rmufc, rmufa
  real*8, dimension(ilg,ntile,lev,2)       :: rmurc, rmura, tmdrc, tmdra
!
  integer :: i,  m,  k,  km1,  l,  lp1
  real    :: extopt, extopta
  real*8  :: c11, c12, c21, c22, pmod, t11, t12, t21, t22, b11, b12, &
             b21, b22, d1, d2, uu1, uu2, du1, du2, a11, a12, a21, a22, r11, r12, r21, r22
!
  call swtrlayer(rdfa, tdfa, rdra, tdra, dtr, rdfc, tdfc, rdrc, tdrc, &
                 taua, taur, taug, tauoma, tauomga, tauomga_4str, f1, &
                 f2, tauc, tauomc, tauomgc, tauomgc_4str,             &
                 rmu, cld, cut, il1, il2, ilg, lay)
!
  cumdtr(il1:il2,1:2,1:lev) = 0.0
!
  do i = il1, il2
    do m = 1, ntile
      rmdfa(i,m,1,1:2,1)    =  0.0d0
      rmdfa(i,m,1,1:2,2)    =  0.0d0
      tmdra(i,m,1,1:2)      =  0.0d0
      rmdfc(i,m,1,1:2,1)    =  0.0d0
      rmdfc(i,m,1,1:2,2)    =  0.0d0
      tmdrc(i,m,1,1:2)      =  0.0d0
      cumdtr(i,1:2,1)       =  attens(i)
      scumdtr(i,1:2,1)      =  dble(attens(i))
!
! ... nitializaton for surface. a:clear sky; c:cloudy sky.
!     rmura(i,m,lev,2) = 0. since intgral p3(mu) dmu = 0.
!
      rmura(i,m,lev,1)      =  dble(csalb(i,m))
      rmura(i,m,lev,2)      = -0.25d0 * dble(csalb(i,m))
      rmufa(i,m,lev,1:2,1)  =  rmura(i,m,lev,1:2)
      rmufa(i,m,lev,1:2,2)  =  0.0d0
!
      rmurc(i,m,lev,1)      =  dble(albsur(i,m))
      rmurc(i,m,lev,2)      = -0.25d0 * dble(albsur(i,m))
      rmufc(i,m,lev,1:2,1)  =  rmurc(i,m,lev,1:2)
      rmufc(i,m,lev,1:2,2)  =  0.0d0
    enddo
  enddo
!
! ... add the layers upward from one layer above surface to the level 1.
!
  do k = 2, lev
    km1 = k - 1
    do i = il1, il2
      do m = 1, ntile
!
! ...  for clear sky
!
        extopta          =  taua(i,km1) + taur(i,km1) + taug(i,km1)
        cumdtr(i,1,k)    =  cumdtr(i,1,km1) * exp(- extopta / rmu(i))
        c12              =  - rdfa(i,km1,1,1) * rmdfa(i,m,km1,1,2) - rdfa(i,km1,1,2) * rmdfa(i,m,km1,2,2)
        c11              =  1.0d0 + c12
        c21              =  - rdfa(i,km1,2,1) * rmdfa(i,m,km1,1,1) - rdfa(i,km1,2,2) * rmdfa(i,m,km1,2,1)
        c22              =  1.0d0 + c21
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod

        b11              =  t11 * tdfa(i,km1,1,1) + t12 * tdfa(i,km1,2,1)
        b12              =  t11 * tdfa(i,km1,1,2) + t12 * tdfa(i,km1,2,2)
        b21              =  t21 * tdfa(i,km1,1,1) + t22 * tdfa(i,km1,2,1)
        b22              =  t21 * tdfa(i,km1,1,2) + t22 * tdfa(i,km1,2,2)
!
        scumdtr(i,1,k)   =  scumdtr(i,1,km1) * dtr(i,1,km1)
        d1               =  rdra(i,km1,1) * scumdtr(i,1,km1) + rdfa(i,km1,1,1) * tmdra(i,m,km1,1) + &
                            rdfa(i,km1,1,2) * tmdra(i,m,km1,2)
        d2               =  rdra(i,km1,2) * scumdtr(i,1,km1) + rdfa(i,km1,2,1) * tmdra(i,m,km1,1) + &
                            rdfa(i,km1,2,2) * tmdra(i,m,km1,2)
        uu1              =  d1 * t11 + d2 * t12
        uu2              =  d1 * t21 + d2 * t22
        du1              =  tmdra(i,m,km1,1) + uu1 * rmdfa(i,m,km1,1,1) + uu2 * rmdfa(i,m,km1,1,2)
        du2              =  tmdra(i,m,km1,2) + uu1 * rmdfa(i,m,km1,2,1) + uu2 * rmdfa(i,m,km1,2,2)
        tmdra(i,m,k,1)   =  tdra(i,km1,1) * scumdtr(i,1,km1) + du1 * tdfa(i,km1,1,1) + du2 * tdfa(i,km1,1,2)
        tmdra(i,m,k,2)   =  tdra(i,km1,2) * scumdtr(i,1,km1) + du1 * tdfa(i,km1,2,1) + du2 * tdfa(i,km1,2,2)
!
        a11              =  tdfa(i,km1,1,1) * rmdfa(i,m,km1,1,1) + tdfa(i,km1,1,2) * rmdfa(i,m,km1,2,1)
        a12              =  tdfa(i,km1,1,1) * rmdfa(i,m,km1,1,2) + tdfa(i,km1,1,2) * rmdfa(i,m,km1,2,2)
        a21              =  tdfa(i,km1,2,1) * rmdfa(i,m,km1,1,1) + tdfa(i,km1,2,2) * rmdfa(i,m,km1,2,1)
        a22              =  tdfa(i,km1,2,1) * rmdfa(i,m,km1,1,2) + tdfa(i,km1,2,2) * rmdfa(i,m,km1,2,2)
        rmdfa(i,m,k,1,1) =  rdfa(i,km1,1,1) + a11 * b11 + a12 * b21
        rmdfa(i,m,k,1,2) =  rdfa(i,km1,1,2) + a11 * b12 + a12 * b22
        rmdfa(i,m,k,2,1) =  rdfa(i,km1,2,1) + a21 * b11 + a22 * b21
        rmdfa(i,m,k,2,2) =  rdfa(i,km1,2,2) + a21 * b12 + a22 * b22
!
! ...  for all sky
!
        if (k < nct(i))  then
          cumdtr(i,2,k)      =  cumdtr(i,1,k)
          scumdtr(i,2,k)     =  scumdtr(i,1,k)
          tmdrc(i,m,k,1:2)   =  tmdra(i,m,k,1:2)
          rmdfc(i,m,k,1:2,1) =  rmdfa(i,m,k,1:2,1)
          rmdfc(i,m,k,1:2,2) =  rmdfa(i,m,k,1:2,2)
        else
          extopt             =  extopta + tauc(i,km1)
          cumdtr(i,2,k)      =  cumdtr(i,2,km1) * exp(- extopt / rmu(i))
          c11                =  1.0d0 - rdfc(i,km1,1,1) * rmdfc(i,m,km1,1,1) - rdfc(i,km1,1,2) * rmdfc(i,m,km1,2,1)
          c12                =      - rdfc(i,km1,1,1) * rmdfc(i,m,km1,1,2) - rdfc(i,km1,1,2) * rmdfc(i,m,km1,2,2)
          c21                =      - rdfc(i,km1,2,1) * rmdfc(i,m,km1,1,1) - rdfc(i,km1,2,2) * rmdfc(i,m,km1,2,1)
          c22                =  1.0d0 - rdfc(i,km1,2,1) * rmdfc(i,m,km1,1,2) - rdfc(i,km1,2,2) * rmdfc(i,m,km1,2,2)
          pmod               =  c11 * c22 - c12 * c21
          t11                =  c22 / pmod
          t12                = -c12 / pmod
          t21                = -c21 / pmod
          t22                =  c11 / pmod

          b11                =  t11 * tdfc(i,km1,1,1) + t12 * tdfc(i,km1,2,1)
          b12                =  t11 * tdfc(i,km1,1,2) + t12 * tdfc(i,km1,2,2)
          b21                =  t21 * tdfc(i,km1,1,1) + t22 * tdfc(i,km1,2,1)
          b22                =  t21 * tdfc(i,km1,1,2) + t22 * tdfc(i,km1,2,2)
!
          scumdtr(i,2,k)     =  scumdtr(i,2,km1) * dtr(i,2,km1)
          d1                 =  rdrc(i,km1,1) * scumdtr(i,2,km1) + rdfc(i,km1,1,1) * tmdrc(i,m,km1,1) + &
                                rdfc(i,km1,1,2) * tmdrc(i,m,km1,2)
          d2                 =  rdrc(i,km1,2) * scumdtr(i,2,km1) + rdfc(i,km1,2,1) * tmdrc(i,m,km1,1) + &
                                rdfc(i,km1,2,2) * tmdrc(i,m,km1,2)
          uu1                =  d1 * t11 + d2 * t12
          uu2                =  d1 * t21 + d2 * t22
          du1                =  tmdrc(i,m,km1,1) + uu1 * rmdfc(i,m,km1,1,1) + uu2 * rmdfc(i,m,km1,1,2)
          du2                =  tmdrc(i,m,km1,2) + uu1 * rmdfc(i,m,km1,2,1) + uu2 * rmdfc(i,m,km1,2,2)
          tmdrc(i,m,k,1)     =  tdrc(i,km1,1) * scumdtr(i,2,km1) + du1 * tdfc(i,km1,1,1) + du2 * tdfc(i,km1,1,2)
          tmdrc(i,m,k,2)     =  tdrc(i,km1,2) * scumdtr(i,2,km1) + du1 * tdfc(i,km1,2,1) + du2 * tdfc(i,km1,2,2)
!
          a11                =  tdfc(i,km1,1,1) * rmdfc(i,m,km1,1,1) + tdfc(i,km1,1,2) * rmdfc(i,m,km1,2,1)
          a12                =  tdfc(i,km1,1,1) * rmdfc(i,m,km1,1,2) + tdfc(i,km1,1,2) * rmdfc(i,m,km1,2,2)
          a21                =  tdfc(i,km1,2,1) * rmdfc(i,m,km1,1,1) + tdfc(i,km1,2,2) * rmdfc(i,m,km1,2,1)
          a22                =  tdfc(i,km1,2,1) * rmdfc(i,m,km1,1,2) + tdfc(i,km1,2,2) * rmdfc(i,m,km1,2,2)
          rmdfc(i,m,k,1,1)   =  rdfc(i,km1,1,1) + a11 * b11 + a12 * b21
          rmdfc(i,m,k,1,2)   =  rdfc(i,km1,1,2) + a11 * b12 + a12 * b22
          rmdfc(i,m,k,2,1)   =  rdfc(i,km1,2,1) + a21 * b11 + a22 * b21
          rmdfc(i,m,k,2,2)   =  rdfc(i,km1,2,2) + a21 * b12 + a22 * b22
        endif
      enddo
    enddo
  enddo
!
  do l = lay, 1, -1
    lp1 = l + 1
    do i = il1, il2
      do m = 1, ntile
!
! ... add the layers upward from layer above surface to the level 1.
!
        c11              =  1.0d0 - rmufa(i,m,lp1,1,1) * rdfa(i,l,1,1) - rmufa(i,m,lp1,1,2) * rdfa(i,l,2,1)
        c12              =      - rmufa(i,m,lp1,1,1) * rdfa(i,l,1,2) - rmufa(i,m,lp1,1,2) * rdfa(i,l,2,2)
        c21              =      - rmufa(i,m,lp1,2,1) * rdfa(i,l,1,1) - rmufa(i,m,lp1,2,2) * rdfa(i,l,2,1)
        c22              =  1.0d0 - rmufa(i,m,lp1,2,1) * rdfa(i,l,1,2) - rmufa(i,m,lp1,2,2) * rdfa(i,l,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
        d1               =  rmura(i,m,lp1,1) * dtr(i,1,l) + rmufa(i,m,lp1,1,1) * tdra(i,l,1) + &
                            rmufa(i,m,lp1,1,2) * tdra(i,l,2)
        d2               =  rmura(i,m,lp1,2) * dtr(i,1,l) + rmufa(i,m,lp1,2,1) * tdra(i,l,1) + &
                            rmufa(i,m,lp1,2,2) * tdra(i,l,2)
        uu1              =  d1 * t11 + d2 * t12
        uu2              =  d1 * t21 + d2 * t22
        rmura(i,m,l,1)   =  rdra(i,l,1) + uu1 * tdfa(i,l,1,1) + uu2 * tdfa(i,l,1,2)
        rmura(i,m,l,2)   =  rdra(i,l,2) + uu1 * tdfa(i,l,2,1) + uu2 * tdfa(i,l,2,2)
!
        c11              =  1.0d0 - rdfa(i,l,1,1) * rmufa(i,m,lp1,1,1) - rdfa(i,l,1,2) * rmufa(i,m,lp1,2,1)
        c12              =      - rdfa(i,l,1,1) * rmufa(i,m,lp1,1,2) - rdfa(i,l,1,2) * rmufa(i,m,lp1,2,2)
        c21              =      - rdfa(i,l,2,1) * rmufa(i,m,lp1,1,1) - rdfa(i,l,2,2) * rmufa(i,m,lp1,2,1)
        c22              =  1.0d0 - rdfa(i,l,2,1) * rmufa(i,m,lp1,1,2) - rdfa(i,l,2,2) * rmufa(i,m,lp1,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
!
        r11              =  tdfa(i,l,1,1) * (t11 * rmufa(i,m,lp1,1,1) +  t21 * rmufa(i,m,lp1,1,2)) + &
                            tdfa(i,l,1,2) * (t11 * rmufa(i,m,lp1,2,1) +  t21 * rmufa(i,m,lp1,2,2))
        r12              =  tdfa(i,l,1,1) * (t12 * rmufa(i,m,lp1,1,1) +  t22 * rmufa(i,m,lp1,1,2)) + &
                            tdfa(i,l,1,2) * (t12 * rmufa(i,m,lp1,2,1) +  t22 * rmufa(i,m,lp1,2,2))
        r21              =  tdfa(i,l,2,1) * (t11 * rmufa(i,m,lp1,1,1) +  t21 * rmufa(i,m,lp1,1,2)) + &
                            tdfa(i,l,2,2) * (t11 * rmufa(i,m,lp1,2,1) +  t21 * rmufa(i,m,lp1,2,2))
        r22              =  tdfa(i,l,2,1) * (t12 * rmufa(i,m,lp1,1,1) +  t22 * rmufa(i,m,lp1,1,2)) + &
                            tdfa(i,l,2,2) * (t12 * rmufa(i,m,lp1,2,1) +  t22 * rmufa(i,m,lp1,2,2))
!
        rmufa(i,m,l,1,1) =  rdfa(i,l,1,1) + r11 * tdfa(i,l,1,1) + r12 * tdfa(i,l,2,1)
        rmufa(i,m,l,1,2) =  rdfa(i,l,1,2) + r11 * tdfa(i,l,1,2) + r12 * tdfa(i,l,2,2)
        rmufa(i,m,l,2,1) =  rdfa(i,l,2,1) + r21 * tdfa(i,l,1,1) + r22 * tdfa(i,l,2,1)
        rmufa(i,m,l,2,2) =  rdfa(i,l,2,2) + r21 * tdfa(i,l,1,2) + r22 * tdfa(i,l,2,2)
!
! ... add the layers upward from layer above surface to the level 1. since
!     the surface albedo is diff for cloudy and clear sky, the clear sky
!     result can not be used for the upward path.
!
        if (nct(i) == lev)  then
          rmurc(i,m,l,1:2)   =  rmura(i,m,l,1:2)
          rmufc(i,m,l,1,1:2) =  rmufa(i,m,l,1,1:2)
          rmufc(i,m,l,2,1:2) =  rmufa(i,m,l,2,1:2)
        else
          c11                =  1.0d0 - rmufc(i,m,lp1,1,1) * rdfc(i,l,1,1) - rmufc(i,m,lp1,1,2) * rdfc(i,l,2,1)
          c12                =      - rmufc(i,m,lp1,1,1) * rdfc(i,l,1,2) - rmufc(i,m,lp1,1,2) * rdfc(i,l,2,2)
          c21                =      - rmufc(i,m,lp1,2,1) * rdfc(i,l,1,1) - rmufc(i,m,lp1,2,2) * rdfc(i,l,2,1)
          c22                =  1.0d0 - rmufc(i,m,lp1,2,1) * rdfc(i,l,1,2) - rmufc(i,m,lp1,2,2) * rdfc(i,l,2,2)
          pmod               =  c11 * c22 - c12 * c21
          t11                =  c22 / pmod
          t12                = -c12 / pmod
          t21                = -c21 / pmod
          t22                =  c11 / pmod
          d1                 =  rmurc(i,m,lp1,1) * dtr(i,2,l) + rmufc(i,m,lp1,1,1) * tdrc(i,l,1) + &
                                rmufc(i,m,lp1,1,2) * tdrc(i,l,2)
          d2                 =  rmurc(i,m,lp1,2) * dtr(i,2,l) + rmufc(i,m,lp1,2,1) * tdrc(i,l,1) + &
                                rmufc(i,m,lp1,2,2) * tdrc(i,l,2)
          uu1                =  d1 * t11 + d2 * t12
          uu2                =  d1 * t21 + d2 * t22
          rmurc(i,m,l,1)     =  rdrc(i,l,1) + uu1 * tdfc(i,l,1,1) + uu2 * tdfc(i,l,1,2)
          rmurc(i,m,l,2)     =  rdrc(i,l,2) + uu1 * tdfc(i,l,2,1) + uu2 * tdfc(i,l,2,2)
!
          c11                =  1.0d0 - rdfc(i,l,1,1) * rmufc(i,m,lp1,1,1) - rdfc(i,l,1,2) * rmufc(i,m,lp1,2,1)
          c12                =      - rdfc(i,l,1,1) * rmufc(i,m,lp1,1,2) - rdfc(i,l,1,2) * rmufc(i,m,lp1,2,2)
          c21                =      - rdfc(i,l,2,1) * rmufc(i,m,lp1,1,1) - rdfc(i,l,2,2) * rmufc(i,m,lp1,2,1)
          c22                =  1.0d0 - rdfc(i,l,2,1) * rmufc(i,m,lp1,1,2) - rdfc(i,l,2,2) * rmufc(i,m,lp1,2,2)
          pmod               =  c11 * c22 - c12 * c21
          t11                =  c22 / pmod
          t12                = -c12 / pmod
          t21                = -c21 / pmod
          t22                =  c11 / pmod
!
          r11                =  tdfc(i,l,1,1) * (t11 * rmufc(i,m,lp1,1,1) + t21 * rmufc(i,m,lp1,1,2)) + &
                                tdfc(i,l,1,2) * (t11 * rmufc(i,m,lp1,2,1) + t21 * rmufc(i,m,lp1,2,2))
          r12                =  tdfc(i,l,1,1) * (t12 * rmufc(i,m,lp1,1,1) + t22 * rmufc(i,m,lp1,1,2)) + &
                                tdfc(i,l,1,2) * (t12 * rmufc(i,m,lp1,2,1) + t22 * rmufc(i,m,lp1,2,2))
          r21                =  tdfc(i,l,2,1) * (t11 * rmufc(i,m,lp1,1,1) + t21 * rmufc(i,m,lp1,1,2)) + &
                                tdfc(i,l,2,2) * (t11 * rmufc(i,m,lp1,2,1) + t21 * rmufc(i,m,lp1,2,2))
          r22                =  tdfc(i,l,2,1) * (t12 * rmufc(i,m,lp1,1,1) + t22 * rmufc(i,m,lp1,1,2)) + &
                                tdfc(i,l,2,2) * (t12 * rmufc(i,m,lp1,2,1) + t22 * rmufc(i,m,lp1,2,2))
!
          rmufc(i,m,l,1,1)   =  rdfc(i,l,1,1) + r11 * tdfc(i,l,1,1) + r12 * tdfc(i,l,2,1)
          rmufc(i,m,l,1,2)   =  rdfc(i,l,1,2) + r11 * tdfc(i,l,1,2) + r12 * tdfc(i,l,2,2)
          rmufc(i,m,l,2,1)   =  rdfc(i,l,2,1) + r21 * tdfc(i,l,1,1) + r22 * tdfc(i,l,2,1)
          rmufc(i,m,l,2,2)   =  rdfc(i,l,2,2) + r21 * tdfc(i,l,1,2) + r22 * tdfc(i,l,2,2)
        endif
      enddo
    enddo
  enddo
!
!----------------------------------------------------------------------
!     add downward to calculate the resultant reflectances and
!     transmittance at flux levels.
!----------------------------------------------------------------------
!
  do k = 1, lev
    do i = il1, il2
      do m = 1, ntile
        c11              =  1.0d0 - rmufa(i,m,k,1,1) * rmdfa(i,m,k,1,1) - rmufa(i,m,k,1,2) * rmdfa(i,m,k,2,1)
        c12              =      - rmufa(i,m,k,1,1) * rmdfa(i,m,k,1,2) - rmufa(i,m,k,1,2) * rmdfa(i,m,k,2,2)
        c21              =      - rmufa(i,m,k,2,1) * rmdfa(i,m,k,1,1) - rmufa(i,m,k,2,2) * rmdfa(i,m,k,2,1)
        c22              =  1.0d0 - rmufa(i,m,k,2,1) * rmdfa(i,m,k,1,2) - rmufa(i,m,k,2,2) * rmdfa(i,m,k,2,2)
        pmod             =  c11 * c22 - c12 * c21
        t11              =  c22 / pmod
        t12              = -c12 / pmod
        t21              = -c21 / pmod
        t22              =  c11 / pmod
!
        d1               =  rmura(i,m,k,1) * scumdtr(i,1,k) + rmufa(i,m,k,1,1) * tmdra(i,m,k,1) + &
                            rmufa(i,m,k,1,2) * tmdra(i,m,k,2)
        d2               =  rmura(i,m,k,2) * scumdtr(i,1,k) + rmufa(i,m,k,2,1) * tmdra(i,m,k,1) + &
                            rmufa(i,m,k,2,2) * tmdra(i,m,k,2)
        uu1              =  d1 * t11 + d2 * t12
        du1              =  tmdra(i,m,k,1) + rmdfa(i,m,k,1,1) * uu1 + rmdfa(i,m,k,1,2) * (d1 * t21 + d2 * t22)
        reflt(i,m,1,k)   =  real(uu1)
        trant(i,m,1,k)   =  real((du1 + scumdtr(i,1,k)))
!
! ...  for all sky
!
        if(nct(i) == lev)  then
          reflt(i,m,2,k) =  reflt(i,m,1,k)
          trant(i,m,2,k) =  trant(i,m,1,k)
        else
          c11            =  1.0d0 - rmufc(i,m,k,1,1) * rmdfc(i,m,k,1,1) - rmufc(i,m,k,1,2) * rmdfc(i,m,k,2,1)
          c12            =      - rmufc(i,m,k,1,1) * rmdfc(i,m,k,1,2) - rmufc(i,m,k,1,2) * rmdfc(i,m,k,2,2)
          c21            =      - rmufc(i,m,k,2,1) * rmdfc(i,m,k,1,1) - rmufc(i,m,k,2,2) * rmdfc(i,m,k,2,1)
          c22            =  1.0d0 - rmufc(i,m,k,2,1) * rmdfc(i,m,k,1,2) - rmufc(i,m,k,2,2) * rmdfc(i,m,k,2,2)
          pmod           =  c11 * c22 - c12 * c21
          t11            =  c22 / pmod
          t12            = -c12 / pmod
          t21            = -c21 / pmod
          t22            =  c11 / pmod
!
          d1             =  rmurc(i,m,k,1) * scumdtr(i,2,k) + rmufc(i,m,k,1,1) * tmdrc(i,m,k,1) + &
                            rmufc(i,m,k,1,2) * tmdrc(i,m,k,2)
          d2             =  rmurc(i,m,k,2) * scumdtr(i,2,k) + rmufc(i,m,k,2,1) * tmdrc(i,m,k,1) + &
                            rmufc(i,m,k,2,2) * tmdrc(i,m,k,2)
          uu1            =  d1 * t11 + d2 * t12
          du1            =  tmdrc(i,m,k,1) + rmdfc(i,m,k,1,1) * uu1 + rmdfc(i,m,k,1,2) * (d1 * t21 + d2 * t22)
          reflt(i,m,2,k) =  real(uu1)
          trant(i,m,2,k) = real((du1 + scumdtr(i,2,k)))
        endif
      enddo
    enddo
  enddo
!
  return

  contains

  subroutine swtrlayer(rdfa, tdfa, rdra, tdra, dtr, rdfc, tdfc, rdrc, tdrc, &
                       taua, taur, taug, tauoma, tauomga, tauomga_4str,     &
                       f1, f2, tauc, tauomc, tauomgc, tauomgc_4str, rmu,    &
                       cld, cut, il1, il2, ilg, lay)
!
!----------------------------------------------------------------------
!  single layer solution for 4-strream spherical harmonic expansion
!  li   & ramaswamy jas 1996)
!----------------------------------------------------------------------
!
  implicit none
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  real, intent(in)    :: cut  !< cloud fraction criterion \f$[0]\f$
  real*8, intent(out), dimension(ilg,lay,2,2) :: rdfc    !< clear sky layer diffuse reflection\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2,2) :: tdfc    !< clear sky layer diffuse transmission\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2,2) :: rdfa    !< all sky layer diffuse reflection\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2,2) :: tdfa    !< all sky layer diffuse transmission\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2)   :: rdrc    !< clear sky layer direct reflection\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2)   :: tdrc    !< clear sky layer direct transmission\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2)   :: rdra    !< all sky layer direct reflection\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,lay,2)   :: tdra    !< all sky layer direct transmission\f$[w/m^2]\f$
  real*8, intent(out), dimension(ilg,2,lay)   :: dtr     !<direct transmission\f$[0]\f$
!
  real, intent(in), dimension(ilg,lay) :: taua      !< Aerosol optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: taur      !< Rayleigh optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: taug      !< Gaseous optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: tauc      !< Cloud optical thickness \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: tauoma    !<taua times aerosol single scattering albedo\f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: tauomc    !<tauc times cloud single scattering albedo\f$[0]\f$

  real, intent(in), dimension(ilg,lay)   :: tauomga !<tauoma times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay,2) :: tauomga_4str !<tauoma times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: tauomgc !<tauomc times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay,2) :: tauomgc_4str !<tauomc times three momemts\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: f1      !<tauomga times the square of first momemt\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: f2      !<tauomgc times the square of first momemt\f$[0]\f$
  real, intent(in), dimension(ilg,lay)   :: cld     !<cloud fraction\f$[0]\f$
  real, intent(in), dimension(ilg)       :: rmu     !<1/(cosine of solar zenith angle\f$[0]\f$
!
  real*8 :: dmu,rmu2,rmu3,dmu2,dmu4,x1,x2,extopta,omars,ssalb,g,g2,g3,f,ssalbf,om,tau, &
            omega1,omega2,omega3,a0,a1,a2,a3,a01,a23,a03,b0,b1,b2,b3,a3b2,a0b1,dmub0,dmub3, &
            uu,uu2,vv,alams1,alams2,alam1,alam2,tr1,tr2,delta,eta0,eta1,eta2,eta3,p1,q1, &
            r1,p2,q2,r2,h1,h2,h3,h4,w11,w12,w13,w14,w21,w22,w23,w24,wa,wb,wc,wd,we,wf,det, &
            ya,yb,yc,yd,ye,yf,yg,yh,c1,d1,c2,d2,c1t,c2t,d1t,d2t,ted0,ted1,ted2,ted3, &
            teu0,teu1,teu2,teu3,extopt,omarcs
  integer :: i, k
  real*8, parameter :: bound = 1.0d-20
  real*8, parameter :: d9 = 1.0d0 / 9.0d0
  real*8, parameter :: eps_10 = 1.0d-10
!
  do k = 1, lay
    do i = il1, il2
      dmu             =  1.0d0 / dble(rmu(i))
      rmu2            =  dble(rmu(i) * rmu(i))
      rmu3            =  rmu2 * dble(rmu(i))
      dmu2            =  dmu * dmu
      dmu4            =  dmu2 * dmu2
!
! ...   clear sky calculation
!
      x1              =  dble(taua(i,k) + taur(i,k) + taug(i,k))
      x2              =  dble(tauoma(i,k) + taur(i,k))
      extopta         =  max(x1, bound)
      omars           =  max(x2, bound)
      ssalb           =  omars / extopta
      g               =  dble(tauomga(i,k)) / omars
      g2              =  dble(tauomga_4str(i,k,1)) / omars
      g3              =  dble(tauomga_4str(i,k,2)) / omars
      f               =  dble(f1(i,k)) / omars
!
      ssalbf          =  ssalb * f
      om              =  ssalb * (1.0d0 - f) / (1.0d0 - ssalbf)
      tau             =  extopta * (1.0d0 - ssalbf)
      omega1          =  3.0d0 * om * (g - f) / (1.0d0 - f)
      omega2          =  5.0d0 * om * (g2 - f) / (1.0d0 - f)
      omega3          =  7.0d0 * om * (g3 - f) / (1.0d0 - f)
!
      a0              =  1.0d0 - om + eps_10
      a1              =  3.0d0 - omega1 + eps_10
      a2              =  5.0d0 - omega2 + eps_10
      a3              =  7.0d0 - omega3 + eps_10
      a01             =  a0 * a1
      a23             =  a2 * a3
      a03             =  a0 * a3
!
      b0              =  0.25d0 * om
      b1              = -0.25d0 * omega1 * dble(rmu(i))
      b2              =  0.125d0 * omega2 * (3.0d0 * rmu2 - 1.0d0)
      b3              = -0.125d0 * omega3 * (5.0d0 * rmu3 - 3.0d0 * dble(rmu(i)))
      a3b2            =  a3 * b2
      a0b1            =  a0 * b1
      dmub0           =  dmu * b0
      dmub3           =  dmu * b3
!
      uu              =  a01 + d9 * (a23 + 4.0d0 * a03)
      uu2             =  uu * uu
      vv              =  4.0d0 * d9 * a01 * a23
!
      alams1          =  0.5d0 * (uu + sqrt(uu2 - vv))
      alams2          =  0.5d0 * (uu - sqrt(uu2 - vv))
      alam1           =  sqrt(alams1)
      alam2           =  sqrt(alams2)
!
      dtr(i,1,k)      =  exp(- tau * dmu)
      tr1             =  exp(- alam1 * tau)
      tr2             =  exp(- alam2 * tau)
!
      delta           =  1.0d0 / (9.0d0 * dmu4 - dmu2 * (9.0d0 * a01 + a23 + 4.0d0 * a03) + a01 * a23)
      eta0            = ((a1 * b0 - dmu * b1) * (a23 - 9.0d0 * dmu2) + &
                          2.0d0 * dmu2 * (a3b2 - 2.0d0 * a3 * b0 - 3.0d0 * dmub3)) * delta
      eta1            = ((a0b1 - dmub0) * (a23 - 9.0d0 * dmu2) - &
                          2.0d0 * dmu * a0 * (a3b2 - 3.0d0 * dmub3)) * delta
      eta2            =  0.625d0 * ((a3b2 - 3.0d0 * dmub3) * (a01 - dmu2) -  &
                          2.0d0 * dmu * a3 * (a0b1 - dmub0)) * delta
      eta3            = ((a2 * b3 - 3.0d0 * dmu * b2) * (a01 - dmu2) + &
                          dmu2 * (6.0d0 * a0b1 - 4.0d0 * a0 * b3 - 6.0d0 * dmub0)) * delta
!
      p1              = -a0 / alam1
      q1              =  0.3125d0 * (a01 / alams1 - 1.0d0)
      r1              = -1.5d0 * (a01 / alam1 - alam1) / a3
      p2              = -a0 / alam2
      q2              =  0.3125d0 * (a01 / alams2 - 1.0d0)
      r2              = -1.5d0 * (a01 / alam2 - alam2) / a3
!
      h1              = -(0.5d0 * eta0 - eta1 + eta2)
      h2              = -(-0.125d0 * eta0 + eta2 - eta3)
      h3              = -(0.5d0 * eta0 + eta1 + eta2) * dtr(i,1,k)
      h4              = -(-0.125d0 * eta0 + eta2 + eta3) * dtr(i,1,k)
!
      w11             =  0.5d0 - p1 + q1
      w12             = (0.5d0 + p1 + q1) * tr1
      w13             =  0.5d0 - p2 + q2
      w14             = (0.5d0 + p2 + q2) * tr2
      w21             =  -0.125d0 + q1 - r1
      w22             = (-0.125d0 + q1 + r1) * tr1
      w23             =  -0.125d0 + q2 - r2
      w24             = (-0.125d0 + q2 + r2) * tr2
!
      wa              =  w11 * w22 - w21 * w12
      wb              =  w14 * w23 - w24 * w13
      wc              =  w11 * w23 - w21 * w13
      wd              =  w11 * w24 - w21 * w14
      we              =  w12 * w23 - w22 * w13
      wf              =  w12 * w24 - w22 * w14
!
      det             =  1.0d0 / (2.0d0 * wa * wb - wc * wc + wd * wd + we * we - wf * wf)
      ya              = ( w22 * wb - w23 * wc + w24 * wd) * det
      yb              = (-w12 * wb + w13 * wc - w14 * wd) * det
      yc              = ( w23 * we - w24 * wf - w21 * wb) * det
      yd              = (-w13 * we + w14 * wf + w11 * wb) * det
      ye              = ( w21 * wc - w22 * we - w24 * wa) * det
      yf              = (-w11 * wc + w12 * we + w14 * wa) * det
      yg              = ( w23 * wa - w21 * wd + w22 * wf) * det
      yh              = (-w13 * wa + w11 * wd - w12 * wf) * det
!
      c1              =  ya * h1 + yb * h2 + yc * h3 + yd * h4
      d1              =  yc * h1 + yd * h2 + ya * h3 + yb * h4
      c2              =  ye * h1 + yf * h2 + yg * h3 + yh * h4
      d2              =  yg * h1 + yh * h2 + ye * h3 + yf * h4
!
      c1t             =  c1 * tr1
      c2t             =  c2 * tr2
      d1t             =  d1 * tr1
      d2t             =  d2 * tr2
!
      teu0            =  c1 + d1t + c2 + d2t + eta0
      teu1            =  p1 * (c1 - d1t) + p2 * (c2 - d2t) + eta1
      teu2            =  q1 * (c1 + d1t) + q2 * (c2 + d2t) + eta2
      teu3            =  r1 * (c1 - d1t) + r2 * (c2 - d2t) + eta3
!
      rdra(i,k,1)     = (teu0 + 2.0d0 * (teu1 + teu2)) * dmu
      rdra(i,k,2)     =  2.0d0 * (-0.125d0 * teu0 + teu2 + teu3) * dmu
!
      ted0            =  c1t + d1 + c2t + d2 + eta0 * dtr(i,1,k)
      ted1            =  p1 * (c1t - d1) + p2 * (c2t - d2) + eta1 * dtr(i,1,k)
      ted2            =  q1 * (c1t + d1) + q2 * (c2t + d2) + eta2 * dtr(i,1,k)
      ted3            =  r1 * (c1t - d1) + r2 * (c2t - d2) + eta3 * dtr(i,1,k)
!
      tdra(i,k,1)     = (ted0 - 2.0d0 * (ted1 - ted2))  * dmu
      tdra(i,k,2)     =  2.0d0 * (-0.125d0 * ted0 + ted2 - ted3) * dmu
!
      c1t             =  ya * tr1
      d1t             =  yc * tr1
      c2t             =  ye * tr2
      d2t             =  yg * tr2
!
      teu0            =  ya + d1t + ye + d2t
      teu1            =  p1 * (ya - d1t) + p2 * (ye - d2t)
      teu2            =  q1 * (ya + d1t) + q2 * (ye + d2t)
      teu3            =  r1 * (ya - d1t) + r2 * (ye - d2t)
      rdfa(i,k,1,1)   =  0.5d0 * teu0 + teu1 + teu2
      rdfa(i,k,2,1)   = -0.125d0 * teu0 + teu2 + teu3
!
      ted0            =  c1t + yc + c2t + yg
      ted1            =  p1 * (c1t - yc) + p2 * (c2t - yg)
      ted2            =  q1 * (c1t + yc) + q2 * (c2t + yg)
      ted3            =  r1 * (c1t - yc) + r2 * (c2t - yg)
      tdfa(i,k,1,1)   =  0.5d0 * ted0 - ted1 + ted2
      tdfa(i,k,2,1)   = -0.125d0 * ted0 + ted2 - ted3
!
      c1t             =  yb * tr1
      d1t             =  yd * tr1
      c2t             =  yf * tr2
      d2t             =  yh * tr2

      teu0            =  yb + d1t + yf + d2t
      teu1            =  p1 * (yb - d1t) + p2 * (yf - d2t)
      teu2            =  q1 * (yb + d1t) + q2 * (yf + d2t)
      teu3            =  r1 * (yb - d1t) + r2 * (yf - d2t)
      rdfa(i,k,1,2)   =  0.5d0 * teu0 + teu1 + teu2
      rdfa(i,k,2,2)   = -0.125d0 * teu0 + teu2 + teu3
!
      ted0            =  c1t + yd + c2t + yh
      ted1            =  p1 * (c1t - yd) + p2 * (c2t - yh)
      ted2            =  q1 * (c1t + yd) + q2 * (c2t + yh)
      ted3            =  r1 * (c1t - yd) + r2 * (c2t - yh)
      tdfa(i,k,1,2)   =  0.5d0 * ted0 - ted1 + ted2
      tdfa(i,k,2,2)   = -0.125d0 * ted0 + ted2 - ted3
!
! ...   all sky calculation
!
      if (cld(i,k) < cut)  then
        rdrc(i,k,1:2)     =  rdra(i,k,1:2)
        tdrc(i,k,1:2)     =  tdra(i,k,1:2)
        dtr(i,2,k)        =  dtr(i,1,k)
        rdfc(i,k,1:2,1:2) =  rdfa(i,k,1:2,1:2)
        tdfc(i,k,1:2,1:2) =  tdfa(i,k,1:2,1:2)
      else
        extopt        =  dble(tauc(i,k)) + extopta
        omarcs        =  max(dble(tauomc(i,k) + taur(i,k)), bound)
        ssalb         =  omarcs / extopt
        g             =  dble(tauomgc(i,k)) / omarcs
        g2            =  dble(tauomgc_4str(i,k,1)) / omarcs
        g3            =  dble(tauomgc_4str(i,k,2)) / omarcs
        f             =  dble(f2(i,k)) / omarcs
!
        ssalbf        =  ssalb * f
        om            =  ssalb * (1.0d0 - f) / (1.0d0 - ssalbf)
        tau           =  extopt * (1.0d0 - ssalbf)
        omega1        =  3.0d0 * om * (g - f) / (1.0d0 - f)
        omega2        =  5.0d0 * om * (g2 - f) / (1.0d0 - f)
        omega3        =  7.0d0 * om * (g3 - f) / (1.0d0 - f)
!
        a0            =  1.0d0 - om + eps_10
        a1            =  3.0d0 - omega1 + eps_10
        a2            =  5.0d0 - omega2 + eps_10
        a3            =  7.0d0 - omega3 + eps_10
        a01           =  a0 * a1
        a23           =  a2 * a3
        a03           =  a0 * a3
!
        b0            =  0.25d0 * om
        b1            = -0.25d0 * omega1 * rmu(i)
        b2            =  0.125d0 * omega2 * (3.0d0 * rmu2 - 1.0d0)
        b3            = -0.125d0 * omega3 * (5.0d0 * rmu3 - 3.0d0 * rmu(i))
        a3b2          =  a3 * b2
        a0b1          =  a0 * b1
        dmub0         =  dmu * b0
        dmub3         =  dmu * b3
!
        uu            =  a01 + d9 * (a23 + 4.0d0 * a03)
        uu2           =  uu * uu
        vv            =  4.0d0 * d9 * a01 * a23
!
        alams1        =  0.5d0 * (uu + sqrt(uu2 - vv))
        alams2        =  0.5d0 * (uu - sqrt(uu2 - vv))
        alam1         =  sqrt(alams1)
        alam2         =  sqrt(alams2)
!
        dtr(i,2,k)    =  exp(- tau * dmu)
        tr1           =  exp(- alam1 * tau)
        tr2           =  exp(- alam2 * tau)

        delta         =  1.0d0 / (9.0d0 * dmu4 - dmu2 * (9.0d0 * a01 + a23 + 4.0d0 * a03) + a01 * a23)
!
        eta0          = ((a1 * b0 - dmu * b1) * (a23 - 9.0d0 * dmu2) + &
                          2.0d0 * dmu2 * (a3b2 - 2.0d0 * a3 * b0 - 3.0d0 * dmub3)) * delta
        eta1          = ((a0b1 - dmub0) * (a23 - 9.0d0 * dmu2) - &
                          2.0d0 * dmu * a0 * (a3b2 - 3.0d0 * dmub3)) * delta
        eta2          =  0.625d0 * ((a3b2 - 3.0d0 * dmub3) * (a01 - dmu2) - &
                          2.0d0 * dmu * a3 * (a0b1 - dmub0)) * delta
        eta3          = ((a2 * b3 - 3.0d0 * dmu * b2) * (a01 - dmu2) + dmu2 * &
                         (6.0d0 * a0b1 - 4.0d0 * a0 * b3 - 6.0d0 * dmub0)) * delta
!
        p1            = -a0 / alam1
        q1            =  0.3125d0 * (a01 / alams1 - 1.0d0)
        r1            = -1.5d0 * (a01 / alam1 - alam1) / a3
        p2            = -a0 / alam2
        q2            =  0.3125d0 * (a01 / alams2 - 1.0d0)
        r2            = -1.5d0 * (a01 / alam2 - alam2) / a3
!
        h1            = -(0.5d0 * eta0 - eta1 + eta2)
        h2            = -(-0.125d0 * eta0 + eta2 - eta3)
        h3            = -(0.5d0 * eta0 + eta1 + eta2) * dtr(i,2,k)
        h4            = -(-0.125d0 * eta0 + eta2 + eta3) * dtr(i,2,k)
!
        w11           =  0.5d0 - p1 + q1
        w12           = (0.5d0 + p1 + q1) * tr1
        w13           =  0.5d0 - p2 + q2
        w14           = (0.5d0 + p2 + q2) * tr2
        w21           =  -0.125d0 + q1 - r1
        w22           = (-0.125d0 + q1 + r1) * tr1
        w23           =  -0.125d0 + q2 - r2
        w24           = (-0.125d0 + q2 + r2) * tr2
!
        wa            =  w11 * w22 - w21 * w12
        wb            =  w14 * w23 - w24 * w13
        wc            =  w11 * w23 - w21 * w13
        wd            =  w11 * w24 - w21 * w14
        we            =  w12 * w23 - w22 * w13
        wf            =  w12 * w24 - w22 * w14
!
        det           =  1.0d0 / (2.0d0 * wa * wb - wc**2 + wd**2 + we**2 - wf**2)
        ya            = ( w22 * wb - w23 * wc + w24 * wd) * det
        yb            = (-w12 * wb + w13 * wc - w14 * wd) * det
        yc            = ( w23 * we - w24 * wf - w21 * wb) * det
        yd            = (-w13 * we + w14 * wf + w11 * wb) * det
        ye            = ( w21 * wc - w22 * we - w24 * wa) * det
        yf            = (-w11 * wc + w12 * we + w14 * wa) * det
        yg            = ( w23 * wa - w21 * wd + w22 * wf) * det
        yh            = (-w13 * wa + w11 * wd - w12 * wf) * det
!
        c1            =  ya * h1 + yb * h2 + yc * h3 + yd * h4
        d1            =  yc * h1 + yd * h2 + ya * h3 + yb * h4
        c2            =  ye * h1 + yf * h2 + yg * h3 + yh * h4
        d2            =  yg * h1 + yh * h2 + ye * h3 + yf * h4
!
        c1t           =  c1 * tr1
        c2t           =  c2 * tr2
        d1t           =  d1 * tr1
        d2t           =  d2 * tr2
!
        teu0          =  c1 + d1t + c2 + d2t + eta0
        teu1          =  p1 * (c1 - d1t) + p2 * (c2 - d2t) + eta1
        teu2          =  q1 * (c1 + d1t) + q2 * (c2 + d2t) + eta2
        teu3          =  r1 * (c1 - d1t) + r2 * (c2 - d2t) + eta3
!
        rdrc(i,k,1)   = (teu0 + 2.0d0 * (teu1 + teu2)) * dmu
        rdrc(i,k,2)   =  2.0d0 * (-0.125d0 * teu0 + teu2 + teu3) * dmu
!
        ted0          =  c1t + d1 + c2t + d2 + eta0 * dtr(i,2,k)
        ted1          =  p1 * (c1t - d1) + p2 * (c2t - d2) +  eta1 * dtr(i,2,k)
        ted2          =  q1 * (c1t + d1) + q2 * (c2t + d2) + eta2 * dtr(i,2,k)
        ted3          =  r1 * (c1t - d1) + r2 * (c2t - d2) + eta3 * dtr(i,2,k)
!
        tdrc(i,k,1)   = (ted0 - 2.0d0 * (ted1 - ted2))  * dmu
        tdrc(i,k,2)   =  2.0d0 * (-0.125d0 * ted0 + ted2 - ted3) * dmu
!
        c1t           =  ya * tr1
        d1t           =  yc * tr1
        c2t           =  ye * tr2
        d2t           =  yg * tr2
!
        teu0          =  ya + d1t + ye + d2t
        teu1          =  p1 * (ya - d1t) + p2 * (ye - d2t)
        teu2          =  q1 * (ya + d1t) + q2 * (ye + d2t)
        teu3          =  r1 * (ya - d1t) + r2 * (ye - d2t)
        rdfc(i,k,1,1) =  0.5d0 * teu0 + teu1 + teu2
        rdfc(i,k,2,1) = -0.125d0 * teu0 + teu2 + teu3
!
        ted0          =  c1t + yc + c2t + yg
        ted1          =  p1 * (c1t - yc) + p2 * (c2t - yg)
        ted2          =  q1 * (c1t + yc) + q2 * (c2t + yg)
        ted3          =  r1 * (c1t - yc) + r2 * (c2t - yg)
        tdfc(i,k,1,1) =  0.5d0 * ted0 - ted1 + ted2
        tdfc(i,k,2,1) = -0.125d0 * ted0 + ted2 - ted3
!
        c1t           =  yb * tr1
        d1t           =  yd * tr1
        c2t           =  yf * tr2
        d2t           =  yh * tr2
!
        teu0          =  yb + d1t + yf + d2t
        teu1          =  p1 * (yb - d1t) + p2 * (yf - d2t)
        teu2          =  q1 * (yb + d1t) + q2 * (yf + d2t)
        teu3          =  r1 * (yb - d1t) + r2 * (yf - d2t)
        rdfc(i,k,1,2) =  0.5d0 * teu0 + teu1 + teu2
        rdfc(i,k,2,2) = -0.125d0 * teu0 + teu2 + teu3
!
        ted0          =  c1t + yd + c2t + yh
        ted1          =  p1 * (c1t - yd) + p2 * (c2t - yh)
        ted2          =  q1 * (c1t + yd) + q2 * (c2t + yh)
        ted3          =  r1 * (c1t - yd) + r2 * (c2t - yh)
        tdfc(i,k,1,2) =  0.5d0 * ted0 - ted1 + ted2
        tdfc(i,k,2,2) = -0.125d0 * ted0 + ted2 - ted3
      endif
    enddo
  enddo
!
  return
  end subroutine swtrlayer

 end subroutine swtran4st
!> \file
!> 4-stream spherical harmonic expansion for sw radiative transfer doubling/adding based
!! on Zhang & Li (JAS 2013), single layer solution based on Li & Ramaswamy (JAS 1996).
