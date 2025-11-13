  subroutine radiation(ccsd, ccsf, casd, casf, flxu, flxd, flxuc, flxdc, fsr, fst, fsrc, fstc, fsdb, fsfb, fdl, &
                       csdb, csfb, fdlc, csd, csf, shtj, t, q, o3, dz, & 
                       tsurf, salb, csal, ang, jday, ncldy, cldt, cld, clw_sub, cic_sub, rel_sub, rei_sub, &
                       exta, exoma, exomga, il1, il2, ilg, lay, lev)
      use psizes_19, only : nbs, nbl
      parameter (maxng = 10, nxloc = 100, ipph_re=0, &
                ipph=0, ioverlap=2, mcica=1, cldwatmin=0, ivers = 1, &
                idist=1, idimr = 17, max_sam = 1, ntile=1,           &
                iradforce = 0, irad_flux_profs = 1, iactive_rt = 1,  &
                itilerad=0, isvaod=1, ipam=1, imam=0, ncld=1, nstrm =4)

      implicit real (a-h,o-z), integer (i-n)
      !implicit none

      real ccsd(ilg), ccsf(ilg), casd(ilg), casf(ilg)
      real fsg(ilg), fsd(ilg), fsf(ilg), fsv(ilg), fsi(ilg),  &
           albpla(ilg), fdl(ilg), ful(ilg), hrs(ilg,lay), hrl(ilg,lay), &
           cst(ilg), csb(ilg), clt(ilg), clb(ilg), par(ilg), olr(ilg),  &
           olrc(ilg), csd(ilg), csf(ilg), fso(ilg), flg(ilg), fdlc(ilg),&
           fsr(ilg,lev,nbs), fst(ilg,lev,nbs), fsrc(ilg,lev,nbs), fstc(ilg,lev,nbs)
      real, dimension(ilg,lay) ::   hrsc        !<variable description\f$[units]\f$
      real, dimension(ilg,lay) ::   hrlc        !<variable description\f$[units]\f$
      real, dimension(ilg,lay,3) :: hrsh2o      !<variable description\f$[units]\f$
      real, dimension(ilg,lay,12) :: hrso3       !<variable description\f$[units]\f$
      real, dimension(ilg,lay,9) :: hrlh2o      !<variable description\f$[units]\f$
      real, dimension(ilg,lay,2) :: hrso2       !<variable description\f$[units]\f$
      real, dimension(ilg,lay,2) :: hrlco2      !<variable description\f$[units]\f$
      real, dimension(ilg,lay)   :: hrsco2      !<variable description\f$[units]\f$
      real, dimension(ilg,lay)   :: hrsch4      !<variable description\f$[units]\f$
      real, dimension(ilg,lay)   :: hrlo3       !<variable description\f$[units]\f$
      real, dimension(ilg,lay)   :: hrlch4      !<variable description\f$[units]\f$
      real, dimension(ilg,lay)   :: hrln2o      !<variable description\f$[units]\f$
      real, dimension(ilg,lay,nbs)   :: exta   !< Aerosol shortwave extinction optical depth
      real, dimension(ilg,lay,nbs)   :: exoma  !< Aerosol optical depth times SSA
      real, dimension(ilg,lay,nbs)   :: exomga !< Aerosol optical depth times asymmetry factor
      real, dimension(ilg,lay,nbs)   :: fa     !<
      real, dimension(ilg,lay,nbl)   :: absa   !< Aerosol longwave absorption coefficient

      real,  dimension(ilg,lay + 2) :: fsau !< all sky solar upward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: fsad !< all sky solar downward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: flau !< all sky thermal upward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: flad !< all sky thermal downward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: fscu !< clear sky solar upward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: fscd !< clear sky solar downward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: flcu !< clear sky thermal upward flux profile \f$[w/m^2]\f$
      real,  dimension(ilg,lay + 2) :: flcd !< clear sky thermal downward flux profile \f$[w/m^2]\f$


      real shtj(ilg,lev)  !<the pressure ratio at the upper level of a layer\f$[mb]\f$
      real tfull(ilg,lev) !<temperature at a layer upper level\f$[k]\f$
      real shj(ilg,lay)   !<the pressure ratio at the center of a layer\f$[mb]\f$
      real dshj (ilg,lay) !<difference between two sh\f$[mb]\f$
      real dz   (ilg,lay) !<layer geometeric thickness\f$[m]\f$
      real t    (ilg,lay) !<temperature at the middle of a layer\f$[k]\f$
      real o3   (ilg,lay) !<o3 mixing ratio\f$[gram/gram]\f$
      real q    (ilg,lay) !<h2o mixing ratio\f$[gram/gram]\f$
      real tsurf(ilg)     !<surface temperature
      real ang(ilg)       !<solar zenith angle
      integer jday (ilg)  !<Julian day
      real co2  (ilg,lay) !<co2 mixing ratio\f$[gram/gram]\f$
      real ch4  (ilg,lay) !<ch4 mixing ratio\f$[gram/gram]\f$
      real an2o (ilg,lay) !<n2o mixing ratio\f$[gram/gram]\f$
      real f11  (ilg,lay) !<f11 mixing ratio\f$[gram/gram]\f$
      real f12  (ilg,lay) !<f12 mixing ratio\f$[gram/gram]\f$
      real f113 (ilg,lay) !<f113 mixing ratio\f$[gram/gram]\f$
      real f114 (ilg,lay) !<f114 mixing ratio\f$[gram/gram]\f$
      real cld  (ilg,lay) !<cloud fraction\f$[0]\f$
!      real rh   (ilg,lay) !<relative humidity
      real anu  (ilg,lay) !<nu factor for gammar distribution used in cloud horizonatl variability\f$[0]\f$
      real eta  (ilg,lay) !<valume fraction of bc in cloud droplet (for semi-direct effect)\f$[m^3/m^3]\f$
      real tbnd (ilg,lay) !<temperature beyond the max/min set\f$[k]\f$
      real aerin(ilg,lay,idimr) !<input aerosol mixing ratio\f$[gram/gram]\f$
      real pressg(ilg)    !<surface presure\f$[pascal]\f$
      real gt    (ilg)    !<grand temperature\f$[k]\f$
      real o3top (ilg)    !<averaged o3 from toa to model top  level\f$[0]\f$
      real rmu   (ilg)    !<cosine of solar zenith ang ngle\f$[0]\f$
      real em0   (ilg)    !<surface emisivity\f$[0]\f$
      real cldt  (ilg)    !<total cloud fraction for a model grid\f$[0]\f$
      real radj  (ilg)    !<not used\f$[0]\f$
      real vtau  (ilg)    !<volcano aerosol optical depth\f$[0]\f$
      real trop  (ilg)    !<tropopause level\f$[mb]\f$
      real dp(ilg,lay)
      integer jlat   ! Deji  it is not been used.
      integer kount  !< Deji it not used Current model timestep \f$[unitless]\f$

      real, dimension(ilg,lay,nbs)   :: vstau !< stratospheric aerosol optical depth
      real, dimension(ilg,lay,nbs)   :: vsssa !< stratospheric aerosol single scattering albedo
      real, dimension(ilg,lay,nbs)   :: vsg   !< stratospheric aerosol assymetric factor
      real, dimension(ilg,lay,nbl)   :: vsabs !< stratospheric aerosol absorption coefficient
      real, dimension(ilg,lev,nbs) :: fsrsw, fsrcsw, fstsw, fstcsw

      real sw_ext_sa(ilg,lay,nbs) !<sw volcano aerosol extinction coeff.\f$[m^2/gram]\f$
      real sw_ssa_sa(ilg,lay,nbs) !<sw volcano aerosol single scattering albedo\f$[0]\f$
      real sw_g_sa  (ilg,lay,nbs) !<sw volcano aerosol asymmetry factor\f$[0]\f$
      real lw_abs_sa(ilg,lay,nbs) !<lw volcano aerosol absorptance\f$[m^2/gram]\f$
      real solv(12)               !<factor for solar vaiability\f$[0]\f$
      real fsdb (ilg,nbs)!<band mean direct downward flux at surface\f$[w/m^2]\f$
      real fsfb (ilg,nbs)!<band mean diffuse downward flux at surface\f$[w/m^2]\f$
      real csdb (ilg,nbs)!<band mean direct clear sky flux at surface\f$[w/m^2]\f$
      real csfb (ilg,nbs)!<band mean diffuse clear sky flux at surface\f$[w/m^2]\f$
      real fssb (ilg,nbs)!<band mean all sky downward at surface\f$[w/m^2]\f$
      real fsscb(ilg,nbs)!<band mean clear sky downward at surface\f$[w/m^2]\f$
      real wrka(ilg,nbs)!<Vertical integral of aerosol optical depth for each wavelength band
      real wrkb(ilg,nbs)!< Vertical integral of aerosol absorption depth for each wavelength band
      real, dimension(ilg,ntile) :: fsgt, fsdt, fsft, fsvt, fsit, &
                                    fdlt, flgt, fdlct, csbt, clbt, &
                                    csdt, csft
      real, dimension(ilg,ntile,nbs) :: fsdbt, fsfbt, csdbt,   &
                                        csfbt, fssbt, fsscbt

      real :: rrsq
      real rel_sub(ilg,lay,nxloc) !<effective radius of water cloud\f$[um]\f$
      real rei_sub(ilg,lay,nxloc) !<effective radius of ice cloud\f$[um]\f
      real clw_sub(ilg,lay,nxloc) !<liquid water content\f$[gram/m^3]\f$
      real cic_sub(ilg,lay,nxloc) !<ice water content\f$[gram/m^3]\f$
      integer ncldy(ilg)          !the number of subgrid with cloud<\f$[0]\f

      real salb(ilg,nbs)  !<all sky surface albedo\f$[0]\f$
      real csal(ilg,nbs)  !<clear sky surface albedo\f$[0]\f$

      real salbt(ilg,ntile,nbs) !<tiled all sky surface albedo\f$[units]\f$
      real csalt(ilg,ntile,nbs) !<tiled clear sky surface albedo\f$[units]\f$
      real faret(ilg,ntile)     !<fraction weight of each tile\f$[units]\f$
      real em0t (ilg,ntile)     !<tiled surface emisivity\f$[units]\f$
      real gtt  (ilg,ntile)     !<tiled grand temperature\f$[k]\f$
      logical :: lrefract     !< account for refraction in solar rt (.true.) or not (.false.)
      logical :: l_zero_lwsrc !< zero thermal sources (.true.) or not (.false.)

      integer, dimension(ilg,4) :: iseedrow     !<variable description\f$[units]\f$

      real, dimension(ilg,max_sam) :: colnum_sw, colnum_lw
      integer :: icolumn

      real, dimension(ilg,lay,nbs):: taucs, omcs, gcs, gcs2
      real, dimension(ilg,lay,nbl):: taucl, omcl, gcl, gcl2
      real lwf(ilg,lev), flxu(ilg,lev), flxd(ilg,lev), flxuc(ilg,lev), flxdc(ilg,lev)
      real flxdsw(ilg,lev), flxusw(ilg,lev), flxdswc(ilg,lev), flxuswc(ilg,lev)

      real fslo(ilg), fsamoon(ilg), flamoon(ilg)

      real p(ilg,lay), pfull(ilg,lev)
      logical lcsw, lclw

!      real wcdw(ilg,lay), wcdi(ilg,lay), wclw(ilg,lay), wcli(ilg,lay), &
!           radeqv(ilg,lay), radeqvw(ilg,lay), radeqvi(ilg,lay)

      real n2o_ppm
      common /radcon/ solar_c
      common /trace / rmco2, rmch4, rmn2o, rmo2, rmf11, rmf12, rmf113, rmf114

      data lcsw /.true./, lclw /.false./
      lrefract = .false.
      l_zero_lwsrc = .false.
      pi = 3.1415926
      rrsq = ((1. + 0.016 * cos((360/365.25*(jday(1)-4))*(pi/180)))** 2)/((1. - 0.016 ** 2) ** (3./2.))

      solv(:) = 0.

      do i = il1, il2
        tfull(i,lev) = tsurf(i)
        colnum_lw(i,max_sam) = 1
      enddo
    
      o3top(:) = 0.3853E-06

      vstau(:,:,:) = 0.0
      vsssa(:,:,:) = 0.0
      vsg(:,:,:) = 0.0
      vsabs(:,:,:) = 0.0
      absa(:,:,:)   = 0

      do k = 1, lay
        do i = il1, il2
          fa(i,k,:)   =  exomga(i,k,:) **2 / exoma(i,k,:)
        enddo
      enddo

      do i = il1, il2
        pressg(i) = 100.0
        gt(i)     = tfull(i,lev) 
        rmu(i)    = cos(ang(i) * 3.14159 / 180.0)
        em0(i) = 1.0
      enddo

      do k = 1, lay
        do i = il1, il2
          shj(i,k)   = sqrt(shtj(i,k) * shtj(i,k+1))
          if (k .eq. 1)                                             then
            tfull(i,k) = 0.5 * (3.0 * t(i,1) - t(i,2))
          else
            tfull(i,k) = 0.5 * (t(i,k-1) + t(i,k))
          endif
        enddo
      enddo
      do k = 1, lay
        do i = il1, il2
          x      =  0.01 * pressg(i)
          p(i,k) =  shj (i,k) * x
        enddo
      enddo

    do k = 1, lay
    kp1 = k + 1
    do i = il1, il2
      dp(i,k) =  0.0102 * pressg(i) * (shtj(i,kp1) - shtj(i,k))
    end do
    end do

      solar_c  = 1361.
      co2(:,:)  = 420. * 1.e-06     * 1.5188126
      ch4(:,:)  = 1.90 * 1.e-06      * 0.5522955
      an2o(:,:) = 0.33400 * 1.e-06  * 1.5188126
      f11(:,:)  = 0.28e-3 * 1.e-06    * 4.7418019  
      f12(:,:)  = 0.53e-3 * 1.e-06 * 4.1736279
      f113(:,:) = 0.!0.084e-3 * 1.e-06  * 5.2440456
      f114(:,:) = 0.!0.015e-3 * 1.e-06  * 5.8301691

      aerin(:,:,:) = 0.
      sw_ext_sa(:,:,:) = 0.
      sw_ssa_sa(:,:,:) = 0.
      sw_g_sa(:,:,:) = 0.




!////////////////////////For SW Calculations//////////////////////////////////
      call raddriv11(flxu, flxd, hrs, hrl, hrsc, hrlc, &
                     hrsh2o, hrso3, hrsco2, hrso2, hrsch4, &
                     hrlh2o, hrlo3, hrlco2, hrlch4, hrln2o, &
                     fsau,fsad,flau,flad, &
                     fscu,fscd,flcu,flcd, &
                     fsg, fsd, fsf, fsv, fsi, fdl, &
                     flg, fdlc, csb, clb, fsr, fst, fstc, fsrc, olr, olrc, & !
                     par, csd, csf, fslo, fsamoon, flamoon, fso, &
                     fsdb, fsfb, csdb, csfb, fssb, fsscb, &
                     wrka, wrkb, &
                     fsgt, fsdt, fsft, fsvt, fsit, &
                     fdlt, flgt, fdlct, csbt, clbt, &
                     csdt, csft, & !
                     fsdbt, fsfbt, csdbt, csfbt, fssbt, fsscbt, &
                     shtj, tfull, shj, dp, dz, t, o3, q, &
                     co2, ch4, an2o, f11, f12, anu, eta, imam, &
                     salb, csal, pressg, gt, o3top, rmu, solv, &
                     vstau, vsssa, vsg, vsabs, &
                     exta, exoma, exomga, fa, absa, &
                     em0, cldt, cld, tbnd, &
                     rel_sub, rei_sub, clw_sub, cic_sub, ncldy, &
                     salbt, csalt, em0t, gtt, faret, &
                     colnum_sw, colnum_lw, max_sam, iseedrow, &
                     lcsw, lclw, il1, il2, ilg, lay, lev, &
                     nxloc, ntile, cldwatmin, ivers, jlat, kount, &
                     nbl, nbs, rrsq, solar_c, pi, &
                     mcica, iradforce, irad_flux_profs, iactive_rt, &
                     itilerad, lrefract, l_zero_lwsrc, nstrm)

       		      

                      
                      
       do i = il1, il2
	 ccsd(i) = csd(i)
	 ccsf(i) = csf(i)
         casd(i) =  fsd(i)
         casf(i) =  fsf(i)
	 do k = 1, lev
	    flxdsw(i,k) = fsad(i,k+1)
	    flxusw(i,k) = fsau(i,k+1)
	    flxdswc(i,k) = fscd(i,k+1)
	    flxuswc(i,k) = fscu(i,k+1)
	    do n = 1, nbs
		fsrsw(i,k,n) = fsr(i,k,n)
	        fsrcsw(i,k,n) = fsrc(i,k,n)
		fstsw(i,k,n) = fst(i,k,n)
	        fstcsw(i,k,n) = fstc(i,k,n)
	    enddo 
	 enddo
       enddo 



!/////////////////////////For LW Calculations/////////////////////////////////
      l_zero_lwsrc = .true.
      lclw = .true.
      lcsw = .false.
      call raddriv11(flxu, flxd, hrs, hrl, hrsc, hrlc, &
                     hrsh2o, hrso3, hrsco2, hrso2, hrsch4, &
                     hrlh2o, hrlo3, hrlco2, hrlch4, hrln2o, &
                     fsau,fsad,flau,flad, &
                     fscu,fscd,flcu,flcd, &
                     fsg, fsd, fsf, fsv, fsi, fdl, &
                     flg, fdlc, csb, clb, fsr, fst, fstc, fsrc, olr, olrc, &
                     par, csd, csf, fslo, fsamoon, flamoon, fso, &
                     fsdb, fsfb, csdb, csfb, fssb, fsscb, &
                     wrka, wrkb, &
                     fsgt, fsdt, fsft, fsvt, fsit, &
                     fdlt, flgt, fdlct, csbt, clbt, &
                     csdt, csft, &
                     fsdbt, fsfbt, csdbt, csfbt, fssbt, fsscbt, &
                     shtj, tfull, shj, dp, dz, t, o3, q, &
                     co2, ch4, an2o, f11, f12, anu, eta, imam, &
                     salb, csal, pressg, gt, o3top, rmu, solv, &
                     vstau, vsssa, vsg, vsabs, &
                     exta, exoma, exomga, fa, absa, &
                     em0, cldt, cld, tbnd, &
                     rel_sub, rei_sub, clw_sub, cic_sub, ncldy, &
                     salbt, csalt, em0t, gtt, faret, &
                     colnum_sw, colnum_lw, max_sam, iseedrow, &
                     lcsw, lclw, il1, il2, ilg, lay, lev, &
                     nxloc, ntile, cldwatmin, ivers, jlat, kount, &
                     nbl, nbs, rrsq, solar_c, pi, &
                     mcica, iradforce, irad_flux_profs, iactive_rt, &
                     itilerad, lrefract, l_zero_lwsrc, nstrm)
           

!  fdlc is the solar doward flux from 4-100um infrared range
fsr = fsrsw 
fsrc = fsrcsw 
fst = fstsw 
fstc = fstcsw            
        do i = il1, il2
	  if (fdl(i) < 0) fdl(i) = 0.0 

          ccsf(i) = ccsf(i) + fdlc(i)!adding lw flux to diffuse for clear-sky
          casf(i) = casf(i) + fdl(i)
	  do k = 1, lev !adding lw flux to sw flux
	    if (flad(i,k+1) < 0) flad(i,k+1) = 0.0
!print*, flad(i,k+1), flxd(i,k)
	    flxd(i,k) = flxdsw(i,k) + flad(i,k+1)
	    flxu(i,k) = flxusw(i,k) + flau(i,k+1)
	    flxdc(i,k) = flxdswc(i,k) + flcd(i,k+1)
	    flxuc(i,k) = flxuswc(i,k) + flcu(i,k+1)
	    fsr(i,k,4) = fsrsw(i,k,4) + flau(i,k+1)
	    fsrc(i,k,4) = fsrcsw(i,k,4) + flcu(i,k+1)
	    fst(i,k,4) = fstsw(i,k,4) + flad(i,k+1)
	    fstc(i,k,4) = fstcsw(i,k,4) + flcd(i,k+1)

	  enddo       
	enddo       

   
    return
end subroutine radiation
