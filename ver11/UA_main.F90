   program UA
   use albedo_motran, only : albedom
   use psizes_19, only : nbs, nbl

   parameter (lay = 40, lev = lay + 1, ilg = 1, il1 = 1, il2 = 1, modlay=34, nxloc = 100)
   !   implicit real (a-h,o-z), integer (i-n)
   
   real, dimension (ilg,lev) :: flxu, flxd, flxuc, flxdc
   real, dimension (ilg,lev,4) ::  fsr, fsrc, fst, fstc
   real ccsd(ilg) !<surface clear-sky direct flux with LW
   real ccsf(ilg) !<surface clear-sy diffuse flux with LW
   real casd(ilg) !<surface all-sky direct flux with LW
   real casf(ilg) !<surface all-sky diffuse flux with LW
   real fdlc(ilg) !< lw surface downward flux clear-sky
   real fdl(ilg) !<lw surface downard flux all-sky
   real fsdb (ilg,nbs)!<band mean all-sky direct downward flux at surface\f$[w/m^2]\f$
   real fsfb (ilg,nbs)!<band mean all-sky diffuse downward flux at surface\f$[w/m^2]\f$
   real csdb (ilg,nbs)!<band mean clear-sky direct downward flux at surface\f$[w/m^2]\f$
   real csfb (ilg,nbs)!<band mean clear-sky diffuse downward flux at surface\f$[w/m^2]\f$
   real csd(ilg) !<SW total surface direct flux
   real csf(ilg) !<SW total surface diffuse flux
   real salb(ilg,nbs)    !<all sky surface albedo\f$[0]\f$
   real csal(ilg,nbs)    !<clear sky surface albedo\f$[0]\f$
   real shtj(ilg,lev)    !<the pressure ratio at the upper level of a layer\f$[mb]\f$
   real t    (ilg,lay)   !<temperature at the middle of a layer\f$[k]\f$
   real q    (ilg,lay)   !<h2o mixing ratio\f$[gram/gram]\f$
   real o3   (ilg,lay)   !<o3 mixing ratio\f$[gram/gram]\f$
   real rh   (ilg,lev)   !<relative humidity 100%
   real dz   (ilg,lay) !<layer geometeric thickness\f$[m]\f$
   real rel_in (ilg,lay)   !< effective radius of water cloud[um]
   real rei_in (ilg,lay)   !< effective radius of ice cloud[um]
   real clw_in (ilg,lay)   !< liquid water content [gram/m^3]
   real cic_in (ilg,lay)   !< ice water content [gram/^2]
   real rel_sub (ilg,lay,nxloc)   !< effective radius of water cloud[um]
   real rei_sub (ilg,lay,nxloc)   !< effective radius of ice cloud[um]
   real clw_sub (ilg,lay,nxloc)   !< liquid water content [gram/m^3]
   real cic_sub (ilg,lay,nxloc)   !< ice water content [gram/^3]
   real cld(ilg,lay)          !cloud fraction by layer<\f$[0]\f
   real cldt(ilg)          !<total cloud fraction
   real cod_vis(ilg)       !<cloud optical depth for vis over entire column
   integer ncldy(ilg)      !the number of subgrid with cloud<\f$[0]\f

   real hts  (ilg,lev)
   real ang(ilg)      !<solar zenith angle
   real tsurf(ilg)       !<surface temperature
   real exta(ilg,lay,nbs)   !< Aerosol shortwave extinction optical depth
   real exoma(ilg,lay,nbs)  !< Aerosol optical depth times SSA
   real exomga(ilg,lay,nbs) !< Aerosol optical depth times asymmetry factor
   real ssa_vis(ilg) !<Aerosol Single Scattering Albedo for band 1
   real asym_vis(ilg) !<Aerosol Asymmetry Factor for band 1
   real aod(ilg) !<AOD
   real waod(ilg)
   real vis(ilg) !visiblity based on aod
   real modaero(ilg,modlay) !<temporary aerosol profile based on MODTRAN code
   integer ialb(ilg)     !< surface scene types
   integer iaerotype(ilg)  !<Aerosol Type
   integer iseasn(ilg)    !<equal to 1 for spring/summer and 0 for fall/winter
   integer jday(ilg)     !<Julian Day


!//////////////////READ IN DATA////////////
   open (unit = 11, file ='input_file.dat',status = 'unknown')
   open (unit = 100, file ='output.dat',status = 'unknown')   
   do i = il1, il2
     read(11, *) ang(i), aod(i), ialb(i), iseasn(i), iaerotype(i), jday(i), cldt(i)
     do k = 1, lay
       read (11,*) hts(i,k),rh(i,k), shtj(i,k), t(i,k), q(i,k), o3(i,k), &
                   clw_in(i,k), cic_in(i,k), rel_in(i,k), rei_in(i,k), cld(i,k)
     enddo
       read (11,*) hts(i,lev),rh(i,lev), shtj(i,lev)
   enddo
   tsurf(:) = 300
!t(:,:) = 300.0
!///////////CALCULATE AEROSOL EXTINCTION, ABSORPTION, AND ASYMMETRY PROFILES///// 
  if (iaerotype(1) == 0) then !assuming pristine for entire file
     exta(:,:,:) = 0.0
     exoma(:,:,:) = 0.0
     exomga(:,:,:) = 0.0
     ssa_vis(:) = 0.0
     asym_vis(:) = 0.0
  else
    do i = il1,il2  
       vis(i)=getvis(aod(i),iseasn(i))
       do k= 1,modlay
         modaero(i,k)=aerprf(k,vis(i),iaerotype(i),iseasn(i))
       enddo
    enddo
    call aerobandpro(exta,exoma,exomga,ssa_vis,asym_vis, &
	             aod,ilg,modaero,lay,lev,modlay,iaerotype,hts,rh)
  endif

!////////////SURFACE ALEBDO/SCENCE TYPE PROCESSING////////////////////
   do i = il1, il2
     do j = 1, 4
       k =  ialb(i)
       salb(i,j) = albedom(j, K)
       csal(i,j) = albedom(j, K)
     enddo
  enddo

!////////////Use same cloud inputs for all subcolumns///////////////////
clw_sub(:,:,:) = 0.0
cic_sub(:,:,:) = 0.0
rel_sub(:,:,:) = 0.0
rei_sub(:,:,:) = 0.0

do i = il1, il2
  do k= 1,lay
    dz(i,k) = (hts(i,k) - hts(i,k+1))*1000.
    do j = 1, nxloc
!      if (j .le. cld(i,k)*nxloc) then
        clw_sub(i,k,j) = clw_in(i,k) 
        cic_sub(i,k,j) = cic_in(i,k) 
        rel_sub(i,k,j) = rel_in(i,k)
        rei_sub(i,k,j) = rei_in(i,k)
!      endif
    enddo
  enddo
  ncldy(i) = nint((cldt(i))*nxloc) 
enddo

!print*, rel_sub(1,36,:)
!clw_sub(:,:,:) = 0.0
!cic_sub(:,:,:) = 0.0
!rel_sub(:,:,:) = 0.0
!rei_sub(:,:,:) = 0.0
!cldt = 0.0
!exta(:,:,:) = 0.0
!exoma(:,:,:) = 0.0
!exomga(:,:,:) = 0.0
!///////////////////////////CALL MAIN CODE//////////////////////////////
   call radiation(ccsd, ccsf, casd, casf, flxu, flxd, flxuc, flxdc, fsr, fst, fsrc, fstc, fsdb, fsfb, fdl, &
		  csdb, csfb, fdlc, csd, csf, shtj, t, q, o3, dz, &
		  tsurf, salb, csal, ang, jday, ncldy, cldt, cld, clw_sub, cic_sub, rel_sub, rei_sub, &
		  exta, exoma, exomga, il1, il2, ilg, lay, lev)

!/////////////////////////WRITE OUTPUT/////////////////////////////////
  do i = il1, il2
	write(100,*) '     Lev(mb)','   ALL-SWUP VIS', '  ALL-SWUP NIR', '  ALL-SWDN VIS', '  ALL-SWDN NIR',&
			'  CLR-SWUP VIS', '  CLR-SWUP NIR', '  CLR-SWDN VIS', '  CLR-SWDN NIR'
     do k = 1, lev
       write(100, 2023) shtj(i,k), fsr(i,k,1), sum(fsr(i,k,2:4)), fst(i,k,1), sum(fst(i,k,2:4)), &
				   fsrc(i,k,1), sum(fsrc(i,k,2:4)), fstc(i,k,1), sum(fstc(i,k,2:4))
     enddo

	write(100,*) '------- All-Sky:'
	write(100,*) 'Surface  Direct(Wm2):  ', 'Surface Diffuse(Wm2):'
	write(100, 2024)  casd(i), casf(i)

	write(100,*) '----- Clear-Sky:'
	write(100,*) 'Surface  Direct(Wm2):  ', 'Surface Diffuse(Wm2):'
	write(100, 2024)  ccsd(i), ccsf(i)
	
	!write(100,*) 'Cloud Optical Depth Band 1:'
	!write(100, 2026) cod_vis(i)
        write(100,*) '-----------------------------------------------------------------------------------------------------------'
   enddo
   


!print*, cod_vis, casd+casf
!print*, casd(1)+casf(1)
 2023 format (9(F12.3, 2x))
 2024 format (2(F21.3, 2x))
 2025 format (5(F12.3, 2x))
 2026 format (1(F12.3, 2x))

   stop
   end  program UA
