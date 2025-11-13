subroutine aerobandpro(exta,exoma,exomga,ssa_vis, asym_vis, &
	               aod,ilg,modaero,lay, lay2, modlay,iaerotype,hts,rh)
      use psizes_19, only : nbs
      !parameter(lay2=41)

      implicit real (a-h,o-z), integer (i-n)

  real exta(ilg,lay,nbs) 
  real exoma(ilg,lay,nbs) 
  real exomga(ilg,lay,nbs) 
  
  real hts(ilg,lay2)
  real rh(ilg,lay2)   
  integer iaerotype(ilg)
  real modaero(ilg,modlay)

  integer :: n
  integer :: ii
  integer :: k
  integer :: ht1
  integer :: ht2
  real  :: ext1
  real  :: ext2
  real , dimension(lay2) :: aeroprof2

  real dis
  real aod1
  real aod2
  real aod_bdry
  real sca
  real ssa_vis(ilg)
  real asym_vis(ilg)
  real aod(ilg)
  real od550(lay2)

  real wspd
  integer nwspd
  integer iwav55
  integer wavel
  real extc(788)
  real absc(788)
  real asymc(788)
  real slope
  real ext55
  real wind(4)

  real , dimension(788,4) :: boundext
  real , dimension(788,4) :: boundabs
  real , dimension(788,4) :: boundsym
  ! 788 spectral properties at each layer with different RH range for each number at end of variable
  real extprof(lay2,788)
  real scatprof(lay2,788)
  real symprof(lay2,788)
  ! Four bands with X1,x2,x3,x4 with four RH ranges
  real xext(lay2,4)
  real xsct(lay2,4)
  real xsym(lay2,4)
  real yy(nbs)

  real rhfrac
  real extint(lay2,nbs)
  real absorb(lay2,nbs)
  real asym(lay2,nbs)

! Loading AEROSOL SPECTRAL (.2 -300 um) Extinction, Absorption, and Asymmetry 
  integer I_HGT
  COMMON/PRFD/I_HGT(34)
  REAL VX0
  COMMON/EXTWAV/VX0(788)
   REAL RUREXT,RURABS,RURSYM,URBEXT,URBABS,URBSYM,                   &
     &  OCNEXT,OCNABS,OCNSYM
  COMMON/EXTD1/RUREXT(788,4),RURABS(788,4),RURSYM(788,4),  &
     &             URBEXT(788,4),URBABS(788,4),URBSYM(788,4),  &
     &             OCNEXT(788,4),OCNABS(788,4),OCNSYM(788,4)
   REAL TROEXT,TROABS,TROSYM,                                        &
     &  FG1EXT,FG1ABS,FG1SYM,FG2EXT,FG2ABS,FG2SYM,BSTEXT,BSTABS,BSTSYM, &
     &  AVOEXT,AVOABS,AVOSYM,FVOEXT,FVOABS,FVOSYM,DMEEXT,DMEABS,DMESYM
  COMMON/EXTD2/TROEXT(788,4),TROABS(788,4),TROSYM(788,4),  &
     &             FG1EXT(788),FG1ABS(788),FG1SYM(788),        &
     &             FG2EXT(788),FG2ABS(788),FG2SYM(788),        &
     &             BSTEXT(788),BSTABS(788),BSTSYM(788),        &
     &             AVOEXT(788),AVOABS(788),AVOSYM(788),        &
     &             FVOEXT(788),FVOABS(788),FVOSYM(788),        &
     &             DMEEXT(788),DMEABS(788),DMESYM(788)
   REAL DSEXT,DEABS,DSASYM
  COMMON/DESAER/DSEXT(788,4),DSABS(788,4),DSASYM(788,4)
   REAL WT
  COMMON/BANDWT/WT(788)

do n=1,ilg
!Linear Interpolation of New MERRA2 Hieghts from MODTRAN Heights
  do k = 1,lay2 
    if (hts(n,k).ge. 1.0e12) then
      aeroprof2(k)=0.
    else
    do ii = 1,modlay
       if (hts(n,k).ge.I_HGT(ii)) then
          ht1=I_HGT(ii)
          ht2=I_HGT(ii+1)
          ext1=modaero(n,ii)
          ext2=modaero(n,ii+1)
       endif
    enddo
       aeroprof2(k)= ext1 +(hts(n,k)-ht1)*((ext2-ext1)/(ht2-ht1))
    endif
  enddo


!Scale ext profile to double check it matches AOD
aod1=0. !below 4km
aod2=0. !above 4km
do k=1,lay2
   if (k .ge. 2) then
    dis=hts(n,k-1)-hts(n,k)
    od550(k)= dis*aeroprof2(k)
   endif
   if (hts(n,k) .gt. 4.) then
   aod2=aod2+od550(k)
   else 
   aod1=aod1+od550(k)
  endif
enddo

aod_bdry=aod(n)-aod2
if (aod_bdry .gt. 0) then
  sca=aod_bdry/aod1
  do k=1,lay2
    if (hts(n,k) .le. 4.) then
      od550(k)=od550(k)*sca
    endif
    if (k .ge. 2) then
      dis=hts(n,k-1)-hts(n,k)
      aeroprof2(k)=od550(k)/dis
    endif
  enddo
else ! for really small AOD values.. i.e. below upper atmosphere set extinction/AOD values
 sca=aod(n)/(aod1+aod2)
  do k=1,lay2
     od550(k)=od550(k)*sca
    if (k .ge. 2) then
     dis=hts(n,k-1)-hts(n,k)
     aeroprof2(k)=od550(k)/dis
    endif
  enddo
endif

aod1=0
do k=2,lay2
 dis = hts(1,k-1)-hts(1,k)
 sca = dis *aeroprof2(k)
 aod1 = aod1+sca
enddo
!print*,aod(n),aod1

!Set Boundray layer aerosol type
if (iaerotype(n).eq. 0) then ! no aerosols 
   boundext(:,:)=0.
   boundabs(:,:)=0.
   boundsym(:,:)=0.
elseif (iaerotype(n).eq. 1) then ! Ocean
   boundext=OCNEXT
   boundabs=OCNABS
   boundsym=OCNSYM
elseif (iaerotype(n).eq. 2) then ! Rural
   boundext=RUREXT
   boundabs=RURABS
   boundsym=RURSYM
elseif (iaerotype(n).eq. 3) then !Urban
   boundext=URBEXT
   boundabs=URBABS
   boundsym=URBSYM
elseif (iaerotype(n).eq. 4) then ! Desert
   !Note 4 values represent different wind speeds(0,10,20, and 30 m/s) but only 10 m/s is used below
      !INTERPOLATE THE RADIATIVE PROPERTIES AT WIND SPEED WSPD FOR DESERT AEROSOLS ONLY
      WIND=(/ 0.,10.,20.,30. /)
      WSPD=10 !Wind Speed
      NWSPD=INT(WSPD/10)+1
      IWAV55=0
      SLOPE=(WSPD-WIND(NWSPD))/(WIND(NWSPD+1)-WIND(NWSPD))
      DO WAVEL=1,788
          IF(VX0(WAVEL).LE.0.55)IWAV55=WAVEL
     !  EXTINCTION COEFFICIENT:
          EXTC(WAVEL)=DSEXT(WAVEL,NWSPD)                              &
     &      *(DSEXT(WAVEL,NWSPD+1)/DSEXT(WAVEL,NWSPD))**SLOPE
     !  ABSORPTION COEFFICIENT:
          ABSC(WAVEL)=DSABS(WAVEL,NWSPD)                              &
     &      *(DSABS(WAVEL,NWSPD+1)/DSABS(WAVEL,NWSPD))**SLOPE
     !  ASYMMETRY PARAMETER:
          ASYMC(WAVEL)=DSASYM(WAVEL,NWSPD)                             &
     &      *(DSASYM(WAVEL,NWSPD+1)/DSASYM(WAVEL,NWSPD))**SLOPE
      ENDDO
     ! Determine 550nm Extinction
     SLOPE=(.55-VX0(IWAV55))/(VX0(IWAV55+1)-VX0(IWAV55))
     EXT55=EXTC(IWAV55)
     IF(EXT55.GT.0. .AND. EXTC(IWAV55+1).GT.0.)THEN
         EXT55=EXT55*(EXTC(IWAV55+1)/EXT55)**SLOPE
     ELSE
         EXT55=EXT55+SLOPE*(EXTC(IWAV55+1)-EXT55)
     ENDIF
     ! Normalize extinction and absorption by dividing 550 extinction
     DO WAVEL=1,788
         EXTC(WAVEL)=EXTC(WAVEL)/EXT55
         ABSC(WAVEL)=ABSC(WAVEL)/EXT55
     ENDDO
   DO ii=1,4
     boundext(:,ii)=EXTC
     boundabs(:,ii)=ABSC
     boundsym(:,ii)=ASYMC
   ENDDO
elseif (iaerotype(n).eq. 5) then !Tropospheric
   boundext=TROEXT
   boundabs=TROABS
   boundsym=TROSYM
else 
   boundext=TROEXT
   boundabs=TROABS
   boundsym=TROSYM
endif


!Make Spectral Extc/Absorp/Asym Profiles
!Extinction and absorption data normalized by dividing both by 550nm Extinction
do k = 1,lay2
 do ii =1,394
   if (hts(n,k).le. 2) then
     if (iaerotype(n) .eq. 4. .or. iaerotype(n) .eq. 0) then     
       extprof(k,ii) = aeroprof2(k) * boundext(ii,2)
       scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * boundabs(ii,2)
       symprof(k,ii) = scatprof(k,ii) * boundsym(ii,2)
     else
       if (rh(n,k) .lt. 70.) then ! Ocean, Rural, Urban,
         rhfrac = (4.6051702 - log(100 - rh(n,k)))/1.2039728
         extprof(k,ii) = aeroprof2(k) * boundext(ii,1) * exp(rhfrac * log(boundext(ii,2)/boundext(ii,1)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * boundabs(ii,1) * exp(rhfrac * log(boundabs(ii,2)/boundabs(ii,1)))
         symprof(k,ii) = scatprof(k,ii) * boundsym(ii,1) * exp(rhfrac * log(boundsym(ii,2)/boundsym(ii,1)))
       elseif (rh(n,k) .lt. 80.) then
         rhfrac = (3.4011974 - alog(100 - rh(n,k)))/0.40546511
         extprof(k,ii) = aeroprof2(k) * boundext(ii,2) * exp(rhfrac * log(boundext(ii,3)/boundext(ii,2)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * boundabs(ii,2) * exp(rhfrac * log(boundabs(ii,3)/boundabs(ii,2)))
         symprof(k,ii) = scatprof(k,ii) * boundsym(ii,2) * exp(rhfrac * log(boundsym(ii,3)/boundsym(ii,2)))
       else
         if (rh(n,k) .ge. 99.) then
           rhfrac=1.
         else
           rhfrac = (2.9957323 - log(100 - rh(n,k)))/2.9957323
         endif
         extprof(k,ii) = aeroprof2(k) * boundext(ii,3) * exp(rhfrac * log(boundext(ii,4)/boundext(ii,3)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * boundabs(ii,3) * exp(rhfrac * log(boundabs(ii,4)/boundabs(ii,3)))
         symprof(k,ii) = scatprof(k,ii) * boundsym(ii,3) * exp(rhfrac * log(boundsym(ii,4)/boundsym(ii,3)))
       endif
     endif 
   elseif (hts(n,k) .le. 11) then
     if (iaerotype(n).eq. 4. .or. iaerotype(n) .eq. 0) then
       extprof(k,ii) = aeroprof2(k) * TROEXT(ii,2)
       scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * TROABS(ii,2)
       symprof(k,ii) = scatprof(k,ii) * TROSYM(ii,2)
     else
       if (rh(n,k) .lt. 70.) then ! Ocean, Rural, Urban,
         rhfrac = (4.6051702 - log(100 - rh(n,k)))/1.2039728
         extprof(k,ii) = aeroprof2(k) * TROEXT(ii,1) * exp(rhfrac * log(TROEXT(ii,2)/TROEXT(ii,1)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * TROABS(ii,1) * exp(rhfrac * log(TROABS(ii,2)/TROABS(ii,1)))
         symprof(k,ii) = scatprof(k,ii) * TROSYM(ii,1) * exp(rhfrac * log(TROSYM(ii,2)/TROSYM(ii,1)))
       elseif (rh(n,k) .lt. 80.) then
         rhfrac = (3.4011974 - alog(100 - rh(n,k)))/0.40546511
         extprof(k,ii) = aeroprof2(k) * TROEXT(ii,2) * exp(rhfrac * log(TROEXT(ii,3)/TROEXT(ii,2)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * TROABS(ii,2) * exp(rhfrac * log(TROABS(ii,3)/TROABS(ii,2)))
         symprof(k,ii) = scatprof(k,ii) * TROSYM(ii,2) * exp(rhfrac * log(TROSYM(ii,3)/TROSYM(ii,2)))
       else
         if (rh(n,k) .ge. 99.) then
           rhfrac=1.
         else
           rhfrac = (2.9957323 - log(100 - rh(n,k)))/2.9957323
         endif
         extprof(k,ii) = aeroprof2(k) * TROEXT(ii,3) * exp(rhfrac * log(TROEXT(ii,4)/TROEXT(ii,3)))
         scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * TROABS(ii,3) * exp(rhfrac * log(TROABS(ii,4)/TROABS(ii,3)))
         symprof(k,ii) = scatprof(k,ii) * TROSYM(ii,3) * exp(rhfrac * log(TROSYM(ii,4)/TROSYM(ii,3)))
       endif  
     endif
   elseif (hts(n,k) .le. 100) then
     extprof(k,ii) = aeroprof2(k) * bstext(ii)
     scatprof(k,ii) = extprof(k,ii) - aeroprof2(k) * bstabs(ii)
     symprof(k,ii) = scatprof(k,ii) * bstsym(ii)
   endif
 enddo
enddo


!Convert to 4 bands
xext(:,:) = 0.0
xsct(:,:) = 0.0
xsym(:,:) = 0.0
yy(:) = 0.0
do k = 1,lay2 
   do ii = 1, 394 
     if (VX0(ii).le. 0.689) then
        xext(k,1) = xext(k,1) + extprof(k,ii) * wt(ii)
        xsct(k,1) = xsct(k,1) + scatprof(k,ii) * wt(ii)
        xsym(k,1) = xsym(k,1) + symprof(k,ii) * wt(ii)
        if (k == 1) yy(1) = yy(1) + wt(ii)
     endif
     if (VX0(ii).gt. 0.689 .and. VX0(ii) .le. 1.19) then
        xext(k,2) = xext(k,2) + extprof(k,ii) * wt(ii)
        xsct(k,2) = xsct(k,2) + scatprof(k,ii) * wt(ii)
        xsym(k,2) = xsym(k,2) + symprof(k,ii) * wt(ii)
        if (k == 1) yy(2)  = yy(2) + wt(ii)
     endif
     if (VX0(ii) .gt. 1.19 .and. VX0(ii) .le. 2.38) then 
        xext(k,3) = xext(k,4) + extprof(k,ii) * wt(ii)
        xsct(k,3) = xsct(k,4) + scatprof(k,ii) * wt(ii)
        xsym(k,3) = xsym(k,4) + symprof(k,ii) * wt(ii)
        if (k == 1) yy(3)  = yy(3) + wt(ii)
     endif
     if (VX0(ii) .gt. 2.89 .and. VX0(ii) .le. 4.) then  
        xext(k,4) = xext(k,4) + extprof(k,ii) * wt(ii)
        xsct(k,4) = xsct(k,4) + scatprof(k,ii) * wt(ii)
        xsym(k,4) = xsym(k,4) + symprof(k,ii) * wt(ii)
        if (k == 1) yy(4)  = yy(4) + wt(ii)
     endif
   enddo
enddo

!check AOD model calculations agianst AERONET measurements
!wavelength vs index:
!1040nm> 163;   870nm> 135;   675nm> 96;   550nm> 71;   500nm> 61;   440nm> 49;   380nm> 37;   340nm> 29; 
!aod1=0
!do k=2,lay2
! dis = hts(1,k-1)-hts(1,k)
! sca = dis *extprof(k,71) !71 corresponds to 55nm
! aod1 = aod1+sca
!enddo
!print*,aod1



ssa_vis(n) = 0! (xsct(38,1)/xext(38,1))
asym_vis(n) = 0!(xsym(38,1)/xext(38,1))
ii=0
do k = 1, lay2
  if (hts(n,k) .le. 2.) then
     ssa_vis(n) = ssa_vis(n) + (xsct(k,1)/xext(k,1))
     asym_vis(n) = asym_vis(n) + (xsym(k,1)/xext(k,1))
     ii = ii + 1
  endif
enddo
ssa_vis(n) = ssa_vis(n)/ii
asym_vis(n) = asym_vis(n)/ii

!print*, ssa_vis, asym_vis

!Get Optical Depth(Exta), OD*SSA(Exoma), OD*Assym(Exomga)
do k = 2, lay2 
    do i=1,4
       x = (hts(n,k-1) - hts(n,k)) / yy(i)
       exta(n,k-1,i)   = x * xext(k,i) 
       exoma(n,k-1,i)  = x * xsct(k,i) 
       exomga(n,k-1,i) = x * xsym(k,i) 
   enddo
enddo
if (iaerotype(n) .eq. 0) then
  exta(:,:,:) = 0.
  exoma(:,:,:) = 0.
  exomga(:,:,:) = 0.
endif
!print*,sum(exta(n,:,1)),aod(n)
!print*,exta,exoma,exomga
enddo! n/ilg

end subroutine aerobandpro
