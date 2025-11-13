PRO MODEL_INPUT
station='SGP'

if station eq "SGP" then begin
 lat_stat=36.5 & lon_stat=-97.5
 sctyp=50 & AEROTYPE=2
 dst=9. & dend=13.
endif
if station eq 'BAR' then begin
 lat_stat=71.5 & lon_stat=-156.5
 sctyp=60  & AEROTYPE=1 ;60 or 46
 dst=18 & dend=20 
endif
if station eq 'TAM' then begin
 lat_stat=22.5 & lon_stat=5.5
 sctyp=56 & AEROTYPE=4
 dst=15. & dend=18.;18.
endif
if station eq 'MNM' then begin
 lat_stat=24.5 & lon_stat=153.5
 sctyp= 57 & AEROTYPE=1
 dst=20. & dend=22.
endif
if station eq 'SYO' then begin
 lat_stat=-69.5 & lon_stat=39.5
 sctyp= 55 & AEROTYPE=1
 dst=2. & dend=8.
endif
if station eq 'PVC' then begin
  lat_stat=42.5 & lon_stat=-70.5
  sctyp=57 & AEROTYPE=1
  dst=4. & dend=5.;1-2  &  4-5 for 062013
stat=station
endif
if station eq 'ENA' then begin
 lat_stat=39.5 & lon_stat=-28.5
 sctyp=57 & AEROTYPE=1
 dst=19. & dend=19. ;13-13 for 032017 or 19-19 for082018
 ;!!!MAKE SURE TO CHANGE MODIS FILE DEPENDING ON YEAR FOR ENA
endif
;AEROTYPE=0
for day=dst,dend do begin
;//////////////////////LOAD MERRA2//////////////////
path='/home/jbrendecke/rad_stations/'+station+'/MERRA2/
;DATA: https://disc.gsfc.nasa.gov/datasets/M2I3NPASM_5.12.4/summary
dd=strcompress(string(fix(dst+(day-dst))),/remove_all)
f0=findfile(path+'*'+dd+'.SUB.nc')
cdfid=ncdf_open(f0[dst-dst])
  ncdf_varget,cdfid,'time',time
  ncdf_varget,cdfid,'lon',lon
  ncdf_varget,cdfid,'lat',lat
  ncdf_varget,cdfid,'lev',lev
  ncdf_varget,cdfid,'H',height
  ncdf_varget,cdfid,'QV',qv
  ncdf_varget,cdfid,'PHIS',elv
  ncdf_varget,cdfid,'PS',ps
  ncdf_varget,cdfid,'SLP',slp
  ncdf_varget,cdfid,'T',temperature
  ncdf_varget,cdfid,'RH',relh
  ncdf_varget,cdfid,'O3',o3
ncdf_close, cdfid
nlev=n_elements(lev)


;/////////////////LOAD MODIS AOD///////////////////
 pathc='/home/jbrendecke/rad_stations/'+station+'/ceres/'
 file_cer=findfile(pathc+'*.nc')

 cdfid=ncdf_open(file_cer[1])
    var_id=ncdf_varid(cdfid,'time')
    ncdf_varget,cdfid,var_id, time_cer
    var_id=ncdf_varid(cdfid,'lon')
    ncdf_varget,cdfid,var_id,lon_cer
    var_id=ncdf_varid(cdfid,'lat')
    ncdf_varget,cdfid,var_id,lat_cer
    var_id=ncdf_varid(cdfid,'ini_aod55_1h')
    ncdf_varget,cdfid,var_id,aod_modis
 ncdf_close, cdfid
 time_cer=time_cer-floor(time_cer[0])+1
 indt=where(time_cer ge dst and time_cer lt dend+1)
 indlat=where(lat_cer eq lat_stat)
 if lon_stat lt 0 then lon_stat2=lon_stat+360. else lon_stat2=lon_stat
 indlon=where(lon_cer eq lon_stat2)
 AOD_DAY=aod_modis[indlat,indlon,*]
 AOD_DAY=AOD_DAY[indt]

;~~~~~~~~~~~~~~~~~~~~~~~~~
;station lat/lon
londiff=min(abs(lon-lon_stat),nnn)
ilon=nnn;where(lon eq lon_stat2,/null)
londiff=min(abs(lat-lat_stat),nnn)
ilat=nnn;where(lat eq lat_stat,/null)

;;if using AERONET AOD
if station eq 'BARq' then begin
  patha='/home/jbrendecke/rad_stations/BAR/aeronet_aod/continous/
  filea=findfile(patha+'*.csv')
  dat2=read_csv(filea[0],header=1,n_table_header=1,table_header=tabvari2,missing_value=-9999)
  timeaod=dat2.field1
  aod_a=dat2.field2
  ind=where(timeaod ge dst and timeaod lt dend+1)
  AOD_DAY=aod_a[ind]
endif
if station eq "ENA" then begin
  filea=findfile('/home/jbrendecke/rad_stations/ENA/aeronet_aod/*.dat')
  restore,filea[0]
  timeaod=fdata.time
  aod_a=fdata.aod550
  ind=where(timeaod ge dst and timeaod lt dend+1)
  AOD_DAY=aod_a[ind]
  AOD_DaY[where(finite(aod_day) eq 0)]=mean(AOD_DAY,/nan)
endif

ind=where(AOD_DAY lt 0.0421,/null);because MODTRAN cant handle really small AOD values
AOD_DAY[ind]=0.042

;//////////RUN CKD CODE FOR EACH HOUR//////////////
for ihr=0,23 do begin;0,23
    hour=0.5
    hour=hour+ihr

    timeind=[0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7]
    itime=timeind[floor(hour)]
	
    ;get julian day
    hr=round(floor(hour))
    minute=round((hour-hr)*60.)
    yymmdd=strmid(f0,strlen(f0)-15,8)
    test=float(yymmdd);check point to make sure yymmdd are all numbers
    if hour lt 10 then utc='0'+strcompress(string(hr)+string(minute)+'00',/remove_all) else $
      utc=strcompress(string(hr)+string(minute)+'00',/remove_all)
    
    set_jday,fix(strmid(yymmdd,0,4)),fix(strmid(yymmdd,4,2)),fix(strmid(yymmdd,6,2)),jday0
    jday=jday0+(hour/24.)

    ;get SZA
    sza1=zenith((jday0-1+(hour/24.)),lat_stat,lon_stat)*(180./!pi)



    if sza1 lt 90. then begin ;CKD code spits out "******" for night time values
      
       AOD=aod_day[ihr+((day-dst)*24.)]	

       if jday gt 75 and jday lt 265 then SUMMER=string(1,format='(i10)') else SUMMER=string(0,format='(i10)')
	;-------------------------------------------------------------------------
        ;//////////Write out MERRA2 Profile///////////	
	;Pressure, Hieght, RH should have 41 levels
	;Temperature, QV, O3 should have 40 levels
	;variables for specific location/time
	height_ = reform(height[ilon,ilat, *,itime])/1000.
	indx1 = where(height_ lt 200.,/null,len1)
	height_ =height_[indx1]
	lev_ = lev[indx1]
	temp_ = reform(temperature[ilon,ilat, *, itime])
	indx2 = where(temp_ lt 350.,/null,len2)
	temp_ = temp_[indx2]
	O3_ = reform(O3[ilon,ilat, indx2,itime])
 	QV_ = reform(QV[ilon, ilat, indx2, itime])
	relh_ = reform(relh[ilon,ilat, indx2,itime])*100.
	elv_ = reform(elv[ilon,ilat,itime])/(1000.*9.8)
	ps_ = reform(ps[ilon,ilat,itime])/100.
	slp_ = reform(slp[ilon,ilat,itime])/100.
	if station eq 'SYO' then elv_=0.03
	if elv_ lt .01 then elv_=0.0

	;certain cases the height/pres has more values than other variables
	if len1 ne len2 then begin
	   delh = height_[1]-height_[0]
	   temp0 = temp_[0]+delh *6.5
	   relh0 = relh_[0]
	   o30 = o3_[0]
	   qv0 = qv_[0];(.611*exp((17.27*temp0)/

	   temp_= [temp0,temp_]
	   relh_ = [relh0, relh_]
	   o3_ = [o30, o3_]
	   qv_ = [qv0, qv_]
	endif
	;make sure all time steps start at same height
	if ps_ gt lev_[0] then begin
	  delh = height_[0] - elv_
	  ht0 = elv_
	  lev0 = ps_
	  if station eq 'SYO' then  lev0= slp_* ((288.-6.5*elv_)/288.)^5.2561;
	  temp0 = temp_[0]+delh*6.5
	  relh0 = relh_[0]	  
   	  o30 = o3_[0]
	  qv0 = qv_[0]
	
	  height_ = [ht0,height_]
	  lev_ = [lev0, lev_]
	  temp_= [temp0,temp_]
	  relh_ = [relh0, relh_]
	  o3_ = [o30, o3_]
	  qv_ = [qv0, qv_]
	endif else begin
	  delh = elv_ - height_[0]
	  height_[0] = elv_
	  lev_[0] = ps_
	  temp_[0] = temp_[0] - delh*6.5
	endelse

        if n_elements(height_) ne n_elements(lev_) or n_elements(height_) ne n_elements(temp_)$
	or n_elements(height_) ne n_elements(o3_) or n_elements(height_) ne n_elements(relh_)$
	or n_elements(height_) ne n_elements(qv_) then begin
	   print,'Error!'
	   stop
	endif

	;interpolate into 41 layers
	lyr = 41
	pres=fltarr(lyr) & temp=fltarr(lyr) & qvx=fltarr(lyr) & htx=fltarr(lyr) & o3x=fltarr(lyr) & rhx=fltarr(lyr)
        htx[0:lyr-1]=interpol(height_,lyr)        
	pres[0:lyr-1]=interpol(lev_,lyr)
        temp[0:lyr-1]=interpol(temp_,lyr)
	qvx[0:lyr-1]=interpol(qv_,lyr)
        rhx[0:lyr-1]=interpol(relh_,lyr)
        o3x[0:lyr-1]=interpol(o3_,lyr)
	if htx[0] ge htx[1] or pres[1] ge pres[0] then print,"ERROR: ZERO THICKNESS IN FIRST LAYER"

	;use layer boundaries to get mid-layer average for temp, QV, and O3
	tempm=fltarr(lyr-1)
	o3m=fltarr(lyr-1)
	qvm=fltarr(lyr-1)
	for i=1,lyr-1 do begin
	    tempm[i-1]=(temp[i-1]+temp[i])/2
	    qvm[i-1]=(qvx[i-1]+qvx[i])/2
	    o3m[i-1]=(o3x[i-1]+o3x[i])/2
	endfor

	;reverse and convert to a strings
	h=strarr(lyr) & p=strarr(lyr) & t=strarr(lyr-1) & q=strarr(lyr-1) & rh=strarr(lyr) & ozone=strarr(lyr-1)
	for i=0,lyr-2 do begin
	   h[i]=string(htx[lyr-i-1],format='(f10.4)')
	   p[i]=string(pres[lyr-i-1],format='(f10.4)')
	   t[i]=string(tempm[lyr-i-2],format='(f10.4)')
	   q[i]=string(qvm[lyr-i-2],format='(e15.4)')
	   rh[i]=string(rhx[lyr-i-1],format='(f10.4)')
	   ozone[i]=string(o3m[lyr-i-2],format='(e15.4)')
	endfor
	h[i]=string(htx[lyr-i-1],format='(f10.4)')
	p[i]=string(pres[lyr-i-1],format='(f10.4)')
	rh[i]=string(rhx[lyr-i-1],format='(f10.4)')
	tskin=string(300.,format='(f10.1)')
	

	;///////////Write the file and call the code////////////////////
	openw,lun,'input_file.dat',/get_lun
	printf, lun, string(sza1,format='(f10.4)'),string(AOD,format='(f10.4)'), $
	string(sctyp,format='(i10)'),summer, string(aerotype,format='(i10)'), string(jday0,format='(I5)')
	for i=0,lyr-2 do begin
	  printf,lun,h[i],rh[i],p[i],t[i],q[i],ozone[i]
	endfor
	printf,lun,h[i],rh[i],p[i],tskin
	free_lun,lun

	spawn,'gfortran *.F90'

	spawn,'./a.out'

	spawn,'cp output.dat '+station+'_UTC_'+yymmdd+'_'+utc+'.data'
	spawn,'mv '+station+'_UTC_'+yymmdd+'_'+utc+'.data /home/jbrendecke/rad_stations/'+station+'/CKD/model_raw/'

     endif else begin

	openw,lun,station+'_UTC_'+yymmdd+'_'+utc+'.data' ,/get_lun
 	  printf,lun,'     Lev(mb)    SW Up(Wm2)     SW Dw(Wm2)'
	for i=0,nlev-1 do begin 
	  printf,lun,'   ---------        0.0000       0.0000'
	endfor
	printf,lun,'                Direct(Wm2) Diffuse(Wm2)'
       printf,lun,'   ---------        0.0000       0.0000'
	printf,lun,'Surface Downward Band Flux(Wm2): '
	printf,lun,'       Band1         Band2         Band3         Band4'
	printf,lun,'      0.0000        0.0000        0.0000        0.0000'
	printf,lun,' -----------------------------------------------------------'
	free_lun, lun
	spawn,'mv '+station+'_UTC_'+yymmdd+'_'+utc+'.data /home/jbrendecke/rad_stations/'+station+'/CKD/model_raw/'
     endelse
endfor
endfor
 

stop
END;=========================================================================
