#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sept 10 12:45:19 2024
#load imput for CCCma with the inclusion of clouds from CERES SYN1
#run loop based on specified cases of cloud from 2013-2023
# version 2 runs with version 10 of CCCma, different output
@author: jbrendecke
""" 
import numpy as np
import xarray as xr
import subprocess
import os
import pickle
import glob
import warnings
warnings.filterwarnings("ignore")
import sys
#==========================================================================================================
# get the doy
def set_jday(iyr, imo, iday):
#    ; Determine Julian Day from Calendar Day
#    ; INPUT:
#    ;; ------------------------------
#    ;; iyr = 4-digit year (INTEGER)
#    ;; imo = 2-digit month (INTEGER)
#    ;; iday = 2-digit day (INTEGER)
#    ;; 
#    ;; OUTPUT:
#    ;; ------------------------------
#    ;; iddd = Julian Day (0-365) (INTEGER)

    idymon = [[0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],  # non-leap year
              [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]]  # leap year

    leap = iyr - (iyr // 4) * 4

    if leap == 0:  # Leap Year
        irow = 1
    else:  # Non-leap Year
        irow = 0

    iddd = idymon[irow][imo-1] + iday
    return iddd
#==========================================================================================================
# calculate the SZA
def zenith(utc, lat, lon):
        # INPUTS: 
        #	utc	flt or fltarr: the time (UTC) as days since Jan 1 00:00 of the year
        #		e.g. the time 0.5 is noon of Jan 1st.
        #		(presently it also works for the following year, if the first year
        #		 had 365 days. I.e it works for the NOXAR times for 1995 and 1996 but
        #		 not yet for 1997
        #	lat	flt or fltarr: the latitude in deg (pos. for northern hemisphere)
        #	lon	flt or fltarr: the longitude in deg (pos. for eastern longitudes)
        # OUTPUTS
        #	the solar zenith angle in radian
    #dim=min(len(utc),len(lat),len(lon))
    #utc=utc[0:dim] if len(utc) > dim else utc
    #lat=lat[0:dim] if len(lat) > dim else lat
    #lon=lon[0:dim] if len(lon) > dim else lon
    
    # calculate the number of days since the 1.1. of the specific year
    daynum=np.floor(utc)+1
    #calculate the relative SUN-earth distance for the given day
    #resulting from the elliptic orbit of the earth
    eta=2.*np.pi*daynum/365.
  
    #calculate the solar declination for the given day
    # the declination varies due to the fact, that the earth rotation axis
    # is not perpendicular to the ecliptic plane
    delta=0.006918 \
    -0.399912*np.cos(eta)-0.006758*np.cos(2*eta)-0.002697*np.cos(3*eta) \
    +0.070257*np.sin(eta)+0.000907*np.sin(2*eta)+0.001480*np.sin(3*eta)
    
    #equation of time, used to compensate for the earth's elliptical orbit
    # around the sun and its axial tilt when calculating solar time
    # eqt is the correction in hours
    et=2.*np.pi*daynum/366.
    eqt=0.0072*np.cos(et)-0.0528*np.cos(2.*et)-0.0012*np.cos(3.*et) \
    -0.1229*np.sin(et)-0.1565*np.sin(2.*et)-0.0041*np.sin(3.*et)
    
    #calculate SZA
    dtr=np.pi/180.
    time=(utc+1-daynum)*24 #time in hours
    omega=(360/24)*(time+lon/15+eqt-12)*dtr
    sinh=np.sin(delta)*np.sin(lat*dtr)+np.cos(delta)*np.cos(lat*dtr)*np.cos(omega)
    solel=np.arcsin(sinh)
    sza=np.pi/2-solel
    
    return sza
#==========================================================================================================
#cloud profile using constant Nc, and Increasing LWC/IWC, and Rc/Ri
def cloud_profile_lw_liq_sgp(hgt,pres,Rc,tau ,CB_pres,CT_pres,CF_total):
    #based on LWC lineraly/adiabaticly increasing
    #values are computed at each hgt/pres level
    #Inputs:
        #hgt: profile of heights from surface to TOA [km]
        #pres: profile of pressure form surface to TOA [hPa]
        #Rc: single value, either liquid or ice particle size [um]
        #Rc_slp: slope of in change of Rc with height [um/km]
        #LWP: total LWP or IWP [g m-2]
        #CB_pres: Cloud Base pressure level [hPa]
        #CT_pres: Cloud Top pressure level [hPa]
        #CF_total: Total Cloud Fraction
    #Outputs:
        #LWC_profile: Profile of LWC or IWC 
        #Rc_profile: Profile of Rc or Ri
        #CF_profile: Profile of Cloud Fraction
        
    CF_profile = np.zeros(len(hgt))    
    Rc_profile = np.zeros(len(hgt))  
    LWC_profile = np.zeros(len(hgt))
    
    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))
    
    lwp = (tau * 5 * Rc*1e-6 *100**3) / 9

    if Rc < 5:
        Rc = 5
    if Rc > 20:
        Rc =20
    if (cldlevb == cldlevt): #single cloud layer
         cldthick = hgt[cldlevt+1]-hgt[cldlevb]    
         Rc_profile[cldlevb:cldlevb+2] = Rc
         CF_profile[cldlevb:cldlevb+2] = CF_total
         LWC_profile[cldlevb:cldlevb+2] = lwp / (cldthick*1000)
        
    else: #multiple cloud layers
        cldthick = (hgt[cldlevt]-hgt[cldlevb])*1000 #m
        #CLOUD FRACTION  //////////////////////////////////////////////////////       
        CF_profile[cldlevb:cldlevt+1] =CF_total
        #PARTICLE SIZE ////////////////////////////////////////////////////////
        Rc /= 1e6 #m
        rc_st = 5e-6 #m
        lwp /= 1000 #kg/m2
        rho = 1000 #kg/m3
        
        #change were MODIS re is placed(rather than midway)
        #rc_slp = (Rc - rc_st)/(cldthick *(1/3))
        #Rc = rc_slp*(cldthick *(1/2)) + rc_st
        #if Rc > 20/1e6:
        #    Rc = 20/1e6
        
        Nc = (3*(lwp/cldthick)) / (4 * np.pi  * rho * Rc**3)
        lwcm = (4/3) * np.pi * rho * Rc**3 * Nc 
        lwc_st = (4/3) * np.pi * rho * (rc_st)**3 * Nc
        lwc_slp = (lwcm - lwc_st)/(cldthick * .5)
        
        #lwc_top= lwc_slp*cldthick + lwc_st
        #rc_top = ((3 * lwc_top) /(4 * np.pi * Nc * rho))**(1/3)
        rc_slp = (Rc- rc_st)/(cldthick * 0.5)
        rc_top = (rc_slp*cldthick +rc_st)
        
        
        if rc_top*1e6 > 20:
            rc_top = 20e-6
            #lwc_top = (4/3) * np.pi * rho * (rc_top)**3 * Nc
            #lwc_slp = (lwc_top - lwcm)/(cldthick *.5)
            #lwc_st = (lwc_slp*(-cldthick *.5))+lwcm
            
            #rc_slp = (rc_top - Rc)/(cldthick * 0.5)
            #rc_st = Rc - (rc_slp*(cldthick * 0.5))

        #print(f'Thick: {cldthick:.3f}, Nc: {Nc/100**3:.3f}, LWP: {LWP*1000:.3f}, Rc:{Rc*1E6:.3f}')
        for i in range(cldlevb,cldlevt+1):
            LWC_profile[i] =  (hgt[i] - hgt[cldlevb])*1000 *lwc_slp + lwc_st
            #Rc_profile[i] = (((3 * LWC_profile[i]) / (4 * np.pi * Nc  * rho))**(1/3)) *1e6 #non-linear
            Rc_profile[i] = ((hgt[i] - hgt[cldlevb])*1000 *rc_slp + rc_st)*1e6 #linear
        
        Rc_profile[Rc_profile > 20] = 20
        
        #decrease top portions of cloud due to entrainment
        if cldlevt - cldlevb > 3:
            lthick2= cldthick * 0.75
            cldlev75 = np.argmin(abs((hgt[cldlevb] + lthick2/1000) - (hgt)))
            rc75 = Rc_profile[cldlev75]
            lwc75 = LWC_profile[cldlev75]
            
            rc_top = rc75 * 1.05 #decrease in size
            lwc_top = (4/3) * np.pi * rho * (rc_top/1e6)**3 * Nc
            
            lwctop= (LWC_profile[cldlevt]/lwc_top)


            rc_slpt = (rc_top - rc75)/((hgt[cldlevt] - hgt[cldlev75])*1000)
            lwc_slpt = (lwc_top - lwc75)/((hgt[cldlevt] - hgt[cldlev75])*1000)
            for i in range(cldlev75, cldlevt+1):
                Rc_profile[i] = ((hgt[i] - hgt[cldlev75])*1000 * rc_slpt + rc75)
                LWC_profile[i] = ((hgt[i] - hgt[cldlev75])*1000  *lwc_slpt + lwc75)
                
            Rc_profile[Rc_profile > 20] = 20

 
            lwp_new =0
            lthick = (hgt[cldlev75] - hgt[cldlevb])*1000
            integral = (lthick*(lwc_slp *lthick + (2 * LWC_profile[cldlevb])))/2
            lwp_new += integral
            lthick = (hgt[cldlevt]-hgt[cldlev75])*1000
            integral = (lthick*(lwc_slpt *lthick + (2 * LWC_profile[cldlev75])))/2 
            lwp_new += integral
            
            lwpdiff = lwp - lwp_new            
            for i in range(cldlevb,cldlevt+1):
                LWC_profile[i] += lwpdiff/cldthick
            
            
            lwp_check = 0
            lwc_slp2 = (LWC_profile[cldlev75] - LWC_profile[cldlevb])/((hgt[cldlev75] - hgt[cldlevb])*1000)
            lwc_slpt2 = (LWC_profile[cldlevt] - LWC_profile[cldlev75])/((hgt[cldlevt] - hgt[cldlev75])*1000)
            for i in range(cldlevb, cldlev75):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slp2 *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
            for i in range(cldlev75, cldlevt):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slpt2 *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
                
            
            if abs(lwp*1000 - lwp_check*1000) > 0.1:
                print('CLOUD PROFILE FUNCTION IS INCORRECT!!!')
                print(lwp_check*1000, lwp*1000)
        
        else:  
            lwctop=0
            lwp_check = 0
            for i in range(cldlevb, cldlevt):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slp *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
                
            if abs(lwp*1000 - lwp_check*1000) > 0.1:
                print('CLOUD PROFILE FUNCTION IS INCORRECT!!!')
                print(lwp_check*1000,lwp*1000)   

        LWC_profile *= 1000
        
        if max(Rc_profile) > 20.01:
            print('Error Over; Rc = ', Rc_profile[(Rc_profile > 20.01)])  
            
        xx = Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)]
        if len(xx) > 0:
            print('Error Under; Rc = ', Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)])  
    return LWC_profile, Rc_profile, CF_profile, lwctop
#==========================================================================================================
def cloud_profile_lw_liq_ena(hgt,pres,Rc,tau ,CB_pres,CT_pres,CF_total):
    #based on LWC lineraly/adiabaticly increasing
    #values are computed at each hgt/pres level
    #Inputs:
        #hgt: profile of heights from surface to TOA [km]
        #pres: profile of pressure form surface to TOA [hPa]
        #Rc: single value, either liquid or ice particle size [um]
        #Rc_slp: slope of in change of Rc with height [um/km]
        #LWP: total LWP or IWP [g m-2]
        #CB_pres: Cloud Base pressure level [hPa]
        #CT_pres: Cloud Top pressure level [hPa]
        #CF_total: Total Cloud Fraction
    #Outputs:
        #LWC_profile: Profile of LWC or IWC 
        #Rc_profile: Profile of Rc or Ri
        #CF_profile: Profile of Cloud Fraction
        
    CF_profile = np.zeros(len(hgt))    
    Rc_profile = np.zeros(len(hgt))  
    LWC_profile = np.zeros(len(hgt))
    
    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))
    
    lwp = (tau * 5 * Rc*1e-6 *100**3) / 9

    if Rc < 5:
        Rc = 5
    if Rc > 20:
        Rc =20
        
    if (cldlevb == cldlevt): #single cloud layer
         cldthick = hgt[cldlevt+1]-hgt[cldlevb]    
         Rc_profile[cldlevb:cldlevb+2] = Rc
         CF_profile[cldlevb:cldlevb+2] = CF_total
         LWC_profile[cldlevb:cldlevb+2] = lwp / (cldthick*1000)
        
    else: #multiple cloud layers
        cldthick = (hgt[cldlevt]-hgt[cldlevb])*1000 #m
        #CLOUD FRACTION  //////////////////////////////////////////////////////       
        CF_profile[cldlevb:cldlevt+1] =CF_total
        #PARTICLE SIZE ////////////////////////////////////////////////////////
        Rc /= 1e6 #m
        rc_st = Rc#7e-6 #m
        lwp /= 1000 #kg/m2
        rho = 1000 #kg/m3
        
        #change were MODIS re is placed(rather than midway)
        rc_slp = (Rc - rc_st)/(cldthick *(1/3))
        Rc = rc_slp*(cldthick *(1/2)) + rc_st
        if Rc > 20/1e6:
            Rc = 20/1e6
        
        Nc = (3*(lwp/cldthick)) / (4 * np.pi  * rho * Rc**3)
        lwcm = (4/3) * np.pi * rho * Rc**3 * Nc 
        lwc_st = (4/3) * np.pi * rho * (rc_st)**3 * Nc
        lwc_slp = (lwcm - lwc_st)/(cldthick * .5)
        
        #lwc_top= lwc_slp*cldthick + lwc_st
        #rc_top = ((3 * lwc_top) /(4 * np.pi * Nc * rho))**(1/3)
        rc_slp = (Rc- rc_st)/(cldthick * 0.5)
        rc_top = (rc_slp*cldthick +rc_st)
        
        
        if rc_top*1e6 > 20:
            rc_top = 20e-6
            #lwc_top = (4/3) * np.pi * rho * (rc_top)**3 * Nc
            #lwc_slp = (lwc_top - lwcm)/(cldthick *.5)
            #lwc_st = (lwc_slp*(-cldthick *.5))+lwcm
            
            #rc_slp = (rc_top - Rc)/(cldthick * 0.5)
            #rc_st = Rc - (rc_slp*(cldthick * 0.5))

        #print(f'Thick: {cldthick:.3f}, Nc: {Nc/100**3:.3f}, LWP: {LWP*1000:.3f}, Rc:{Rc*1E6:.3f}')
        for i in range(cldlevb,cldlevt+1):
            LWC_profile[i] =  (hgt[i] - hgt[cldlevb])*1000 *lwc_slp + lwc_st
            #Rc_profile[i] = (((3 * LWC_profile[i]) / (4 * np.pi * Nc  * rho))**(1/3)) *1e6 #non-linear
            Rc_profile[i] = ((hgt[i] - hgt[cldlevb])*1000 *rc_slp + rc_st)*1e6 #linear
        
        Rc_profile[Rc_profile > 20] = 20
        
        #decrease top portions of cloud due to entrainment
        if cldlevt - cldlevb > 3000:#!!!change back to <3
            lthick2= cldthick * 0.75
            cldlev75 = np.argmin(abs((hgt[cldlevb] + lthick2/1000) - (hgt)))
            rc75 = Rc_profile[cldlev75]
            lwc75 = LWC_profile[cldlev75]
            
            rc_top = rc75 * .85 #decrease in size
            lwc_top = (4/3) * np.pi * rho * (rc_top/1e6)**3 * Nc
            
            lwctop= (LWC_profile[cldlevt]/lwc_top)


            rc_slpt = (rc_top - rc75)/((hgt[cldlevt] - hgt[cldlev75])*1000)
            lwc_slpt = (lwc_top - lwc75)/((hgt[cldlevt] - hgt[cldlev75])*1000)
            for i in range(cldlev75, cldlevt+1):
                Rc_profile[i] = ((hgt[i] - hgt[cldlev75])*1000 * rc_slpt + rc75)
                LWC_profile[i] = ((hgt[i] - hgt[cldlev75])*1000  *lwc_slpt + lwc75)
                
            Rc_profile[Rc_profile > 20] = 20

            #add back lwc to so that lwp is constrained
            lwp_new =0
            lthick = (hgt[cldlev75] - hgt[cldlevb])*1000
            integral = (lthick*(lwc_slp *lthick + (2 * LWC_profile[cldlevb])))/2
            lwp_new += integral
            lthick = (hgt[cldlevt]-hgt[cldlev75])*1000
            integral = (lthick*(lwc_slpt *lthick + (2 * LWC_profile[cldlev75])))/2 
            lwp_new += integral
            
            lwpdiff = lwp - lwp_new            
            for i in range(cldlevb,cldlevt+1):
                LWC_profile[i] += lwpdiff/cldthick
            
            #double check new lwp matches satellite lwp
            lwp_check = 0
            lwc_slp2 = (LWC_profile[cldlev75] - LWC_profile[cldlevb])/((hgt[cldlev75] - hgt[cldlevb])*1000)
            lwc_slpt2 = (LWC_profile[cldlevt] - LWC_profile[cldlev75])/((hgt[cldlevt] - hgt[cldlev75])*1000)
            for i in range(cldlevb, cldlev75):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slp2 *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
            for i in range(cldlev75, cldlevt):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slpt2 *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
                
            #print(lwp*1000, lwp_check*1000, lwp_new*1000)
            
            if abs(lwp*1000 - lwp_check*1000) > 0.1:
                print('CLOUD PROFILE FUNCTION IS INCORRECT!!!')
                print(lwp_check*1000, lwp*1000)
        
        else:  
            lwctop=0
            lwp_check = 0
            for i in range(cldlevb, cldlevt):
                lthick = (hgt[i+1]-hgt[i])*1000
                integral = (lthick*(lwc_slp *lthick + (2 * LWC_profile[i])))/2 
                lwp_check += integral
                
            if abs(lwp*1000 - lwp_check*1000) > 0.1:
                print('CLOUD PROFILE FUNCTION IS INCORRECT!!!')
                print(lwp_check*1000,lwp*1000)   

        LWC_profile *= 1000
        
        if max(Rc_profile) > 20.01:
            print('Error Over; Rc = ', Rc_profile[(Rc_profile > 20.01)])  
            
        xx = Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)]
        if len(xx) > 0:
            print('Error Under; Rc = ', Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)])  
    return LWC_profile, Rc_profile, CF_profile
#==========================================================================================================
def cloud_profile_constant(hgt, pres, Rc, tau, CB_pres, CT_pres, CF, ICE):
    #Coputes constant RC and LWC with height
    #Inputs:
        #hgt: profile of heights from surface to TOA [km]
        #pres: profile of pressure form surface to TOA [hPa]
        #Rc: single value, either liquid or ice particle size [um]
        #LWP: total LWP or IWP [g m-2]
        #CB_pres: Cloud Base pressure level [hPa]
        #CT_pres: Cloud Top pressure level [hPa]
    #Outputs:
        #LWC_profile: Profile of LWC or IWC 
        #Rc_profile: Profile of Rc or Ri
        
    Re_profile = np.full(len(hgt), np.nan)  
    WC_profile = np.full(len(hgt), np.nan)
    CF_profile = np.full(len(hgt), np.nan) 

    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))
    
    if ICE:
        #de = Rc*2
        #wp = tau*((0.259*de) + (0.819e-3*de**2) - (0.880e-6*de**3))
        wp = (tau *2.1 * Rc*1e-6 * 93**3) /3
        
    else:
        wp = (tau * 5 * Rc*1e-6 *100**3) / 9   
    
    if (cldlevb == cldlevt):
        CF_profile[cldlevb:cldlevb+2] =CF
        cldthick = (hgt[cldlevt+1]-hgt[cldlevb])*1000 
        Re_profile[cldlevb:cldlevt+2] = Rc
        WC_profile[cldlevb:cldlevt+2] = wp / cldthick
    else:
        CF_profile[cldlevb:cldlevt+1] = CF
        cldthick = (hgt[cldlevt] - hgt[cldlevb])*1000
        Re_profile[cldlevb:cldlevt+1] = Rc
        WC_profile[cldlevb:cldlevt+1] = wp / cldthick
            
    check = 0
    for i in range(cldlevb, cldlevt):
        lthick = (hgt[i+1]-hgt[i])*1000
        check += lthick *WC_profile[i]
        
    if abs(check -wp) > 0.1:
        print('CLOUD PROFILE NOT CORRECT')
        #print(check,wp)
        
    return WC_profile, Re_profile, CF_profile
#==========================================================================================================
def ice_cloud_profile(hgt,pres, rc, tau, CB_pres, CT_pres, CF):
    #based on LWC lineraly/adiabaticly increasing
    #values are computed at each hgt/pres level
    #Inputs:
        #hgt: profile of heights from surface to TOA [km]
        #pres: profile of pressure form surface to TOA [hPa]
        #Rc: single value, either liquid or ice particle size [um]
        #Rc_slp: slope of in change of Rc with height [um/km]
        #LWP: total LWP or IWP [g m-2]
        #CB_pres: Cloud Base pressure level [hPa]
        #CT_pres: Cloud Top pressure level [hPa]
        #CF_total: Total Cloud Fraction
    #Outputs:
        #LWC_profile: Profile of LWC or IWC 
        #Rc_profile: Profile of Rc or Ri
        #CF_profile: Profile of Cloud Fraction
        
    CF_profile = np.zeros(len(hgt))    
    Rc_profile = np.zeros(len(hgt))  
    IWC_profile = np.zeros(len(hgt))
    
    cldlevb = np.argmin(abs(CB_pres-pres))
    cldlevt = np.argmin(abs(CT_pres-pres))

    iwp = (tau *2* rc*1e-6 *96**3) / 3
    #ver2 rho=94.5
    #ver3 rho=95
    #ver4 rho=94
    #ver5 rho=96
    #ver6 rho=96 but 41 layers
    #ver6b rho=96, 41 layers, both hgh & midh layers
    #ver7 rho=96, 41 layers, select cases, and rc_st = *1.5
    #ver8 rho=96, 41 layers, select cases, and rc_st = *2
    #ver9 rho=96, 41 layers, select cases, and rc_st = *1/1.5/2
    #ver10 rho=96, 41 layers, select cases, and rc_st = *1.25
    #ver11 rho=96, 41 layers, select cases, and rc_st = *1.1

    
    
    
    if rc < 5:
        rc = 5
    if rc > 60:
        rc =60
        
    if (cldlevb == cldlevt): #single cloud layer
         cldthick = hgt[cldlevt+1]-hgt[cldlevb]    
         Rc_profile[cldlevb:cldlevb+2] = rc
         CF_profile[cldlevb:cldlevb+2] = CF
         IWC_profile[cldlevb:cldlevb+2] = iwp / (cldthick*1000)
        
    else: #multiple cloud layers
        cldthick = (hgt[cldlevt]-hgt[cldlevb])*1000 #m
        #CLOUD FRACTION  //////////////////////////////////////////////////////       
        CF_profile[cldlevb:cldlevt+1] =CF
        #PARTICLE SIZE ////////////////////////////////////////////////////////
        if (iwp >= 100) & (iwp < 300):
            scale = 1.5
        elif (iwp > 300):
            scale = 2.0
        else:
            scale = 1.0
        
        scale=1.0
            
        rc /= 1e6 #m
        rc_st = rc * scale#1.5 #m
        iwp /= 1000 #kg/m2
        rho = 930 #kg/m3
        
        rc_slp = (rc - rc_st)/(cldthick *(1/2))
        Nc = (3*(iwp/cldthick)) / (4 * np.pi  * rho * rc**3)
        iwcm = (4/3) * np.pi * rho * rc**3 * Nc
        iwc_st = (4/3) * np.pi * rho * (rc_st)**3 * Nc
        iwc_slp = (iwcm - iwc_st)/(cldthick * .5)
        
        for i in range(cldlevb,cldlevt+1):
            IWC_profile[i] =  (hgt[i] - hgt[cldlevb])*1000 *iwc_slp + iwc_st
            #Rc_profile[i] = (((3 * LWC_profile[i]) / (4 * np.pi * Nc  * rho))**(1/3)) *1e6 #non-linear
            Rc_profile[i] = ((hgt[i] - hgt[cldlevb])*1000 *rc_slp + rc_st)*1e6 #linear
        
        #Rc_profile[Rc_profile > 60] = 60
        
        
        #double check new lwp matches satellite lwp
        iwp_check = 0
        for i in range(cldlevb, cldlevt):
            lthick = (hgt[i+1]-hgt[i])*1000
            integral = (lthick*(iwc_slp *lthick + (2 * IWC_profile[i])))/2 
            iwp_check += integral
            
                 
        if abs(iwp*1000 - iwp_check*1000) > 0.1:
            print('CLOUD PROFILE FUNCTION IS INCORRECT!!!')
            print(iwp_check*1000, iwp*1000)
    
        IWC_profile *= 1000
        
        if max(Rc_profile) > 60.01:
            print('Error Over; Rc = ', Rc_profile[(Rc_profile > 60.01)])  
            
        xx = Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)]
        if len(xx) > 0:
            print('Error Under; Rc = ', Rc_profile[(Rc_profile > 0) & (Rc_profile < 5)])  
    return IWC_profile, Rc_profile, CF_profile
#==========================================================================================================


stations =['SGP']
# file_ = '/home/jbrendecke/Ice_cloud_overcast_cases.pkl'
# with open(file_, 'rb') as f:
#     fdata = pickle.load(f)
# times = fdata['time']
# stations_ = fdata['station']
# stations_ = np.insert(stations_,45,'RUN')
# times = np.insert(times,45,np.datetime64('2020-02-22T02:22'))


xo=0
for station in stations:
    SAVE = True

    # ind_ = np.where(station == stations_)[0]
    # times_ = times[ind_]
    
    
    if station == 'SGP':
        lat_stat=36.5; lon_stat=-97.5
        sctyp_main=52; aerotype=2
    if station == 'ENA':
        #lat_stat =39.5; lon_stat=-28.5
        lat_stat =39.0916; lon_stat=-28.0257
        #lat_stat =39.5916; lon_stat=-28.5257
        sctyp_main=57; aerotype =1 #57
    
    
    if lon_stat < 0:
        lon_stat2 = lon_stat+360
    else:
        lon_stat2=lon_stat
    
    #/////////////LOAD SYN1 Cloud and AOD data///////////////////////////////
    path =f'/home/jbrendecke/rad_stations/{station}/ceres/SYN1_cloud/'
    filesyn1 = glob.glob(path+'*_Ed4.2_*.nc')
    nfiles=len(filesyn1)
    
    ds = xr.open_dataset(filesyn1[0])
    ds = ds.sel(lat = lat_stat, lon = lon_stat2, method='nearest')
    #ds = ds.where(ds['cldicerad_low_1h'].notnull, drop =True)
    time_syn=np.array(ds.time, dtype='datetime64[m]')
    lon_syn = np.array(ds.lon)
    lat_syn = np.array(ds.lat)
    cf_tt_syn = np.squeeze(ds['cldarea_total_1h'].data)/100.
    cf_hg_syn = np.squeeze(ds['cldarea_high_1h'].data)/100.
    cf_mh_syn = np.squeeze(ds['cldarea_mid_high_1h'].data)/100.
    cf_ml_syn = np.squeeze(ds['cldarea_mid_low_1h'].data)/100.
    cf_lw_syn = np.squeeze(ds['cldarea_low_1h'].data)/100.
    cod_tt_syn = np.squeeze(ds['cldtau_total_1h'].data)
    cld_bas_hg_syn = np.squeeze(ds['cldpress_base_high_1h'].data)
    cld_bas_mh_syn = np.squeeze(ds['cldpress_base_mid_high_1h'].data)
    cld_bas_ml_syn = np.squeeze(ds['cldpress_base_mid_low_1h'].data)
    cld_bas_lw_syn = np.squeeze(ds['cldpress_base_low_1h'].data)
    cld_bas_tt_syn = np.squeeze(ds['cldpress_base_total_1h'].data)
    cld_top_hg_syn = np.squeeze(ds['cldpress_top_high_1h'].data)
    cld_top_mh_syn = np.squeeze(ds['cldpress_top_mid_high_1h'].data)
    cld_top_ml_syn = np.squeeze(ds['cldpress_top_mid_low_1h'].data)
    cld_top_lw_syn = np.squeeze(ds['cldpress_top_low_1h'].data)
    cld_top_tt_syn = np.squeeze(ds['cldpress_top_total_1h'].data)
    rc_tt_syn = np.squeeze(ds['cldwatrad_total_1h'].data)
    rc_hg_syn = np.squeeze(ds['cldwatrad_high_1h'].data)
    rc_mh_syn = np.squeeze(ds['cldwatrad_mid_high_1h'].data)
    rc_ml_syn = np.squeeze(ds['cldwatrad_mid_low_1h'].data)
    rc_lw_syn = np.squeeze(ds['cldwatrad_low_1h'].data)
    ri_tt_syn = np.squeeze(ds['cldicerad_total_1h'].data)
    ri_hg_syn = np.squeeze(ds['cldicerad_high_1h'].data)
    ri_mh_syn = np.squeeze(ds['cldicerad_mid_high_1h'].data)
    ri_ml_syn = np.squeeze(ds['cldicerad_mid_low_1h'].data)
    ri_lw_syn = np.squeeze(ds['cldicerad_low_1h'].data)
    phase_tt_syn = np.squeeze(ds['cldphase_total_1h'].data)
    phase_hg_syn = np.squeeze(ds['cldphase_high_1h'].data)
    phase_mh_syn = np.squeeze(ds['cldphase_mid_high_1h'].data)
    phase_ml_syn = np.squeeze(ds['cldphase_mid_low_1h'].data)
    phase_lw_syn = np.squeeze(ds['cldphase_low_1h'].data)
    aod = np.squeeze(ds['ini_aod55_1h'].data)
    snowid =np.squeeze(ds['aux_snow_1h'].data)/100.
    #swd_ceres = np.squeeze(ds['adj_atmos_sw_down_all_surface_1h'].data)
    swd_ceres = np.squeeze(ds['ini_sfc_sw_down_all_1h'].data)
    #dir_ceres = np.squeeze(ds['adj_sfc_sw_direct_all_1h'].data)
    #dif_ceres = np.squeeze(ds['adj_sfc_sw_diff_all_1h'].data)
    swtoa_obs = np.squeeze(ds['toa_sw_all_1h'].data)
    swtoa_ini = np.squeeze(ds['ini_toa_sw_all_1h'].data)
    
    
    #time_syn = time_syn - int(time_syn[0]) + jday_1stofmonth
    ntime=len(time_syn)
    year = ds.time.dt.year.data
    jday = (ds.time.dt.dayofyear + (.5+ds.time.dt.hour)/24.).data
    
    month= ds.time.dt.month.data
    day= ds.time.dt.day.data
    hr = ds.time.dt.hour.data
    minute = ds.time.dt.minute.data
    minute[:]= 30
    
    #%%
    swd_toa_vis =[]; swd_toa_nir=[]
    swu_toa_all_vis =[]; swu_toa_all_nir =[]; swu_toa_clr_vis =[]; swu_toa_clr_nir =[];
    swu_20km_all_vis =[]; swu_20km_all_nir =[]; swu_20km_clr_vis =[]; swu_20km_clr_nir =[];
    swd_sfc_all_vis =[]; swd_sfc_all_nir =[]; swd_sfc_clr_vis =[]; swd_sfc_clr_nir =[];
    swu_sfc_all_vis =[]; swu_sfc_all_nir =[]; swu_sfc_clr_vis =[]; swu_sfc_clr_nir =[];
    direct_sfc_all =[]; diffuse_sfc_all =[]
    direct_sfc_clr =[]; diffuse_sfc_clr =[]
    direct_all_ceres =[]; diffuse_all_ceres=[]
    cod_ceres=[]; aod_out=[]; sza_out=[]; swd_sur_ceres=[]; ice_re=[]; liq_re=[];
    cld_top=[]; cld_base=[]; cb_h=[]; ct_h=[]; phase=[];
    swtoa_obs_ceres=[]; swtoa_ini_ceres=[]
    n=[]; nclevs=[]; lwcdiff=[]
    nnn=[]
    
    for ihr in range(ntime):
        #only do low-level liquid cloud with CF ~100%
        #if (cf_lw_syn[ihr] > 0.99) & (np.isnan(ri_lw_syn[ihr])) & (np.isnan(rc_hg_syn[ihr])) \
        #    & (np.isnan(rc_mh_syn[ihr])) & (np.isnan(rc_ml_syn[ihr])) & (~np.isnan(rc_lw_syn[ihr])) \
        #    & (np.isnan(ri_tt_syn[ihr])) & (rc_tt_syn[ihr] < 21) & (year[ihr] >= 2013) & (snowid[ihr] < 10):
        #if (cf_lw_syna[ihr] > 0.99) & (np.isnan(ri_lw_syna[ihr])) & (np.isnan(rc_hg_syna[ihr])) \
        #    & (np.isnan(rc_mh_syna[ihr])) & (np.isnan(rc_ml_syna[ihr])) & (~np.isnan(rc_lw_syna[ihr])) \
        #    & (np.isnan(ri_tt_syna[ihr])) & (rc_tt_syna[ihr] < 21) & (year[ihr] >= 2013) & (snowida[ihr] < 10):
         
        #only do ice-clouds and CF ~100%
        if  ((ri_tt_syn[ihr] > 0) & (np.isnan(rc_tt_syn[ihr])) & (cf_lw_syn[ihr] == 0.) & \
            (cf_ml_syn[ihr] == 0) & (cf_mh_syn[ihr] == 0) & (cf_hg_syn[ihr] > .99) & (snowid[ihr] < 0.2 )) \
            or \
            ((ri_tt_syn[ihr] > 0) & (np.isnan(rc_tt_syn[ihr])) & (cf_lw_syn[ihr] == 0.) & \
            (cf_ml_syn[ihr] == 0) & (cf_mh_syn[ihr] > .99) & (cf_hg_syn[ihr] == 0) & (snowid[ihr] < 0.2 )): #& \

        
            sza = zenith(jday[ihr]-1, lat_stat,lon_stat)*(180/np.pi)
            if sza < 70:
                xo +=1
                nnn.append(1)
                month_s = f'{month[ihr]:02}'
                day_s = f'{day[ihr]:02}'
                hr_s = f'{hr[ihr]:02}'
                min_s = f'{minute[ihr]:02}'
                
                utc = hr_s + min_s +'00'
                yymmdd = str(year[ihr])+month_s+day_s
                jday_int = int(jday[ihr])
                
                if (jday_int>75) & (jday_int<265):
                    summer=f'{1:9d}'
                else:
                    summer=f'{0:9d}'
    
                n.append(ihr)
                #time.sleep(0.1)
                #//////AOD////////
                AOD = aod[ihr]
                aod_out.append(AOD)
                
                #/////SNOW on surface/////
                if snowid[ihr] > 90:
                    sctyp = 55
                else:
                    sctyp = sctyp_main
                
                cod_ceres.append(cod_tt_syn[ihr])
                sza_out.append(sza)
                swd_sur_ceres.append(swd_ceres[ihr])
                #direct_all_ceres.append(dir_ceres[ihr])
                #diffuse_all_ceres.append(dif_ceres[ihr])
                swtoa_obs_ceres.append(swtoa_obs[ihr])
                swtoa_ini_ceres.append(swtoa_ini[ihr])
               
                
                #///////////////LOAD MERRA2 DATA/////////////////////////////////
                pathm=f'/mothership/data/MERRA2/MERRA2_Month_Data/{yymmdd[:6]}/'    
                dd = yymmdd[-2:]
                merra2 = glob.glob(f'{pathm}*'+dd+'.SUB.nc')[0]
                
                dsm = xr.open_dataset(merra2)
                dsm['time'] = dsm['time'] + np.timedelta64(90, 'm')
                dsm = dsm.sel(lat = lat_stat, lon = lon_stat2, time=time_syn[ihr], method='nearest')
                indx1 = np.where(~np.isnan(dsm['H']))[0]
                indx2 = np.where(~np.isnan(dsm['T']))[0]
                
                height = np.squeeze(dsm['H'].data)/1000.
                lev = np.squeeze(dsm['lev'].data)
                height = height[indx1]
                lev = lev[indx1]
                dsm = dsm.where(dsm['T'].notnull, drop =True)
                timeM2 = np.squeeze(dsm['time'].data)
                QV = np.squeeze(dsm['QV'].data)
                temp = np.squeeze(dsm['T'].data)
                O3 = np.squeeze(dsm['O3'].data)
                relh = np.squeeze(dsm['RH'].data)*100
                ps = np.squeeze(dsm['PS'].data[0])/100
                slp = np.squeeze(dsm['SLP'].data[0])/100
                elv = np.squeeze(dsm['PHIS'].data[0])/(1000*9.8)
                nlev = len(lev)
                
    
                #/////////////Write out MERRA2 Profile//////
                #Pressure, Hieght, RH should have 41 levels
                #Temperature, QV, O3 should have 40 layers
                
                if station == 'SYO':
                    elv = 0.03
                if elv < .01 :
                    elv = 0.
    
                #certain cases the height/pres has more values than other variables
                if len(indx1) != len(indx2):
                    delh = height[1] - height[0]
                    temp0 = temp[0]+delh*6.5
                    relh0 = relh[0]
                    O30 = O3[0]
                    QV0 = QV[0]
                    
                    temp = np.concatenate(([temp0],temp))
                    relh = np.concatenate(([relh0],relh))
                    O3 = np.concatenate(([O30],O3))
                    QV = np.concatenate(([QV0],QV))
                
                #make sure all times steps start at the same height
                if ps > lev[0]:
                    delh = height[0] - elv
                    ht0 = elv
                    lev0 = ps
                    if station == 'SYO':
                        lev0 = slp *((288-6.5*elv)/288)**5.2561
                    temp0 = temp[0]+delh*6.5
                    relh0 = relh[0]
                    O30 = O3[0]
                    QV0 = QV[0]
                    
                    height = np.concatenate(([ht0],height))
                    lev = np.concatenate(([lev0],lev))
                    temp = np.concatenate(([temp0],temp))
                    relh = np.concatenate(([relh0],relh))
                    O3 = np.concatenate(([O30],O3))
                    QV = np.concatenate(([QV0],QV))
                else:
                    delh = elv - height[0]
                    height[0] = elv
                    lev[0] = ps
                    temp[0] = temp[0] - delh*6.5
                   
                if len(height) != len(lev) | len(height) != len(temp) | \
                len(height) != len(O3) | len(height) != len(relh) | len(height) != len(QV):
                    print('ERROR: Number of levels differ for variables')
                    continue
                
                
                lyr = 41 #41
                presx = np.zeros(lyr); tempx=np.zeros(lyr); qvx=np.zeros(lyr)
                htx=np.zeros(lyr); o3x=np.zeros(lyr); rhx=np.zeros(lyr)
                
                #add layer below 500mb increase of 20 @ sea level
                ilev = -25 #500mb lev
                lyr_b500 = lyr - len(lev[ilev:])# new num of layers below 500mb
                
                x = np.linspace(0, 1, len(lev[:ilev]))  # Example x-coordinates, adjust if needed.
                new_x = np.linspace(0, 1, lyr_b500)  # New x-coordinates for interpolation
                
                presx[:lyr_b500] = np.interp(new_x, x, lev[:ilev])  # Assuming NP.INTERP
                htx[:lyr_b500] = np.interp(new_x, x, height[:ilev])
                tempx[:lyr_b500] = np.interp(new_x, x, temp[:ilev])
                qvx[:lyr_b500] = np.interp(new_x, x, QV[:ilev])
                rhx[:lyr_b500] = np.interp(new_x, x, relh[:ilev])
                o3x[:lyr_b500] = np.interp(new_x, x, O3[:ilev])
                
                presx[lyr_b500:] = lev[ilev:]
                htx[lyr_b500:] = height[ilev:]
                tempx[lyr_b500:] = temp[ilev:]
                qvx[lyr_b500:] = QV[ilev:]            
                rhx[lyr_b500:] = relh[ilev:]
                o3x[lyr_b500:] = O3[ilev:]
                
                '''
                #interpolat into 41 levels
                x = np.linspace(0, 1, len(lev))  # Example x-coordinates, adjust if needed.
                new_x = np.linspace(0, 1, lyr)  # New x-coordinates for interpolation
                
                presx[:lyr] = np.interp(new_x, x, lev)  # Assuming NP.INTERP
                htx[:lyr] = np.interp(new_x, x, height)
                tempx[:lyr] = np.interp(new_x, x, temp)
                qvx[:lyr] = np.interp(new_x, x, QV)
                rhx[:lyr] = np.interp(new_x, x, relh)
                o3x[:lyr] = np.interp(new_x, x, O3)
                '''
                #////////////////////////CLOUD PROPERTIES//////////////////////////////        
                #determine cloud properties using MERRA2 Hieghts and CERES SYN Properties
                rc = np.zeros(lyr)
                lwc = np.zeros(lyr)
                ri = np.zeros(lyr)
                iwc = np.zeros(lyr)
                cld_tot =  cf_tt_syn[ihr]
                if (cf_hg_syn[ihr] >.99):
                    cb_h.append(np.interp(cld_bas_hg_syn[ihr], np.flip(presx), np.flip(htx)))
                    ct_h.append(np.interp(cld_top_hg_syn[ihr], np.flip(presx), np.flip(htx)))
                    liq_re.append(rc_tt_syn[ihr])
                    ice_re.append(ri_hg_syn[ihr])
                    cld_top.append(cld_top_hg_syn[ihr])
                    cld_base.append(cld_bas_hg_syn[ihr])
                    phase.append(phase_hg_syn[ihr])
                elif (cf_mh_syn[ihr] >.99):
                    cb_h.append(np.interp(cld_bas_mh_syn[ihr], np.flip(presx), np.flip(htx)))
                    ct_h.append(np.interp(cld_top_mh_syn[ihr], np.flip(presx), np.flip(htx)))
                    liq_re.append(rc_tt_syn[ihr])
                    ice_re.append(ri_mh_syn[ihr])
                    cld_top.append(cld_top_mh_syn[ihr])
                    cld_base.append(cld_bas_mh_syn[ihr])
                    phase.append(phase_mh_syn[ihr])
                    
                
                if (cld_tot > 0.01):
                    cldlevb = np.argmin(abs(cld_bas_tt_syn[ihr]-presx))
                    cldlevt = np.argmin(abs(cld_top_tt_syn[ihr]-presx))
                nclevs.append(cldlevt-cldlevb)
                
                if ~np.isnan(rc_tt_syn[ihr]): 
                    ICE = False
                    lwc, rc, cld_cf = cloud_profile_lw_liq_ena(htx, presx, rc_tt_syn[ihr], cod_tt_syn[ihr], 
                                               cld_bas_tt_syn[ihr], cld_top_tt_syn[ihr], cf_tt_syn[ihr])
                    
                    #lwc, rc, cld_cf = cloud_profile_constant(htx, presx, rc_tt_syn[ihr], cod_tt_syn[ihr],
                    #                            cld_bas_tt_syn[ihr], cld_top_tt_syn[ihr], cf_tt_syn[ihr])
                if ~np.isnan(ri_tt_syn[ihr]):
                    ICE = True
                    iwc, ri, cld_cf = ice_cloud_profile(htx, presx, ri_tt_syn[ihr], cod_tt_syn[ihr], 
                                            cld_bas_tt_syn[ihr], cld_top_tt_syn[ihr], cf_tt_syn[ihr])
                
                #change MEERA2 profile based on MODIS Cloud boundary
                A =2.53E11
                B =5420
                ee = 0.622
                for i in range(1,lyr):
                    if (cld_cf[i] > 0.5):
                        rhx[i] = 100.0
                        #qvx1=(ee * A*np.exp(-B/tempx[i]))/(presx[i]*100) #specific humidity
                        #qvx[i] = qvx1/(1-qvx1) #mixing ratio
                          
                #user layer boundaries to get mid-layer avearage for temp, QV, and O3, Rc/Ri, & CF
                tempm = np.zeros(lyr-1)
                qvm = np.zeros(lyr-1)
                o3m = np.zeros(lyr-1)
                rcm = np.zeros(lyr-1)
                rim = np.zeros(lyr-1)
                lwcm = np.zeros(lyr-1)
                iwcm = np.zeros(lyr-1)
                cfm = np.zeros(lyr-1)
                htm = np.zeros(lyr-1)
                for i in range(1,lyr):
                    tempm[i-1] = (tempx[i-1]+tempx[i])/2
                    qvm[i-1] = (qvx[i-1]+qvx[i])/2
                    o3m[i-1] = (o3x[i-1]+o3x[i])/2
                    htm[i-1] = (htx[i-1]+htx[i])/2
                if cld_tot > 0.01:
                    if cldlevb==cldlevt:
                        cldlevt +=1
                    for i in range(cldlevb,cldlevt): 
                        lwcm[i] = (lwc[i+1]+lwc[i])/2
                        iwcm[i] = (iwc[i+1]+iwc[i])/2
                        rcm[i] = (rc[i+1]+rc[i])/2
                        rim[i] = (ri[i+1]+ri[i])/2
                        cfm[i] = (cld_cf[i+1]+cld_cf[i])/2
     
                
                #flip data to print on input file
                presx = np.flip(presx)
                htx = np.flip(htx)
                htm = np.flip(htm)
                rhx = np.flip(rhx)
                tempm = np.flip(tempm)
                qvm = np.flip(qvm)
                o3m = np.flip(o3m)
                rcm = np.flip(rcm)
                rim = np.flip(rim)
                lwcm = np.flip(lwcm)
                iwcm = np.flip(iwcm)
                cfm = np.flip(cfm)
                       
                
                #write data into the input file
                ckdpath= '/home/jbrendecke/CKD/ver11/'
                os.chdir(ckdpath)
                
                with open('UA_main.F90', 'r') as file:
                    dmain = file.readlines()
                dmain[4] = f'   parameter (lay = {lyr-1}, lev = lay + 1, ilg = 1, il1 = 1, il2 = 1, modlay=34, nxloc = 100)\n'
                with open('UA_main.F90', 'w') as file:
                    file.writelines(dmain)
                    
                with open('input_file.dat', 'w') as file:
                    file.write(f"{sza:10.4f} {AOD:9.4f} {sctyp:9d} {summer} {aerotype:9d} {jday_int:5d} {cld_tot:7.5f}\n")
                    for ilyr in range(lyr-1):
                        file.write(f"{htx[ilyr]:10.4f} {rhx[ilyr]:9.4f} {presx[ilyr]:9.4f} {tempm[ilyr]:9.4f} {qvm[ilyr]:14.4e} {o3m[ilyr]:14.4e}"
                                   f"{lwcm[ilyr]:14.4f} {iwcm[ilyr]:14.4f} {rcm[ilyr]:14.4f} {rim[ilyr]:14.4f} {cfm[ilyr]:12.5f}\n")
                    file.write(f"{htx[ilyr+1]:10.4f} {rhx[ilyr+1]:9.4f} {presx[ilyr+1]:9.4f}\n")           
            
                subprocess.run(['gfortran *.F90'], shell=True)
                subprocess.run(['./a.out'], shell=True)
                
                #subprocess.run(['cp output.dat '+station+'_UTC_'+yymmdd+'_'+utc+'.data'], shell=True)
                #subprocess.run(['mv '+station+'_UTC_'+yymmdd+'_'+utc+'.data /home/jbrendecke/rad_stations/'+station+'/CKD/model_raw_py/'], shell=True)
                
                #//////////////READ OUTPUT FILE////////////////////    
                output_file = ckdpath+'output.dat'
                nlines=49
                ind20km = np.argmin(abs(20 - htx))
                with open(output_file, 'r') as ofile:
                    lines = ofile.readlines()
                    swd_toa_vis.append(float(lines[1].split()[3]))
                    swd_toa_nir.append(float(lines[1].split()[4]))
                    swu_toa_all_vis.append(float(lines[1].split()[1]))
                    swu_toa_all_nir.append(float(lines[1].split()[2]))
                    swu_toa_clr_vis.append(float(lines[1].split()[5]))
                    swu_toa_clr_nir.append(float(lines[1].split()[6]))
                    swu_20km_all_vis.append(float(lines[ind20km].split()[1]))
                    swu_20km_all_nir.append(float(lines[ind20km].split()[2]))
                    swu_20km_clr_vis.append(float(lines[ind20km].split()[5]))
                    swu_20km_clr_nir.append(float(lines[ind20km].split()[6]))
                    swd_sfc_all_vis.append(float(lines[41].split()[3]))
                    swd_sfc_all_nir.append(float(lines[41].split()[4]))
                    swd_sfc_clr_vis.append(float(lines[41].split()[7]))
                    swd_sfc_clr_nir.append(float(lines[41].split()[8]))
                    swu_sfc_all_vis.append(float(lines[41].split()[1]))
                    swu_sfc_all_nir.append(float(lines[41].split()[2]))
                    swu_sfc_clr_vis.append(float(lines[41].split()[5]))
                    swu_sfc_clr_nir.append(float(lines[41].split()[6]))
                    direct_sfc_all.append(float(lines[(nlines-5)].split()[0]))
                    diffuse_sfc_all.append(float(lines[(nlines-5)].split()[1]))
                    direct_sfc_clr.append(float(lines[(nlines-2)].split()[0]))
                    diffuse_sfc_clr.append(float(lines[(nlines-2)].split()[1]))
             
    #%%
    lwcdiff = np.array(lwcdiff)
    lwcdiff =lwcdiff[lwcdiff != 0]
    
    #//////////SAVE OUTPUTS///////////////
    swd_toa = np.array(swd_toa_vis) + np.array(swd_toa_nir)
    swu_toa_clr = np.array(swu_toa_clr_vis) + np.array(swu_toa_clr_nir)
    swu_toa_all = np.array(swu_toa_all_vis) + np.array(swu_toa_all_nir)
    swu_20km_all = np.array(swu_20km_all_vis) + np.array(swu_20km_all_nir)
    swd_sur_clr = np.array(swd_sfc_clr_vis) + np.array(swd_sfc_clr_nir)
    swd_sur_all = np.array(swd_sfc_all_vis) + np.array(swd_sfc_all_nir)
    direct_clr = np.array(direct_sfc_clr)
    diffuse_clr = np.array(diffuse_sfc_clr)
    direct_all = np.array(direct_sfc_all)
    diffuse_all = np.array(diffuse_sfc_all)
    
    cod_ceres = np.array(cod_ceres)
    liq_re = np.array(liq_re)
    ice_re = np.array(ice_re)
    cld_base = np.array(cld_base)
    cld_top = np.array(cld_top)
    cb_h = np.array(cb_h)
    ct_h = np.array(ct_h)
    phase = np.array(phase)
    swd_sur_ceres = np.array(swd_sur_ceres)
    swtoa_obs_ceres = np.array(swtoa_obs_ceres)
    swtoa_ini_ceres = np.array(swtoa_ini_ceres)
    direct_all_ceres =np.array(direct_all_ceres)
    diffuse_all_ceres =np.array(diffuse_all_ceres)
    aod_out = np.array(aod_out)
    sza_out = np.array(sza_out)
    nclevs = np.array(nclevs)
    minute = minute[n]
    hour = hr[n]
    day =day[n]
    month = month[n]
    year = year[n]
    julday = jday[n]
    
    path = f'/home/jbrendecke/rad_stations/{station}/CKD/cloud_syn1/'
    fdata = {'JulianDay':julday, 'Year':year, 'Month':month, 'Day':day, 'Hour':hour, 'Minute':minute,
                 'AOD_ceres':aod_out, 'COD_ceres':cod_ceres, 'nclevs':nclevs,
                 'Re_liq':liq_re, 'Re_ice':ice_re, 'Cloud_bas':cld_base, 'Cloud_top':cld_top, 
                 'Cloud_base_h':cb_h, 'Cloud_top_h':ct_h, 'Phase':phase, 'SZA':sza_out,
                 'SWD_TOA': swd_toa, 'SWU_TOA_all':swu_toa_all, 'SWU_TOA_clr':swu_toa_clr, 'SWU_20km_all':swu_20km_all,
                 'SWU_TOA_Obs_ceres':swtoa_obs_ceres, 'SWU_TOA_ini_ceres':swtoa_ini_ceres,
                 'SWD_surf_all':swd_sur_all, 'SWD_surf_clr':swd_sur_clr, 'SWD_surf_CERES':swd_sur_ceres,
                 'Direct_all':direct_all, 'Diffuse_all':diffuse_all, 'Direct_clr':direct_clr,
                 'Diffuse_clr':diffuse_clr, 'Direct_ceres':direct_all_ceres, 'Diffuse_ceres':diffuse_all_ceres
                 }
    
    if SAVE:
        with open(path+f'CCCma_v11_SYN1_iceovercast_constant_{year[0]}_to_{year[-1]}_v6b_CERES_4_2.pkl', 'wb') as f:
            pickle.dump(fdata, f)
    print(f'Finished: {station} ',xo)





