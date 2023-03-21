import numpy as np
import matplotlib.pyplot as plt
import time

from astropy.io import fits

import os
import sys

import astropy.cosmology as cc
cosmo = cc.Planck15

import astropy.units as uu
import astropy.constants as const

from tqdm import tqdm

import ppxf.ppxf_util as util

from PyAstronomy import pyasl


snap = int(96)
gal = int(2)
redshift = float(0.03)

namegal = ('snap'+str(snap)+"gal"+str(gal))



miii = "./Data/"+namegal+"/"+namegal+'MIII.dat'
exist_miii = os.path.isfile(miii)

ver = str("v0.2")
lib = str("th")

# you need to have the MaStar model. Download it from https://www.icg.port.ac.uk/mastar/
# An updated version, the one used in Nanni+2022 paper, will be soon made public on the same webpage

hdul = fits.open('MaStar_SSP_'+ver+'.fits')
t=hdul[1].data[:,0,0,0] #parameter matrix

Z=hdul[1].data[0,:,0,1] #parameter matrix

s=hdul[1].data[0,0,:,2] #parameter matrix
fluxgrid=hdul[3].data #fluxgrid matrix

wave=hdul[2].data[0,:] #lambda array

r_instrument = hdul[2].data[1,:]
sres = r_instrument
sig2fwhm        = 2.0 * np.sqrt(2.0 * np.log(2.0))

hdul.close()

sin = float(1.3)

wave_obs = wave*(1+redshift)


def spec_mastar(Mstar, Zin_m, tin_m):
    jv_mastar = np.zeros(len(wave))

    if (Zin_m<-1.35 and tin_m<1.):
        for k in range(len(wave)):    
            jv_mastar[k]= -99.
    else:
        p_ = np.unique(np.sort(np.append( t, tin_m)))
        pR =  np.array((np.where(tin_m == p_)))[0]
        pL = pR -1
        tL = t[pL]
        tR = t[pR]
        ht = (tin_m-tL)/(tR-tL)
        
        
        m_ = np.unique(np.sort(np.append( Z, Zin_m)))
        m = np.array(np.where(Zin_m == m_))
        if m>=len(Z)-1:
            mR = len(Z)-1
            mL = mR-1
        else:
            mR = np.array(np.where(Zin_m == m_))[0]
            mL = mR-1
            
        ZL = Z[mL]
        ZR = Z[mR]
        hZ = (Zin_m-ZL)/(ZR-ZL)   

    
        s_ = np.unique(np.sort(np.append( s , sin)))
        _s = np.array(np.where(s_ == sin))[0] 
        sR = _s
        sL = sR -1
        alphaL = s[sL]
        alphaR = s[sR]
        hs = (sin-alphaL)/(alphaR -alphaL)    
        
        hpL = 1.-ht
        hpR = ht
        hmL = 1.-hZ
        hmR = hZ
        hsL = 1.-hs
        hsR = hs

        jv_mastar[:] =             hpL * hmL * hsL * fluxgrid[pL, mL, sL, :]+            hpR* hmL*hsL* fluxgrid[pR, mL, sL, :]+            hpR*hmR*hsL*fluxgrid[pR, mR, sL, :]+            hpR*hmR*hsR*fluxgrid[pR,mR,sR,:]+            hpR*hmL*hsR*fluxgrid[pR,mL,sR,:]+            hpL*hmR*hsL*fluxgrid[pL,mR,sL,:]+            hpL*hmR*hsR*fluxgrid[pL,mR,sR, :]+            hpL*hmL*hsR*fluxgrid[pL,mL,sR,:]

    return (jv_mastar*Mstar)/np.asarray(4*np.pi*np.asarray(cosmo.luminosity_distance(redshift).to(uu.cm))**2)


Nlambda = 1800
NZrel = 5
NlogC = 6
Nlogp = 5

_Zrelv = np.ones((NZrel))
_logCv= np.ones((NlogC))
_logpv= np.ones((Nlogp))
_j0vv= np.ones((NZrel,NlogC,Nlogp,Nlambda))
_j1vv= np.ones((NZrel,NlogC,Nlogp,Nlambda))
_lambdav = np.ones((Nlambda))

_Zrelv[0] = 0.05
_Zrelv[1] = 0.20    
_Zrelv[2] = 0.40     
_Zrelv[3] = 1.00
_Zrelv[4] = 2.00
Zrelnamev = ["Z005", "Z020", "Z040", "Z100", "Z200"]

_logCv[0] = 4.0
_logCv[1] = 4.5
_logCv[2] = 5.0
_logCv[3] = 5.5
_logCv[4] = 6.0     
_logCv[5] = 6.5     
logCnamev = ["C40","C45", "C50", "C55", "C60", "C65"]

_logpv[0] = 4.0  
_logpv[1] = 5.0 
_logpv[2] = 6.0    
_logpv[3] = 7.0   
_logpv[4] = 8.0      
logpnamev = ["p4", "p5", "p6","p7","p8" ]


def spec_miii(SFR,Ziii,logC,logp,fPDR):
    Zrel = Ziii
    Zrel = np.max([Ziii,0.05])
    Zrel = np.min([Ziii,2.0])

    ii = np.unique(np.sort(np.append( _Zrelv,Zrel)))
    i = np.array(np.where(Zrel == ii))[0]

    if i>=len(_Zrelv):
        iR = int(np.array(np.where(Zrel == np.max(_Zrelv)))[0])
        iL = int(iR-1)
    else:
        iR = int(i)
        iL = int(iR-1)
    hZrel = (Zrel-_Zrelv[iL])/(_Zrelv[iR]-_Zrelv[iL])


    logC =np.max([logC,4.0])
    logC = np.min([logC,6.5])

    jj = np.unique(np.sort(np.append(_logCv,logC)))
    j =  np.array((np.where(logC == jj)))[0]
    if j>=len(_logCv):
        jR = int(np.array(np.where(logC == np.max(_logCv)))[0])
        jL = int(jR-1)
    else:
        jR = int(j)
        jL = int(jR-1)

    hlogC = (logC-_logCv[jL])/(_logCv[jR]-_logCv[jL])

    logp = np.max([logp,4.0])
    logp = np.min([logp,8.0])

    kk = np.unique(np.sort(np.append(_logpv,logp)))
    k =  np.array((np.where(logp == kk)))[0]
    if k>=len(_logpv):
        kR = int(np.array(np.where(logp == np.max(_logpv)))[0])
        kL = int(kR-1)
    else:
        kR = int(k)
        kL = int(kR-1)  
    hlogp = (logp-_logpv[kL])/(_logpv[kR]-_logpv[kL])
    
    # you need to have the Mappings model 
    # the data can be download from here https://github.com/SKIRT/SKIRT8/tree/master/SKIRT/data/SED/Mappings
    j0LLLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jL]+"_"+logpnamev[kL]+".dat")[:,1]
    j0RLLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jL]+"_"+logpnamev[kL]+".dat")[:,1]
    j0LRLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jR]+"_"+logpnamev[kL]+".dat")[:,1]
    j0RRLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jR]+"_"+logpnamev[kL]+".dat")[:,1]
    j0LLRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jL]+"_"+logpnamev[kR]+".dat")[:,1]
    j0RLRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jL]+"_"+logpnamev[kR]+".dat")[:,1]
    j0LRRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jR]+"_"+logpnamev[kR]+".dat")[:,1]
    j0RRRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jR]+"_"+logpnamev[kR]+".dat")[:,1]
    j1LLLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jL]+"_"+logpnamev[kL]+".dat")[:,2]
    j1RLLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jL]+"_"+logpnamev[kL]+".dat")[:,2]
    j1LRLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jR]+"_"+logpnamev[kL]+".dat")[:,2]
    j1RRLv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jR]+"_"+logpnamev[kL]+".dat")[:,2]
    j1LLRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jL]+"_"+logpnamev[kR]+".dat")[:,2]
    j1RLRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jL]+"_"+logpnamev[kR]+".dat")[:,2]
    j1LRRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jR]+"_"+logpnamev[kR]+".dat")[:,2]
    j1RRRv = np.genfromtxt("./Mappings/Mappings_"+Zrelnamev[iR]+"_"+logCnamev[jR]+"_"+logpnamev[kR]+".dat")[:,2]

    j0 = (1.0-hZrel)*(1.0-hlogC)*(1.0-hlogp)*j0LLLv    + hZrel*(1.0-hlogC)*(1.0-hlogp)*j0RLLv    + (1.0-hZrel)*hlogC*(1.0-hlogp)*j0LRLv    + hZrel*hlogC*(1.0-hlogp)*j0RRLv    + (1.0-hZrel)*(1.0-hlogC)*hlogp*j0LLRv    + hZrel*(1.0-hlogC)*hlogp*j0RLRv    + (1.0-hZrel)*hlogC*hlogp*j0LRRv    + hZrel*hlogC*hlogp*j0RRRv;
    j1 = (1.0-hZrel)*(1.0-hlogC)*(1.0-hlogp)*j1LLLv[k]    + hZrel*(1.0-hlogC)*(1.0-hlogp)*j1RLLv[k]    + (1.0-hZrel)*hlogC*(1.0-hlogp)*j1LRLv[k]    + hZrel*hlogC*(1.0-hlogp)*j1RRLv[k]    + (1.0-hZrel)*(1.0-hlogC)*hlogp*j1LLRv[k]    + hZrel*(1.0-hlogC)*hlogp*j1RLRv[k]    + (1.0-hZrel)*hlogC*hlogp*j1LRRv[k]    + hZrel*hlogC*hlogp*j1RRRv[k];

    miii_ = ((1.0-fPDR)*j0 + fPDR*j1)*SFR*10**-3
    
    wave_miii = np.genfromtxt( "./Mappings/Mappings_"+Zrelnamev[iL]+"_"+logCnamev[jL]+"_"+logpnamev[kL]+".dat")[:,0]
    wave_miii = wave_miii*10**10
    miii = np.interp(wave, wave_miii, miii_)




    return miii/np.asarray(4*np.pi*np.asarray(cosmo.luminosity_distance(redshift).to(uu.cm))**2)


import astropy.cosmology as co
cosmo = co.Planck15
h = float(np.asarray(cosmo.h))

gal_info = np.loadtxt("./mangacat.txt")

# information from TNG50, the half mass stellar radius

HMSR = 3.3009703/h

import astropy.units as uu

dist_A = cosmo.angular_diameter_distance(redshift)
theta_spaxel = 0.5*uu.arcsec
l_spaxel = (theta_spaxel * dist_A).to(uu.pc, uu.dimensionless_angles())
FoV = 30*HMSR*10**3

def round_up_to_odd(f):
    return int(np.ceil(f) // 2 * 2 + 1)


n_pixel = round_up_to_odd(np.array(FoV/l_spaxel))
n_pixel = 300 if n_pixel > 300 else n_pixel
FoV = n_pixel * np.asarray(l_spaxel)

x_edge = np.linspace(-FoV/2., FoV/2., n_pixel+1)
y_edge = np.linspace(-FoV/2., FoV/2., n_pixel+1)

cdelt = float(np.array(l_spaxel.copy()))



data_gal = np.loadtxt("./Data/"+namegal+"/"+namegal+".dat")
x = data_gal[:,0]
y = data_gal[:,1]
z = data_gal[:,2]

Mstar = np.round(data_gal[:,4],2)
Zin = np.round(data_gal[:,5],2)
tin = np.round(data_gal[:,6],2)
data_gal_vel = np.loadtxt("./Data/"+namegal+"/"+namegal+'velocity.dat')
velz = data_gal_vel[:,2]



x_bin_id_list = np.digitize(x, x_edge)
y_bin_id_list = np.digitize(y, y_edge)



if exist_miii == True:
    data_miii = np.loadtxt("./Data/"+namegal+"/"+namegal+'MIII.dat')
    x_miii = data_miii[:,0]
    y_miii = data_miii[:,1]
    z_miii = data_miii[:,2]
    h_miii = data_miii[:,3]
    SFR = data_miii[:,4]
    Z_ = data_miii[:,5]
    logC = data_miii[:,6]
    logp = data_miii[:,7]
    fPDR = data_miii[:,8]

    data_miii = np.loadtxt("./Data/"+namegal+"/"+namegal+'starsMIII.dat')
    MstarIII = data_miii[:,0]
    ZinIII = data_miii[:,1]
    tinIII = data_miii[:,2]
    
    
    data_gal_velIII = np.loadtxt("./Data/"+namegal+"/"+namegal+'velocityMIII.dat')

    velz_III = data_gal_velIII[:,2]


    data_gas = np.loadtxt("./Data/"+namegal+"/"+namegal+'velocityGas.dat')

    x_gas = np.concatenate((data_gas[:,0], x_miii))
    y_gas = np.concatenate((data_gas[:,1], y_miii))
    
    velz_gas = np.concatenate((data_gas[:,4], velz_III))
    x_bin_id_list_gas = np.digitize(x_gas, x_edge)
    y_bin_id_list_gas = np.digitize(y_gas, y_edge)
    
    
    x_bin_id_list_miii = np.digitize(x_miii, x_edge)
    y_bin_id_list_miii = np.digitize(y_miii, y_edge)







def bin_spectrum(x_bin_id, y_bin_id):
    same_bin = np.where( np.logical_and( (x_bin_id_list==x_bin_id), (y_bin_id_list==y_bin_id) ) )
    same_bin = np.array(same_bin).flatten()
    spectra = np.zeros([len(wave),len(same_bin)])
    velz_sb = np.zeros(len(same_bin))
    spec_ = np.zeros([len(wave)])
    


    if (len(same_bin))>= 1:
        for sb in range(len(same_bin)):  
            spectra[:,sb] = spec_mastar(Mstar[same_bin[sb]], Zin[same_bin[sb]], tin[same_bin[sb]])
            velz_sb[sb] = velz[same_bin[sb]]

        mean_velz = np.nanmean(velz_sb)
        sigma_velz = np.nanstd(velz_sb)
        spec_ = pyasl.dopplerShift( wave, np.sum(spectra, axis=1) , mean_velz, vlim = np.array(c)*10**-3)[0]
        spec_[np.isnan(spec_)] = np.sum(spectra, axis=1)[np.isnan(spec_)]
        
        return spec_, sigma_velz
    else:
        return np.zeros(len(wave)), -1


    


def bin_spectrum_miii(x_bin_id, y_bin_id):
    same_bin_miii = np.where( np.logical_and( (x_bin_id_list_miii==x_bin_id), (y_bin_id_list_miii==y_bin_id) ) )

    same_bin_miii = np.array(same_bin_miii).flatten() 
    spectra_miii = np.zeros([len(wave),len(same_bin_miii)])
    spec_ = np.zeros([len(wave)])
    
    if (len(same_bin_miii))>= 1:
        same_bin_gas = np.where( np.logical_and( (x_bin_id_list_gas==x_bin_id), (y_bin_id_list_gas==y_bin_id) ) )
        same_bin_gas = np.array(same_bin_gas).flatten() 
        velz_sb = np.zeros(len(same_bin_gas))

        for sb in range(len(same_bin_miii)):  
            spectra_miii[:,sb] = spec_miii(SFR[same_bin_miii[sb]], Z_[same_bin_miii[sb]], logC[same_bin_miii[sb]], logp[same_bin_miii[sb]], fPDR[same_bin_miii[sb]]) 

        for sb in range(len(same_bin_gas)):  
           
            velz_sb[sb] = velz_gas[same_bin_gas[sb]]
        
        mean_velz = np.nanmean(velz_sb)
        sigma_velz_gas = np.nanstd(velz_sb)
    


        spec_ = pyasl.dopplerShift( wave, np.sum(spectra_miii, axis=1) , mean_velz, vlim = np.array(c)*10**-3)[0]
        spec_[np.isnan(spec_)] = np.sum(spectra_miii, axis=1)[np.isnan(spec_)]
        
        return spec_, sigma_velz_gas
    else:
        return np.zeros(len(wave)), -1



if exist_miii == False:
    arr = np.zeros((len(wave),n_pixel,n_pixel))
    std_vz = np.zeros((n_pixel,n_pixel))

    for i in tqdm(range(n_pixel)):
        for j in range(n_pixel):
            arr[:,i,j], std_vz[i,j] = bin_spectrum(i+1, j+1)
         


else:
    arr = np.zeros((len(wave),n_pixel,n_pixel))
    arr_miii = np.zeros((len(wave),n_pixel,n_pixel))
    std_vz = np.zeros((n_pixel,n_pixel))
    std_vz_gas = np.zeros((n_pixel,n_pixel))


    for i in tqdm(range(n_pixel)):
        for j in range(n_pixel):
            arr[:,i,j], std_vz[i,j] = bin_spectrum(i+1, j+1)
            arr_miii[:,i,j], std_vz_gas[i,j] = bin_spectrum_miii(i+1, j+1)
    
    std_vz_gas[(np.isnan(std_vz_gas))] = 0.1
    std_vz_gas[(std_vz_gas<0.1)] = 0.1


std_vz[(np.isnan(std_vz))] = 0.1
std_vz[(std_vz<0.1)] = 0.1

c = np.asarray(c)*10**-3



if exist_miii == False:
    grid = (arr[2000,:,:]>0.) 
    arr_down = np.zeros(np.shape(arr))
    sig_instrument = c/r_instrument/sig2fwhm

    for i in tqdm(range(n_pixel)):
        for j in range(n_pixel):
            if ((grid[i, j]==True)):
                sigma_gal = std_vz[i,j]
                new_sigma = np.sqrt(sigma_gal**2.0 +sig_instrument**2.0)
                new_fwhm    = wave * new_sigma/c
                arr_down[:,i, j] = util.gaussian_filter1d(arr[:,i, j], new_fwhm)
    
    hdu = fits.PrimaryHDU(arr_down[:, : , : ])
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("./"+namegal+'/'+namegal+'.fits')
    hdulist.close()
         


else:
    grid = (arr[2000,:,:]>0.) 
    
    grid_miii = (arr_miii[800,:,:]>0.) 
    

    arr_down_miii = np.zeros(np.shape(arr))
    sig_instrument = c/r_instrument/sig2fwhm

    for i in tqdm(range(n_pixel)):
        for j in range(n_pixel):
            if ((grid_miii[i, j]==True)):
                sigma_gal = std_vz_gas[i,j]
                new_sigma = np.sqrt(sigma_gal**2.0 +sig_instrument**2.0)
                new_fwhm    = wave * new_sigma/c
                arr_down_miii[:,i, j] = util.gaussian_filter1d(arr_miii[:,i, j], new_fwhm)

    arr_down = np.zeros(np.shape(arr))
    for i in tqdm(range(n_pixel)):
        for j in range(n_pixel):
            if ((grid[i, j]==True)):
                sigma_gal = std_vz[i,j]
                new_sigma = np.sqrt(sigma_gal**2.0 +sig_instrument**2.0)
                new_fwhm    = wave * new_sigma/c
                arr_down[:,i, j] = util.gaussian_filter1d(arr[:,i, j], new_fwhm)



    arr_comb = arr_down_miii+ arr_down

    hdu = fits.PrimaryHDU(arr_comb)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("./"+namegal+'/'+namegal+'.fits')
    hdulist.close()



