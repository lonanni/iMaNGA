#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import numpy as np

from PyAstronomy import pyasl
import astropy.cosmology as cc
cosmo = cc.Planck15
import astropy.units as uu
from astropy.constants import c


# In[2]:


ver = str("v0.2")
lib = str("th")

hdul = fits.open('MaStar_SSP_'+ver+'.fits')
t=hdul[1].data[:,0,0,0] #parameter matrix

Z=hdul[1].data[0,:,0,1] #parameter matrix

s=hdul[1].data[0,0,:,2] #parameter matrix
fluxgrid=hdul[3].data #fluxgrid matrix

wave=hdul[2].data[0,:] #lambda array
hdul.close()

sin = float(1.3)
    


# In[3]:


def spec_mastar(redshift, Mstar, Zin_m, tin_m):
        """
    synthetic MaStar spectrum for a given stellar particle
    
    Parameters
    ----------
    redshift : float
        redshift of the simulated galaxies
        
    Mstar : float
        stellar mass in solar masses M_sun
        
    Zin_m : float
        stellar metallicity [Z/H]
        
    tin_m : float
        stellar age [Gyr]
        
        
        
    Returns
    -------
    array
        stellar spectrum
        
    """
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


# In[4]:


def bin_spectrum(redshift, Mstar, Zin, tin, velz, x_bin_id_list,y_bin_id_list, x_bin_id, y_bin_id):
            """
        return the sum of all the stellar spectra (as returned by spec_mastar) in a spaxel 
   
   Parameters
    ----------
    redshift : float
        redshift of the simulated galaxies
        
    Mstar : float
        stellar mass in solar masses M_sun
        
    Zin_m : float
        stellar metallicity [Z/H]
        
    tin_m : float
        stellar age [Gyr]
        
    velz : float
        stellar velocity along the line-of-sight [km/s]
        
    x_bin_id_list : int
    
        spatial grid 
        
    y_bin_id_list : int
    
        spatial grid 
        
        
    x_bin : int
    
        spatial position 
        
    y_bin : int
    
        spatial position 
        
        
    
        
        
        
    Returns
    -------
    array, array
        stellar spectrum, velocity dispersion
        
    """
        
    same_bin = np.where( np.logical_and( (x_bin_id_list==x_bin_id), (y_bin_id_list==y_bin_id) ) )
    same_bin = np.array(same_bin).flatten()
    
    
    spectra = np.zeros([len(wave),len(same_bin)])
    velz_sb = np.zeros(len(same_bin))
    spec_ = np.zeros([len(wave)])
    
    


    if (len(same_bin))>= 1:
        for sb in range(len(same_bin)):  
            spectra[:,sb] = spec_mastar(redshift, Mstar[same_bin[sb]], Zin[same_bin[sb]], tin[same_bin[sb]])
            velz_sb[sb] = velz[same_bin[sb]]

        mean_velz = np.nanmean(velz_sb)
        sigma_velz = np.nanstd(velz_sb)
        spec_ = pyasl.dopplerShift( wave, np.sum(spectra, axis=1) , mean_velz, vlim = np.array(c)*10**-3)[0]
        spec_[np.isnan(spec_)] = np.sum(spectra, axis=1)[np.isnan(spec_)]
        
        return spec_, sigma_velz
    else:
        return np.zeros(len(wave)), -1


# In[5]:


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


# In[1]:


def spec_miii(redshift, SFR,Ziii,logC,logp,fPDR):
            """
    synthetic MIII spectrum for a given stellar particle
    
    Parameters
    ----------
    redshift : float
        redshift of the simulated galaxies
        
    SFR : float
        star formation rate [M_sun/yr]
        
    logC : float
        compactness of the HII region
        
    logp : float
        pressure of the HII region
        
    fPDR: float
    
        fPDR of the HII region
        
    see Nanni et al 2022 for the definition of these parameters
        
        
        
    Returns
    -------
    array
        stellar spectrum
        
    """
        
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


# In[5]:


def bin_spectrum_miii(redshift, SFR,Ziii,logC,logp,fPDR, velz_miii, x_bin_id_list_miii, y_bin_id_list_miii, x_bin_id, y_bin_id):
    
    
            """
        return the sum of all the stellar spectra (as returned by spec_miii) in a spaxel 
   
   Parameters
    ----------
    redshift : float
        redshift of the simulated galaxies
        
    SFR : float
        star formation rate [M_sun/yr]
        
    logC : float
        compactness of the HII region
        
    logp : float
        pressure of the HII region
        
    fPDR: float
    
        fPDR of the HII region
        
    see Nanni et al 2022 for the definition of these parameters
        
        
    velz_miii : float
        stellar velocity along the line-of-sight [km/s]
        
    x_bin_id_list : int
    
        spatial grid 
        
    y_bin_id_list : int
    
        spatial grid 
        
        
    x_bin : int
    
        spatial position 
        
    y_bin : int
    
        spatial position 
        
        
    
        
        
        
    Returns
    -------
    array, array
        stellar spectrum, velocity dispersion
        
    """
    same_bin_miii = np.where( np.logical_and( (x_bin_id_list_miii==x_bin_id), (y_bin_id_list_miii==y_bin_id) ) )

    same_bin_miii = np.array(same_bin_miii).flatten() 
    spectra_miii = np.zeros([len(wave),len(same_bin_miii)])
    spec_ = np.zeros([len(wave)])
    velz_sb = np.zeros(len(same_bin_miii))
    
    if (len(same_bin_miii))>= 1:


        for sb in range(len(same_bin_miii)):  
            spectra_miii[:,sb] = spec_miii(redshift, SFR[same_bin_miii[sb]], Ziii[same_bin_miii[sb]], logC[same_bin_miii[sb]], logp[same_bin_miii[sb]], fPDR[same_bin_miii[sb]])            
            velz_sb[sb] = velz_miii[same_bin_miii[sb]]
        
        mean_velz = np.nanmean(velz_sb)
        sigma_velz_gas = np.nanstd(velz_sb)
    


        spec_ = pyasl.dopplerShift( wave, np.sum(spectra_miii, axis=1) , mean_velz, vlim = np.array(c)*10**-3)[0]
        spec_[np.isnan(spec_)] = np.sum(spectra_miii, axis=1)[np.isnan(spec_)]
        
        return spec_, sigma_velz_gas
    else:
        return np.zeros(len(wave)), -1




