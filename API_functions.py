#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
import numpy as np
import h5py

baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key" : "YOUR_KEY_FROM_ILLUSTRISTNG"}


# In[2]:


def get(path, params=None, fName='temp'): # gets data from url, saves to file
    """
    Routine to pull data from online
    """
    
    # make HTTP GET request to path
    if (len(headers['api-key'])!=32):
        print("Have you put in your API key? This one isn't working")
        print('Currently it is: ',headers['api-key'])
        print('You can find your API key on the Illustris website:')
        print('http://www.illustris-project.org/data/')
        print('and update it in this program using')
        print("iApi.headers['api-key']='*MYAPIKEY*'")
        print("Or permanently change it in iApi.py")
    r = requests.get(path, params=params, headers=headers)
    
    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    dataFile=fName+'.hdf5'
    # Saves to file, currently disabled
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(dataFile, 'wb') as f:
            f.write(r.content)
        return dataFile # return the filename string

    return r


# In[3]:


def getGalaxy(gal, fields, simulation='Illustris-1', snapshot=135,
              fileName='tempGal', rewriteFile=1, getHalo=0):

    fields=np.array(fields) # converts to array
    order=np.argsort(fields[:,0])
    disorder=np.argsort(order) # needed to unsort the fields later...
    fields=fields[order,:] # orders by particle type
    nFields=order.size
    
    if rewriteFile==1: # redownloads file from the internet
        url='https://www.tng-project.org/api/'+simulation+'/snapshots/'+str(snapshot)+'/subhalos/'+str(gal)+'/cutout.hdf5?'
        if getHalo==1:
            url='https://www.tng-project.org/api/'+simulation+'/snapshots/'+str(snapshot)+'/halos/'+str(gal)+'/cutout.hdf5?'
        
        thisParticle=0
        thisEntry=0
        firstParticle=1
        while thisParticle<6: # cycles through all particle type

            if (int(fields[thisEntry,0])!=thisParticle): # checks there is at least one field for this particle
                thisParticle+=1
                continue
            if firstParticle==1: # first particle requires no ampersand
                firstParticle=0
            else: # all later particles do
                url+='&' 
            particleTypeNames=['gas','dm','error','tracers','stars','bhs']                 
            url+=particleTypeNames[thisParticle]+'=' # adds the name of the particle type
        
            firstEntry=1
            while int(fields[thisEntry,0])==thisParticle:
                if firstEntry==1: #first entry requires no comma
                    firstEntry=0
                else: # all later entries do
                    url+=','
                url+=fields[thisEntry,1] # adds every associated field
                thisEntry+=1
                if thisEntry==nFields:
                    break
            if thisEntry==nFields:
                break
            thisParticle+=1
        dataFile=get(url,fName=fileName)
    # end of "if rewriteFile==1:"
    if rewriteFile == 0: # if we're not redownloading need to set path to the file
        dataFile=fileName+'.hdf5'
    
    # gets a dictionary for unit conversions
    #units=changeUnits.makeParticleDict(simulation=simulation,snapshot=snapshot)
    
    # actually get the data (saved to .hdf5 file)
    data=[] # initially empty list that we will fill up with the data
    with h5py.File(dataFile,'r') as f:
        for i in range(disorder.size):
            thisField=fields[disorder[i],:] # ensures data returned in original order of fields
            
            #!!!!!!!!!
            #            data.append(units[thisField[1]]*np.array(f['PartType'+thisField[0]][thisField[1]]))
            data.append(np.array(f['PartType'+thisField[0]][thisField[1]]))
            # returns all particle data of each field as a numpy array
    return data # returns all the particle fields as a list of numpy arrays in the same order as initial fields


# In[4]:


def getSubhaloField(field, simulation='Illustris-1', snapshot=135,
                    fileName='tempCat', rewriteFile=1):

    
    if rewriteFile==1: # redownloads file from the internet
        url='https://www.tng-project.org/api/'+simulation+'/files/groupcat-'+str(snapshot)+'/?Subhalo='+field
        dataFile=get(url,fName=fileName)
    if rewriteFile == 0: # if we're not redownloading need to set path to the file
        dataFile=fileName+'.hdf5'
        
    with h5py.File(dataFile,'r') as f:
        data=np.array(f['Subhalo'][field])
    return data
    
  


# In[ ]:




