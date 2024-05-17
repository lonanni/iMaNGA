#!/usr/bin/env python
# coding: utf-8

# In[4]:


import requests
import numpy as np




# In[2]:


baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"b65ee99582d40446ede7aa5ed7d79ac4"}

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r





def apply_FOV(grid, datacube):
    grid = np.where(grid>-1, True, False)
    #reshaping to apply to the larger synthetic datacube
    grid_300 = np.zeros(shape=(300,300))*False
    grid_300[75:225, 75:225] = grid 
    FoV_datacube = np.where(grid_300==True, datacube, float("Nan"))
    return grid_300, FoV_datacube





