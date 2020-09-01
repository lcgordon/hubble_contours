# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 09:37:27 2020

@lcgordon


"""
import pandas as pd
import numpy as np
import astroquery
import astropy
import os
import shutil

from astroquery.mast import Observations
from astroquery.mast import Catalogs

%matplotlib inline
from pylab import *
import matplotlib
import matplotlib.pyplot as plt

from astropy.io import ascii
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

import os, fnmatch
import subprocess
from subprocess import call
import re
from numpy import arange

from time import sleep


class contour_producer(object):
    """ documentation """
    
    def __init__(self, path = "./", RA_DEC_list="./ra_declist.csv"):
        print("Producing Contours From List")
        #read in list of coords
        self.coordlist = pd.read_csv(RA_DEC_list)
        self.coordlist = self.coordlist.drop_duplicates()
        #empty object for later
        #distVec = []
        self.path = path
        self.image_directory = self.path + "Contour-Images/"
        os.mkdir(self.image_directory)
        print("Successfully created image directory")

        for n in range(len(self.coordlist)):
            try:
                #pull the coordinate pairs from the csv file
                self.RA_float = self.coordlist.iloc[n,0]
                RA = str(self.RA_float)
                self.DEC_float = self.coordlist.iloc[n,1]
                DEC = str(self.DEC_float)
                self.coords = RA + "_" + DEC
                self.download_dir = self.path + self.coords
            
                #search for and download all data products associated with those coordinates
                print("coordinates:", self.coords, " n:", n)     
                self.table_of_observations = Observations.query_region(self.coords, radius = 0.02)
                self.data_product_list = Observations.get_product_list(self.table_of_observations)
                downloads = Observations.download_products(self.data_product_list, 
                                                           download_dir=self.download_dir,
                                                           obs_collection = "HST",
                                                           dataproduct_type = "image",
                                                           productType = "SCIENCE")
            
                #Filter out all files that are not drz.fits files (drizzle reduced)--->
                #there is probably an easier way to do this using os.endswith()
                filelist = []
                for root, dirs, files in os.walk(self.download_dir):
                    for name in files:
                        #print(name)
                        if name.endswith(("drz.fits")):
                            filelist.append(root + "/" + name)
                #findQSOdrz = 'find . -type f \( -name "*drz.fits" \)'
                #out = subprocess.Popen(findQSOdrz, shell=True, stdin=subprocess.PIPE, 
                 #                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #(stdout, stderr) = out.communicate()
                #filelist = stdout.decode().split()
            
                length = len(filelist)
                print("drz files found:", length)
            
                #produce contour plots for all available drz images
                if length > 0:
                    imageDirName = self.image_directory + self.coords       #making image directory and pulling RA/Dec from file name 
                    os.mkdir(imageDirName)
                    print("created directory {}".format(imageDirName))
                    for i in range(length):
                        file_location = filelist[i]
                        m = re.search('HST/(.+?)/', file_location)
                        if m:
                            file_name = m.group(1)
                            print(file_name)
                        p = re.search('HSThits/(.+?)_', file_location)
                        if p:
                            RA = p.group(1)
                        q = re.search('_(.+?)/mastDownload', file_location)
                        if q:
                            DEC = q.group(1)
                            
                        world_coords = RA + "_" + DEC
            
                        path_name = filelist[i]
                        print(path_name)
            
                        hdu = fits.open(path_name, memmap=False)[1]       #converting world coords to pixel coords
                        wcs = WCS(hdu.header)
                        pixCoords=SkyCoord([RA], [DEC], frame='fk5', unit=u.deg).to_pixel(wcs=wcs, mode='all')
            
                        if 0< pixCoords[0][0] < hdu.header['NAXIS1'] and 0 < pixCoords[1][0] < hdu.header['NAXIS2']:
                            fig = plt.figure()                        #plotting only quasars that lie in image
                            fig.add_subplot(111, projection=wcs)
                            plt.imshow(hdu.data, origin='lower', cmap='gray',clim=(0,.99))
                            plt.xlabel('RA')
                            plt.ylabel('Dec')
                            plt.grid(color='white', ls='solid')
                            plt.contour(hdu.data, levels=[.1,.5,1,4,7,10,20], cmap='cool', alpha=0.5)
            
                            pixCoords=np.concatenate(pixCoords)            #cropping image to only quasar
                            plt.xlim(pixCoords[0]-50, pixCoords[0]+50)
                            plt.ylim(pixCoords[1]-50, pixCoords[1]+50)
                            plt.savefig(fname = file_name ,bbox_inches='tight')
                            pic_name = file_name + ".png"
                            plt.clf()
                            close()
                            shutil.move(pic_name, imageDirName)            #moving the image into its specific directory
                        
                            #distVec.append(((pixCoords[0]-hdu.header['CRPIX1'])/hdu.header['NAXIS1']+(pixCoords[1]-hdu.header['CRPIX2'])/hdu.header['NAXIS2'])/2)
                            
                        else:
                            print("Target not in image")
                        
                        hdu.close()
                        
                    if os.path.isdir(self.download_dir) == True:
                        shutil.rmtree(self.download_dir)               #deletes ALL data to conserve space
                
                else: 
                    print("No data")
                    if os.path.isdir(self.download_dir) == True:
                        #deleting the folder in case no drz files exist
                        shutil.rmtree(self.download_dir)
                
                if os.path.isdir(imageDirName) == True and len(os.listdir(imageDirName)) == 0:
                    #deleting empty image directories
                    os.rmdir(imageDirName)               
            
            except (ValueError, TimeoutError): #in case an error occurs, wait a few seconds, delete the folder created.
                print("error occurred")
                sleep(4)
                if os.path.isdir(self.download_dir) == True:
                    shutil.rmtree(self.download_dir)