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


class ContourProducer(object):
    """ documentation """

    def __init__(self,  path = "./", RA_DEC_list="./ra_declist.csv"):

        # All the class properties we need
        self.coordlist = None # pandas table with coordinates
        self.path = None  # main path, all output will be in subdirectories
        self.image_directory = None  # directory for the image files
        self.download_directory = None  # directory for the downloaded drz files

        # Loading what we can
        self.load_coords(RA_DEC_list)
        self.setup_output_directories(path)
        

       
    def setup_output_directories(self, path, verbose=True):
        """
        Given a path sets up main directory and image directory.
        """

        self.path = path
        
        self.image_directory = os.path.join(self.path,"Contour-Images")
        os.mkdir(self.image_directory)

        self.download_directory = os.path.join(self.path,"Download_DRZs")
        os.mkdir(self.download_directory)

        if verbose:
            print("Successfully created image and download directories")

            
        
    def load_coords(self, ra_dec_file, verbose=True):
        """
        Loads ra/dec corrdinates from a file.
        TODO: info about required format of file

        Parameters
        ----------
        ra_dec_file : str
            The file containing the ra/dec values
        verbose : bool
            Optional, default True. Set to False to suppress print output.
        """

        if verbose:
            print("Loading coordinates from {}".format(ra_dec_file))

        self.coordlist = pd.read_csv(ra_dec_file)
        self.coordlist = self.coordlist.drop_duplicates()
      

    def get_datafiles(self, ra, dec, verbose=True):
        """
        Given an ra and dec create the correct output directory,
        and download associated drizzle files
        """
        coords = "{}_{}".format(ra, dec)

        if verbose:
            print(coords)

        # Querying MAST and downloading all associated drz files
        table_of_observations = Observations.query_criteria(coordinates=coords, radius=0.02,
                                                            obs_collection="HST", dataproduct_type="image")
        
        data_product_list = Observations.get_product_list(table_of_observations)
        manifest = Observations.download_products(data_product_list, extension="drz.fits", obs_collection="HST",
                                                  download_dir=self.download_directory)
        
        return list(manifest["Local Path"])  # We just want the local file locations 
        
        
    def produce_contour_img(self, ra, dec, drz_files, verbose=True):
        """
        Given ra/dec and list of drizzle files, produce contour image
        """

        coords = "{}_{}".format(ra, dec)
        imageDirName = os.path.join(self.image_directory, coords)

        if verbose:
            print("created directory {}".format(imageDirName))

        output_imgs = []
            
        for path_name in drz_files:
            hdulist = fits.open(path_name, memmap=False)
            hdu = hdulist[1]      
            wcs = WCS(hdu.header)
            pixCoords=SkyCoord([ra], [dec], frame='fk5', unit=u.deg).to_pixel(wcs=wcs, mode='all')

            # getting the file name for the output image
            _, img_file = os.path.split(path_name)
            img_file = os.path.join(imageDirName, img_file.replace("fits","png"))
            
            if 0 < pixCoords[0][0] < hdu.header['NAXIS1'] and 0 < pixCoords[1][0] < hdu.header['NAXIS2']:
                if not os.path.exists(imageDirName):
                    os.mkdir(imageDirName) # only make the directory if we need it
                
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
                plt.savefig(fname=img_file, bbox_inches='tight')
                plt.clf()
                close()

                output_imgs.append(img_file)
                           
            else:
                if verbose:
                    print("Target not in image")
                        
            hdulist.close()

        return output_imgs

            
    def clean_download_dir(self):
        """
        This will remove all subdirectories and files from the download directory, 
        so use with caution.
        """
        shutil.rmtree(self.download_directory)
        os.mkdir(self.download_directory)  # This seems hacky but I can't figure out a cleaner way


    def make_all_contours(self, verbose=True):
        """
        Go through the whole coordlist and make contour images
        """

        contour_imgs = []
        
        for n in range(len(self.coordlist)):
            # pull the coordinate pairs from the coord data file
            ra = self.coordlist.iloc[n,0]
            dec = self.coordlist.iloc[n,1]
            
            try:
                filelist = self.get_datafiles(ra, dec, verbose)
                contour_imgs += self.produce_contour_img(ra, dec, filelist, verbose)
                
            except (ValueError, TimeoutError): #in case an error occurs, wait a few seconds, delete the folder created.
                if verbose:
                    print("Error occurred with coord {}, {}".format(ra,dec))
                sleep(4)

            self.clean_download_dir()

            if verbose:
                print("Completed coordinate {}, {}\n".format(ra, dec))

        return contour_imgs

        
    
    
