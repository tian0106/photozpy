"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- apply the calibrations (master bias, master dark and master flat)
"""

from ccdproc import ImageFileCollection, subtract_bias, subtract_dark
from .headers import HeaderManipulation
from ..collection_manager import CollectionManager
from astropy.stats import mad_std, sigma_clipped_stats, median_absolute_deviation
import numpy as np
from astropy.nddata import CCDData
from astropy.io import fits
import astropy.units as u
from pathlib import Path

class Reduction():

    def __init__(self, image_collection):

        """
        Initial inputs for calibration.
        
        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; This image collection isn't the one that will be used for calibration. 
                                                       It will be refreshed when initiated or calibrated.

        Returns
        -------

        """

        self._image_collection = CollectionManager.refresh_collection(image_collection, rescan = True)
        self._location = self._image_collection.location
        

    def apply_bias_correction(self, save_location = ""):
        """
        Apply to bias correction to dark, flat and light type images. It can recognize different types of images.

        Parameters
        ----------

        Returns
        -------
        self._image_collection
        """

        if save_location == "":
            save_location = self._location
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        # get the collection of the master bias
        master_bias_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Master Bias"]})
        master_bias_path = master_bias_collection.files_filtered(include_path = True)  
        if len(master_bias_path) == 0:
            raise ValueError("You haven't a combined bias yet!")
        elif len(master_bias_path) > 1:
            raise ValueError(f"You have more than two master bias frames --> {master_bias_path}!")
            
        master_bias_path = Path(master_bias_path[0])
        master_bias_file = master_bias_path.name
        master_bias_ccd = CCDData.read(master_bias_path)

        # get the collection of the images to be bias corrected
        correct_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Dark", "Light", "Flat"]})
        
        # bias correction
        for ccd, fname in correct_collection.ccds(return_fname=True):
            print(f"Apply {master_bias_file} correction to {fname}")
            ccd = subtract_bias(ccd, master_bias_ccd)
            ccd.write(save_location / fname, overwrite=True)

        # refresh image collection
        correct_collection = CollectionManager.refresh_collection(correct_collection)

        # add bias correction header. Because the CCDData writes mask and uncertainty into the hdul so ImageFileCollection
        # iteration won't work. (it doesn't work on the multi-extension fits files!
        for i in correct_collection.files_filtered(include_path = True):
            with fits.open(i, mode = "update") as hdul:
                hdul[0].header["BIASCORR"] = "Yes"
            hdul.flush()

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)
        print("Bias Correction completed!")
        print("----------------------------------------\n")

        return

    def apply_dark_correction(self, save_location = ""):

        """
        Apply master dark correction to flat and light images

        Parameters
        ----------

        Returns
        -------
        None

        """

        if save_location == "":
            save_location = self._location

        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        # Get the master dark collection
        master_dark_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Master Dark"]})
        master_dark_path = master_dark_collection.files_filtered(include_path = True)
        if len(master_dark_path) == 0:
            raise ValueError("You don't have a master dark frame yet!")
        elif len(master_dark_path) > 1:
            raise ValueError(f"You have more than one master dark frame --> {master_dark_path}!")
        # read the ccd data of the master dark
        master_dark_path = Path(master_dark_path[0])
        master_dark_file = master_dark_path.name
        master_dark_ccd = CCDData.read(master_dark_path)

        # get the image collection of flat and light
        correct_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Light", "Flat"]})        

        # dark correction
        # since the flat exposure is much much smaller than the dark exposure
        # we have to scale the darks
        for ccd, fname in correct_collection.ccds(return_fname = True):
            print(f"Apply {master_dark_file} correction to {fname}")
            ccd = subtract_dark(ccd, master_dark_ccd, 
                                exposure_time = "EXPTIME",exposure_unit = u.second, 
                                scale = True)

            ccd.write(save_location / fname, overwrite=True)
            
        for i in correct_collection.files_filtered(include_path = True):
            with fits.open(i, mode = "update") as hdul:
                hdul[0].header["DARKCORR"] = "Yes"
                hdul.flush()


        
            
        