"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- apply the calibrations (master bias, master dark and master flat)
"""

from ccdproc import ImageFileCollection, subtract_bias
from .headers import HeaderManipulation
from ..collection_manager import CollectionManager
from astropy.stats import mad_std, sigma_clipped_stats, median_absolute_deviation
import numpy as np
from astropy.nddata import CCDData
from astropy.io import fits

class Calibration():

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
        

    def apply_bias_correction(self):
        """
        Apply to bias correction to dark, flat and light type images. It can recognize different types of images.

        Parameters
        ----------

        Returns
        -------
        self._image_collection
        """
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)

        # get the collection of the master bias
        master_bias_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Master Bias"]})
        master_bias_file = master_bias_collection.files
        master_bias_path = master_bias_collection.location/master_bias_file[0]

        # get the collection of the images to be bias corrected
        correct_collection = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": ["Dark", "Light", "Flat"]})
        files_to_correct = correct_collection.files
        
        if len(master_bias_file) == 0:
            raise ValueError("You haven't combined bias yet!")
            
        elif len(master_bias_file) > 1:
            raise ValueError(f"You have more than two master bias frames --> {master_bias_file}")
            
        else:
            print("Start Bias correction....")
            master_bias_ccd = CCDData.read(master_bias_path)
            
            # bias correction
            for ccd, fname in correct_collection.ccds(return_fname=True):
                print(f"Apply {master_bias_file} correction to {fname}")
                ccd = subtract_bias(ccd, master_bias_ccd)
                ccd.write(self._location / fname, overwrite=True)

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

        
            
        