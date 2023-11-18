"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Combine a series of registered images.
- Check image types before combine.
"""

from ccdproc import ImageFileCollection
from ccdproc import combine as ccdpro_combine
from .headers import HeaderCorrection
from .headers import HeaderManipulation
from ..collection_manager import CollectionManager
from astropy.stats import mad_std, sigma_clipped_stats, median_absolute_deviation
import numpy as np

class Combine():
    
    """
    Combine different types of images.
    """
    
    def __init__(self, image_collection, telescope):
        """
        Init the instance.
        
        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        telescope:
        save_location: str; the location to save the combined images
        overwrite; boolean; set True to overwrite the original files if the save_location is "" or the original path
        """

        self._image_collection = CollectionManager.refresh_collection(image_collection)
        self._telescope = telescope
        self._location = self._image_collection.location
        
        
    @staticmethod
    def inv_median(a):
        return 1/np.median(a)

    @staticmethod
    def check_imtype(image_collection, imtype):
        """
        Check if the imtype of the image_collection is the same as the one we input.

        Parameters
        ----------
        image_collection: ccdproc.ImageFileCollection; the collection of images
        imtype: str; the value for IMTYPE header.

        Returns
        -------
        final: boolean; True or False
        """
        
        headers_values = {"IMTYPE": imtype}

        imtype_value = HeaderManipulation.check_headers(image_collection, **headers_values)
        
        return imtype_value

    def combine_bias_or_dark(self, image_type, method = "sigma clip", 
                             sigma_clip_low=5, sigma_clip_high=5, quite = False):
        """
        Combine bias or dark images

        Parameters
        ----------
        image_type; str; Flat or Light
        method:
        sigma_clip_low
        sigma_clip_high

        Returns
        -------
        """
        
        # refresh the full collection
        self._image_collection = CollectionManager.refresh_collection(self._image_collection)

        # Get the collection to combine and the object names
        collection_to_combine = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": [image_type]})
        object_names = HeaderManipulation.get_header_values(collection_to_combine, header = "object")  # It seems only the lower cases work
        # note that object_names inludes all the objects! (targets, bias ,dark and flat)

        for object_name in object_names:
            # this is the list of the absolute path to the fits images to be combine
            image_list_to_combine = list(collection_to_combine.files_filtered(object = object_name, include_path = True))
            image_names_to_combine = list(collection_to_combine.files_filtered(object = object_name, include_path = False))
            image_number = len(image_list_to_combine)

            if image_number == 0:
                print(f"No {object_name} found! Skipping......")

            elif image_number == 1:
                print(f"Only one {object} image! Just copying and renaming the file.")
                
                for hdu in collection_temp.hdus(save_with_name = "_master", save_location = self._location, overwrite = True):
                    hdu.header["IMTYPE"] = f"Master {image_type}"
                    #hdu.header["NCOMBINE"] = int(image_number)  It seems already added into the combiner of ccdproc
                
                fname_stem = Path(image_list_to_combine[0]).stem  # get the filename stem: crab_sdss_g.fits --get--> crab_sdss_g
                fname_suffix = Path(image_list_to_combine[0]).suffix
                name_path = Path(self._location) / f"{fname_stem}_master.{fname_suffix}"  # the path to the saved file
                new_name_path = name_path.with_name(f"Master_{object_name}.{fname_suffix}")  # the new path
                name_path.rename(new_name_path)  # rename the file, notice that you have to use the full path to rename, not just the file name

            elif image_number > 1:
                if quite == False:
                    print(f"Using {image_names_to_combine} to combine the {object_name} image!")
                if method == "sigma clip":
                    combined = ccdpro_combine(image_list_to_combine,  
                                              method = "average",
                                              sigma_clip=True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high,
                                              sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                elif method == "minmax clip":
                    raise ValueError("This function isn't ready yet!")

                else:
                    raise ValueError("Unsupported combining method! What supported now is only sigma clip!")
                
                combined.meta["IMTYPE"] = f"Master {image_type}"
                #combined.meta["NCOMBINE"] = int(number)
                combined.write(self._location / f'Master_{object_name}.fits', overwrite=True)
                print(f"{image_type} combined!")
                print("----------------------------------------\n")


        # it will return an image collection that contains all the fits files including the created master bias
        new_image_collection = CollectionManager.refresh_collection(self._image_collection, rescan = True)
        new_image_collection = ImageFileCollection(location = self._location, glob_include = "*.fits")
        self._image_collection = new_image_collection

        return


    def combine_light_or_flat(self, image_type, method = "sigma clip", scale_funtion = "auto", 
                              sigma_clip_low=5, sigma_clip_high=5, save_location = "", overwrite = True):
        """
        Combine light or flat images

        Parameters
        ----------
        image_type; str; Flat or Light

        Returns
        -------
        """

        filters = telescope.filters

        if save_location == "":
            save_location = self._image_collection.location

        # Get the object names
        object_names = HeaderManipulation.get_header_values(self._image_collection, header = "object")  # It seems only the lower cases work
        old_location = self._image_collection.location
    

        for object_name in object_names:
            for filter in filters:
                # make sure only the fits image with the correct image_type requested. It also filters the object and filters.
                collection_temp = CollectionManager.filter_collection(self._image_collection, **{"IMTYPE": [image_type], "OBJECT": [object_name], "FILTER": [filter]})
                # this is the list of the absolute path to the fits images to be combine
                image_list_to_combine = list(collection_temp.files_filtered(object = object_name, filter = filter, include_path = True))
                image_number = len(image_list_to_combine)
                
                if image_number == 0:
                    print(f"No {object_name} in {filter} filter found! Skipping......")
                    
                elif image_number == 1:
                    print(f"Only one {object} image in {filter} ({image_files})! Just copying and renaming the file.")
                    
                    for hdu in collection_temp.hdus(save_with_name = "_master", save_location = save_location, overwrite = overwrite):
                        hdu.header["IMTYPE"] = f"Master {image_type}"
                        #hdu.header["NCOMBINE"] = int(image_number)  It seems already added into the combiner of ccdproc
                        
                        if image_type == "Light":  # need to update rdnoise and gain, I don't know why it's only for light images
                            rdnoise = hdu.header["RDNOISE"]
                            hdu.header["RDNOISE"] = rdnoise/np.sqrt(image_number)  # update the rdnoise
                            gain = hdu.header["GAIN"]
                            hdu.header["GAIN"] = gain*image_number  # update the gain
                            
                    fname_stem = Path(image_list_to_combine[0]).stem  # get the filename stem: crab_sdss_g.fits --get--> crab_sdss_g
                    fname_suffix = Path(image_list_to_combine[0]).suffix
                    name_path = Path(save_location) / f"{fname_stem}_master.{fname_suffix}"  # the path to the saved file
                    new_name_path = name_path.with_name(f"Master_{object_name}_{filter}.{fname_suffix}")  # the new path
                    name_path.rename(new_name_path)  # rename the file, notice that you have to use the full path to rename, not just the file name
                    print("----------------------------------------------------------\n")

                elif image_number > 1:
                    print(f"Using {image_list_to_combine} to combine the master {object_name} image in {filter} filter!")
                    if method == "sigma clip":
                        if image_type == "Flat":
                            if scale_function == "auto":
                                scale_function = inv_median.__func__()
                            
                            combined = ccdpro_combine(image_list_to_combine, 
                                                      method = "average", scale = scale_function,
                                                     sigma_clip=True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high,
                                                     sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                        
                        elif image_type == "Light":
                            combined = ccdpro_combine(image_list_to_combine, 
                                                      method = "average", 
                                                      sigma_clip = True, sigma_clip_low_thresh=sigma_clip_low, sigma_clip_high_thresh=sigma_clip_high, 
                                                      sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
                    elif method == "minmax clip":
                        raise ValueError("This function isn't ready yet!")

                    else:
                        raise ValueError("Unsupported combining method! What supported now is only sigma clip!")

                    combined.meta["IMTYPE"] = f"Master {image_type}"
                    combined.meta["NCOMBINE"] = int(image_number)
                    if image_type == "Light":
                        rdnoise = combined.meta["RDNOISE"]
                        combined.meta["RDNOISE"] = rdnoise/np.sqrt(image_number)  # update the RDNOISE
                        gain = combined.meta["GAIN"]
                        combined.meta["GAIN"] = gain*image_number  # update the GAIN
                    combined.write(save_location / f"Master_{object_name}_{filter}.fits", overwrite = True)
                    print(f"{object_name} in {filter} filter combined!")
                    print("----------------------------------------------------------\n")

        return
