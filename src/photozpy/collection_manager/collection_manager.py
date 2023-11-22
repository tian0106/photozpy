"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Combine a series of registered images.
- Check image types before combine.
"""

from ccdproc import ImageFileCollection
from pathlib import Path
from itertools import product

class CollectionManager():

    def __init__(self):

        pass

    @staticmethod
    def refresh_collection(image_collection, rescan = False):
        """
        Refresh the image collection by reading the files. You can also rescan all the fits files to include newly created fits files.

        Paremeters
        ----------
        image_collection: ccdproc.ImageFileCollection; the image collection you want to refresh
        rescan: boolean; set to True if you want to rescan all the fits files in the location of image_collection to include new fits files.

        Returns
        -------
        new_image_collection: ccdproc.ImageFileCollection; the refreshed image collection
        """

        if rescan == False:
            new_image_collection = ImageFileCollection(location = Path(image_collection.location), 
                                                   filenames = image_collection.files)
        elif rescan == True:
            new_image_collection = ImageFileCollection(location = Path(image_collection.location), 
                                                       glob_include = "*.fits")

        return new_image_collection


    @staticmethod
    def unwarp_dictionary(dictionary):
    
        """
        Unwarp a dictionary to produce a dictionary list whose element is a dictionary that contains all the combinations of 
        values for the headers.
    
        Parameters
        ----------
        dictionary: dict; the input dictionary
    
        Returns
        -------
        dict_list: list; the list of all dictionaries that exhausts all the combination of the values.
        """
        
        # get all the combinations of values from different headers
        values = list(product(*list(dictionary.values())))
    
        headers = list(dictionary.keys())
    
        dict_list = []
        for i in values:
            dict_ = dict(zip(headers, i))  # assemble headers and values. Note that headers are always the same. We just iterate through the combination of the values.
            dict_list += [dict_]
    
        return dict_list
    
    @staticmethod
    def filter_collection(image_collection, **headers_values):
        """
        Filter the image collection based on the headers and values. One header can have multiple corresponding values.
    
        Paremeters
        ----------
        image_collection:
        headers_values: dict; the header names and corresponding values. The values should be in list.
                              {"header": ["value"]} and the list can contain multiple values.
    
        Returns
        -------
        new_image_collection
        """
    
        dict_list = CollectionManager.unwarp_dictionary(headers_values)
    
        files = []
        for dict in dict_list:
            files_temp = list(image_collection.files_filtered(**dict))
            files += files_temp
    
        # remove duplicate file names
        files = [*set(files)]
    
        new_image_collection = ImageFileCollection(location = Path(image_collection.location), filenames = files)
    
        return new_image_collection