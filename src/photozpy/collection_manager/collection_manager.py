"""
Written by Yong Sheng at Clemson University, 2023 for the photozpy project.
Advisor: Dr. Marco Ajello
Other contributor(s):

Main function: 
- Combine a series of registered images.
- Check image types before combine.
"""

from ccdproc import ImageFileCollection

class CollectionManager():

    def __init__(self):

        pass


    # @staticmethod
    # def filter_collection(image_collection, **headers_values):
    #     """
    #     Filter the image collection based on the headers and values. One header can only have one corresponding value.

    #     Paremeters
    #     ----------
    #     image_collection:
    #     headers_values: dict; the headers their values. 

    #     Returns
    #     -------
    #     new_image_collection
    #     """
        
    #     location = image_collection.location
    #     files = list(image_collection.files_filtered(**headers_values))
    
    #     new_image_collection = ImageFileCollection(location = location, filenames = files)
    
    #     return new_image_collection


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
            new_image_collection = ImageFileCollection(location = image_collection.location, 
                                                   filenames = image_collection.files)
        elif rescan == True:
            new_image_collection = ImageFileCollection(location = image_collection.location, 
                                                       glob_include = "*.fits")

        return new_image_collection


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

        files = []
        for header, values in headers_values.items():
            for value in tuple(values):
                files_temp = list(image_collection.files_filtered(**{header: value}))
                files += files_temp

        new_image_collection = ImageFileCollection(location = image_collection.location, filenames = files)
    
        return new_image_collection