#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 15:00:21 2023

@author: yongsheng
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5
from swifttools.swift_too import Data, ObsQuery, TOORequests
import time
import shutil
import glob


class UVOTZ():
    
    def __init__(self, analysis_root_path, src_meta_file = None):
        
        
        """ The input is a cvs file with 3 colums seperated by a comma. 
            The column names are name, ra, and dec.
            You can generate a template by using the function 
            
            The analysis root path should contain the folders with sourcs names.
        """
        
        self.ana_root_dir = analysis_root_path
        self.full_filter_list = ["ubb", "um2", "uuu", "uvv", "uw1", "uw2"]
        
        if not src_meta_file:
            print("You haven't defined the src meta info yet!")
            df = pd.DataFrame(columns=["name", "ra", "dec"])
            df.to_csv(self.ana_root_dir+"/metadata.csv", sep = ",", index=False)
            print("A template csv file is generated in the analysis_root_path!")
            print("rNote that the csv is seperated by ,")
        else:
            self.src_meta_data = self.ana_root_dir + "/" + src_meta_file
            print(self.src_meta_data)
            
    def str_to_list(self, str):
        return str[1:-2].replace("'", "").split(", ")
        
    def keywords_finder(self, file_path, *keywords):
        ''' 
        This function takes a fits file and any keywords you want to query
        from the fits file. It deals with multi-extension fits file so you 
        will get a dictionary of all the keywords in each extension.
        
        The ext_No and EXTNAME are mandatory.
        
        # Run
        keywords_finder(file_path, "ASPCORR", "HDUCLAS1")
        
        # Output
        {'ext_No': [1, 2],
         'EXTNAME': ['bb649482961I', 'bb649505811I'],
         'ASPCORR': ['DIRECT', 'DIRECT'],
         'HDUCLAS1': ['IMAGE', 'IMAGE']}
        '''
        with fits.open(file_path) as hdul:
            ext_nums = len(hdul) - 1
            
            # initialize the dictoinary
            dict_ext = {}
            dict_ext["ext_No"] = [i for i in np.arange(1, ext_nums+1)]
            dict_ext["EXTNAME"] = [None] * ext_nums
            for i in keywords:
                dict_ext[i] = [None] * ext_nums
            
            # read and record the keywords in the dictionary
            n_ = 1
            while n_ <= ext_nums:
                dict_ext["EXTNAME"][n_-1] = hdul[n_].header["EXTNAME"]
                for i in keywords:
                    dict_ext[i][n_-1] = hdul[n_].header[i]
                n_ += 1
                
        return dict_ext
    
    def add_dict_df(self, dict_, df):
        for key in dict_.keys():
            if key not in df.columns:
                df[key] = "None"
            
        for i in np.arange(df.shape[0]):
            for key, value in zip(dict_.keys(), dict_.values()):
                df.loc[i, key] = str(value)
                
        return df
                   
    def obsquery_input(self, **kwargs):
        
        
        if "obsid" in kwargs.keys():
            return {"name":kwargs["name"],
                    "obsid": kwargs["obsid"]}
        
        elif "targetid" in kwargs.keys():
            return {"name":kwargs["Name"],
                    "targetid": kwargs["targetid"]}
        
        elif "ra" and "dec" in kwargs.keys():
            ra_dec = {key: value for key, value in kwargs.items() if key in {"ra", "dec"}}
            skycoord = SkyCoord(**ra_dec, unit = (u.hourangle, u.deg), frame = "icrs")
            
            return {"name":kwargs["name"],
                    "skycoord": skycoord}
        
    def ObsDownload(self, begin, end, radius = 5/60, uvotmode = "0x30ed"):
        
        # Read the csv file
        df = pd.read_csv(self.src_meta_data, sep = ",")
        
        row_nums = df.shape[0]
        
        obsid_all = []  # the obsid for all the sources
        obsid_path = []  # the obsid path for all the sources
        obsid_time = []  # the obsid time for all the sources
        files = []
        filters = []
        
        for i in np.arange(row_nums):
            
            info_ = dict(df.iloc[i])
            src_name = info_["name"].replace(" ", "_")
            src_ra = info_["ra"]
            src_dec = info_["dec"]
            src_path = self.ana_root_dir + "/" + src_name
            
            # create or delete the folder for the source
            if os.path.exists(src_path):
                shutil.rmtree(src_path)
            os.mkdir(src_path)
                 
            target_info = self.obsquery_input(**info_)
            oq = ObsQuery(radius = radius, begin = begin, end = end, **target_info)
            id_ = "00000000000"  # this is used to avoid downloading the same data file mutiple times by comparing the observation id
            
            src_obsid = []
            src_obsid_path = []
            src_obsid_time = []
            for i in np.flip(np.arange(-len(oq), 0)):
                if oq[i].obsid == id_:
                    print(f"{oq[i].obsid} has been downloaded/examined, skipping ...")
                else:
                    if oq[i].uvot_mode == uvotmode:
                        date_ = oq[i].begin.strftime("%Y-%m-%d %H:%M:%S")
                        print(f"{oq[i].obsid} on {date_} is being downloaded, the uvot mode is {oq[i].uvot_mode}")
                        oq[i].download(uvot = True, outdir = src_path)
                        
                        src_obsid.append(oq[i].obsid)  # record the obsid for this source
                        src_obsid_path.append(src_path + "/" + oq[i].obsid)  # record the path to obsids for this source
                        src_obsid_time.append(date_)
                        
                        id_ = oq[i].obsid
                    else:
                        print(f"{oq[i].obsid} uvod mode {oq[i].uvot_mode} is not the one you requested as {uvotmode}, skipping ...")
                        id_ = oq[i].obsid
                        
            obsid_all.append(src_obsid)
            obsid_path.append(src_obsid_path)
            obsid_time.append(src_obsid_time)
            print(f"Source {src_name} download finished!")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
        df["obsid_all"] = obsid_all
        df["obsid_path"] = obsid_path
        df["obsid_time"] = obsid_time
        df.to_csv(self.ana_root_dir+"/metadata.csv", sep = ",", index=False)
        
        return
    
    def Sum1(self, file, sum_type, **kwargs):
        if sum_type == "final":
            filter_name = file[-13:-10]
            outname = kwargs["summed_dir"] + f"/all{filter_name}.fits"
            os.system(f"uvotimsum infile={file} outfile={outname} | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            #print(f"the summ output is {outname}")
        elif sum_type == "intermediate":
            filter_name = file[-13:-10]
            obsid = file[-24:-13]
            outname = kwargs["summed_dir"] + f"/{filter_name}_{obsid}.fits"
            os.system(f"uvotimsum infile={file} outfile={outname} | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            #print(f"the summ output is {outname}")
        elif sum_type == "exclude=NONE":
            filter_name = file[-8:-5]
            outname = kwargs["summed_dir"] + f"/all{filter_name}.fits"
            os.system(f"uvotimsum exclude=NONE infile={file} outfile={outname} | tee -a uvotimsum_log.txt >/dev/null 2>&1")
            
        else:
            print("Wrong sum type! Please use final or intermediate")
            #break
        return outname

    def ASPCORR_check(self, file):
        dict_exts = self.keywords_finder(file, "ASPCORR")
        asps = dict_exts["ASPCORR"]
        check = ["Yes" for i in asps if "DIRECT" in i]
        if len(check) == len(asps):
            return True
        elif len(check) != len(asps):
            return False
        
    def del_create_dir(self, dir_path):
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.mkdir(dir_path)
        return
    
    def create_dir(self, dir_path):
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        return
    
    def src_region(self, save_dir, ra, dec):
        
        reg_files = [save_dir + "/" + i + ".reg" for i in self.full_filter_list]
        
        for reg in reg_files:
            f = open(reg, "w")
            f.write("# Region file format: DS9 version 4.1\n")
            f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
            f.write("fk5\n")
            f.write(f'circle({ra},{dec},5.000")')
            f.close()
        #return reg_files
    
    def UVOTSum(self):
        
        df = pd.read_csv(self.src_meta_data, sep = ",")
        df["final_fits"] = None
        df["reg_files"] = None
        
        for i in np.arange(df.shape[0]):
            final_fits = []
            obsids = self.str_to_list(df.loc[i,"obsid_all"])
            ra = df.loc[i,"ra"]
            dec = df.loc[i,"dec"]
            print(f"Now summing {df.iloc[i,0]}")
            
            # creat the summed dir
            summed_dir = self.ana_root_dir + "/" + df.loc[i, "name"].replace(" ", "_") + "/summed"
            self.del_create_dir(summed_dir)
            
            if len(obsids) == 1:
                
                obsid_path = self.str_to_list(df.loc[i,"obsid_path"])[0]
                img_files = glob.glob(obsid_path + "/uvot/image/" + "*sk.img.gz")
                for img_file in img_files:
                    if self.ASPCORR_check(img_file):
                        summed_fits = self.Sum1(img_file, sum_type = "final", summed_dir = summed_dir)
                        final_fits.append(summed_fits)
            else:
                for obsid_path in self.str_to_list(df.loc[i,"obsid_path"]):
                    
                    img_files = glob.glob(obsid_path + "/uvot/image/" + "*sk.img.gz")
                    for img_file in img_files:
                        if self.ASPCORR_check(img_file):
                            self.Sum1(img_file, sum_type = "intermediate", summed_dir = summed_dir)
                            
                inter_files = glob.glob(summed_dir + "/*fits")            
                for filter_ in self.full_filter_list:
                    filter_fits = [fits for fits in inter_files if filter_ in fits]
                    if len(filter_fits) == 1:
                        shutil.copy2(filter_fits[0], filter_fits[0][0:-21] + f"/all_{filter_}.fits")
                    else:
                        fappended_file = filter_fits[-1][0:-21] + f"/_{filter_}.fits"  # the file to be appended
                        shutil.copy2(filter_fits[-1], fappended_file)
                        
                        for j in filter_fits[0:-1]:
                            os.system(f"fappend {j} {fappended_file}")
                        
                        summed_fits = self.Sum1(fappended_file, sum_type="exclude=NONE", summed_dir = summed_dir)
                        final_fits.append(summed_fits)
                        #os.remove(fappended_file)
                        
                self.create_dir(summed_dir + "/intermediate")
                _ = summed_dir + "/*_*"
                inter_dir = summed_dir + "/intermediate"
                os.system(f"mv {_} {inter_dir}")
                
            self.src_region(save_dir = summed_dir, ra = ra, dec = dec)
            
            reg_files = [summed_dir + "/" + m[-8:-5] + ".reg" for m in final_fits]
            
            df.loc[i, "final_fits"] = str(final_fits)
            df.loc[i,"reg_files"] = str(reg_files)
            
        df.to_csv(self.ana_root_dir+"/metadata.csv", sep = ",", index=False)
        
    def find_mag(self, result):
        file_read = open(result, "r")
        text = "AB system"
        
        lines = file_read.readlines()
        ab_mag = 0
        ab_err = 0
        
        for i in np.arange(len(lines)):
            if text in lines[i]:
                ab_mag_line = lines[i+1]
                #print(ab_mag_line)
                
        file_read.close()
        
        if ab_mag_line[22:23] == ">":
            ab_mag = round(float(ab_mag_line[24:]),4)
            ab_err = -1
        else:
            ab_mag = round(float(ab_mag_line[22:27]),4)
            stat_err = float(ab_mag_line[32:36])
            sys_err = float(ab_mag_line[48:52])
            ab_err = round(np.sqrt(stat_err**2+sys_err**2),4)
            
        return ab_mag, ab_err
                
        
    def UVOTSource(self):
        
        df = pd.read_csv(self.src_meta_data, sep = ",")
        
        df_results = pd.DataFrame(columns = ['source_name', "uw2", "uw2_err", "um2", "um2_err", "uw1", "uw1_err", "uuu", "uuu_err", "ubb", "ubb_err", "uvv", "uvv_err"])
        
        for i in np.arange(df.shape[0]):
            
            src_name = df.loc[i, "name"].replace(" ", "_")
            dict_new = {'source_name' : src_name,
                        "uw2" : -99, 
                        "uw2_err" : -99,
                        "um2" : -99,
                        "um2_err" : -99,
                        "uw1" : -99,
                        "uw1_err" : -99,
                        "uuu" : -99,
                        "uuu_err" : -99,
                        "ubb" : -99,
                        "ubb_err" : -99,
                        "uvv" : -99,
                        "uvv_err" : -99}
            df_results = df_results.append(dict_new,ignore_index = True)
            
            for f, r in zip(self.str_to_list(df.loc[0,"final_fits"]), self.str_to_list(df.loc[0,"reg_files"])):

                # defining all kinds of inputs and outputs also the column names
                fits_file = f
                src_region = r
                summed_dir = r[:-8]
                filter_ = r[-7:-4]
                filter_err = r[-7:-4] + "_err"
                bkg_region = r[0:-7] + "bg" + r[-7:]
                outfits = summed_dir + "/" + r[-7:-4] + "_Results.fits"
                outtxt = summed_dir + "/" + r[-7:-4] + "_Results.txt"
                
                # run the command
                print(f"Running uvotsource image={fits_file} srcreg={src_region} bkgreg={bkg_region} sigma=3 cleanup=y clobber=y outfile={outfits} | tee {outtxt} >/dev/null")
                os.system(f"uvotsource image={fits_file} srcreg={src_region} bkgreg={bkg_region} sigma=3 cleanup=y clobber=y outfile={outfits} | tee {outtxt} >/dev/null")
                
                # extract the magnitudes
                ab_mag, ab_err = self.find_mag(outtxt)
                
                index = df_results.shape[0]-1
                
                df_results.loc[index, filter_] = ab_mag
                df_results.loc[index, filter_err] = ab_err
                
                df_results.to_csv(self.ana_root_dir+"/Magnitudes.csv", sep = ",", index=False)
                
                
                
                
            
            
        
        
        
                
        
    
            
                        
                        
        
        
        
    