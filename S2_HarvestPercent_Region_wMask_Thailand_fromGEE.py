# This program provides percentage of sugarcane harvested during the harvest season based on Sentinel-2 results from Google Earth Engine.
# It reads in a Sentinel-2 classification image (GeoTIFF) which includes a sugarcane class, 
# and two-weekly NDVI images for the sugarcane harvest season, and NDVI < 0.42 is considered to be harvested.
# The results are output to a .txt file in the working_directory (defined below), named SugarcaneExtentHarvested_(site)_(year).txt
# The masked NDVI images are output as GeoTIFF files to the classification_directory

# Note that the classification and NDVI images are assumed to be the same extent

#### Written in June 2022 ####

from osgeo import gdal
from osgeo import osr, ogr
import numpy as np
import os
import sys
import datetime
import xarray as xr
from os import listdir
from os.path import isfile, join

#################################################################################################
# Provide input information to run program

srs = osr.SpatialReference() # establish encoding
srs.ImportFromEPSG(4326) # WGS'84

working_directory = "Z:\\work\\cjt\\SugarCast\\Thailand\\MultiYear_CFNs\\"

# Region name  - which is 'Kanchanaburi', 'UdonThani' or 'NakhonRatchasima'
site = 'UdonThani' 
if site == 'Kanchanaburi': siteAbrev = 'kan'
if site == 'UdonThani': siteAbrev = 'udon'
if site == 'NakhonRatchasima': siteAbrev = 'nr'

classification_directory = working_directory + site + "\\Classifications\\"
NDVI_directory = working_directory + site + "\\NDVI\\" # NDVI images are in NDVI_directory + year folder + GEE_NDVI

# Specify sugarcane value in classification image
Sugar_class = 4

###############################################################################################

years = ['19-20', '20-21', '21-22']

for year in years:

	# Define input and output details
	year_txt = year[:2]+'to'+year[3:] #'19to20' '20to21' '21to22'
	Image_directory = NDVI_directory + "20"+year+"\\GEE_NDVI\\"
	mask_filename =  classification_directory + siteAbrev + "ML_"+year+".tif"
	outfile = open(working_directory + "SugarcaneExtentHarvested_" + site +"_20"+year_txt+".txt","w")

	# Open classification image
	src_ms = gdal.Open(mask_filename, gdal.GA_ReadOnly)
	gt = src_ms.GetGeoTransform()
	msband = src_ms.GetRasterBand(1)
	mask_band = msband.ReadAsArray()
	yt,xt = mask_band.shape
	msband = None

	# Calculate sugarcane extent
	Harvest_mask = np.zeros_like(mask_band, dtype=float)
	Harvest_mask[mask_band==4]=1 # '4' is sugarcane class for classification
	masked_band = np.ma.masked_where(Harvest_mask==0,Harvest_mask)
	TotalSugarcane_pixels = masked_band.count()
	print('Sugarcane_pixels in NovDec =',TotalSugarcane_pixels)
	outfile.write("Full sugarcane extent (pixels):"+str(TotalSugarcane_pixels)+' \n')
	outfile.write("Year,Month,Day, SugarcaneExtent (pixels), ExtentHarvested (%), % valid data \n")

	# Read in images and sort according to date
	fnames = [f for f in listdir(Image_directory) if isfile(join(Image_directory, f))]
	fnames = sorted(fnames)
	print(fnames)

	for f in fnames:
		if f.endswith(".tif"):
			print ("Reading NDVI file:", f)
			src_ds = gdal.Open(Image_directory+f, gdal.GA_ReadOnly)
			gt = src_ds.GetGeoTransform()
			srcband = src_ds.GetRasterBand(1)
			NDVI_band = srcband.ReadAsArray()
			yt_test, xt_test = NDVI_band.shape
			if ((yt != yt_test) and (xt != xt_test)): 
				print('Classification image and NDVI image are not the same size!')
				quit()
			srcband = None

			# Calculate day of year for harvest date output
			Year = f[5:9] #f[len(f)-19:len(f)-15]
			Month = f[9:11] #f[len(f)-15:len(f)-13]
			Day = f[11:13] #f[len(f)-13:len(f)-11]
			print('YMD=',Year,Month,Day)

			# Calculate % valid pixels in region
			NDVI_mask = np.ones_like(NDVI_band, dtype=float) 
			Good_pixels = np.ma.masked_where((np.isnan(NDVI_band)) | (mask_band==0),NDVI_mask)
			Num_Good = Good_pixels.count()
			Region_pixels = np.ma.masked_where(mask_band==0,NDVI_mask)
			Num_Region = Region_pixels.count()
			Percent_valid = int(100.0*(Num_Good/Num_Region))
			print('% valid in current image=', Percent_valid)

			# Create NDVI mask from NDVI image
			NDVI_mask[NDVI_band<0.42] = 0 # Find harvested area
			NDVI_mask[Harvest_mask==0] = 0 # Remove pixels outside remaining sugarcane area
			Harvest_mask[NDVI_mask==0] = 0 # Remove harvested pixels from sugarcane area
			print('Harvest_mask shape=',Harvest_mask.shape)
			masked_band = np.ma.masked_where(Harvest_mask==0,Harvest_mask)
			Sugarcane_pixels = masked_band.count()
			Harvested_percent = round(100.0 - (100.0*Sugarcane_pixels/TotalSugarcane_pixels),2)
			print('Remaining sugarcane pixels and % harvested from YMD =',Sugarcane_pixels, Harvested_percent, Year, Month, Day)
			outfile.write(str(Year)+','+str(Month)+','+str(Day)+','+str(Sugarcane_pixels)+','+str(Harvested_percent)+','+str(Percent_valid)+' \n')

			print('Writing Harvest_mask file')
			# Output NDVI mask image to GeoTIFF
			Image_Out = classification_directory + "SugarcaneExtent_"+str(Year)+str(Month)+str(Day)+".tif"
			ndwi_dsM = gdal.GetDriverByName('GTiff').Create(Image_Out, xt, yt, 1, gdal.GDT_Byte)
			ndwi_dsM.SetGeoTransform(gt) # specify coordinates
			ndwi_dsM.SetProjection(srs.ExportToWkt()) # export coords to file
			ndwi_dsM.GetRasterBand(1).WriteArray(Harvest_mask) # write band to raster
			ndwi_dsM.FlushCache()  # write to file
			ndwi_dsM = None # save and close

	outfile.close()
print('Finished!')