# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 15:06:59 2021

@author: 
    Mitch Creelman
    Boise State University - CryoGARS Group
"""

#%% Readme:
    
# this script is meant to take 2 separate DEM's and a pre=defined shape file 
    # outlining an intersecting road in order to determine the changes in height 
    # of snow for the entire area and then use the differences along the defined 
    # road to determin systematic offsets
    
#%% Steps to Use:

    # input files
        ## Have 2 separate .tiff files from 2 different observation days
            # type name of flights into lines 34-35 with .tif ending
                # if desired, you can enter as many flight names as you want, 
                # switching between them using lines 64 and 66
        ## have a .shp and .shx file for the road
            # type the name of the road into line 30 with .shp ending
        ## copy the file path where these are located to file_path on line 25
            # ensure it ends with a /
            
    # levers and adjustments
        # go through and look at how you want to fine-tune the comaprison. 
        # These should be fairly self-explanatory save the AOI
            # input northing/easting bounds for area of interest that you would 
            # like to zoom into for a better visualization of the comparison
                # adjusted at lines 80-83
                
    # hit run and watch the magic happen!
            
    
#%% File Names and Directories

# set the file path for where the rasters are located
file_path = "E:/Finals/Rasters/"

# Shape file for road extraction
road = "Bogus_Road_20m Buffer_planar.shp"

# determine the flights (be sure to include the '.tif' ending)
flights = ["0303_U2.tif",
           "0310_U2.tif"] 

    
#%% Levers and Adjustments

color = 'coolwarm' # colormap for difference map. make divergent

## Image options
# set dpi for images
dpi = 1300
road_width = 0.5 # width of road overlay line on rasters
road_color = 'w' # color of road overlay on rasters


# set range of outliers and inliers
mx = 3 # maximum value2
mi = -3 # minimum value

# set the coordinate reference system
ref_sys = "EPSG:32611"

# check for atribute offset. choices: 'aspect' 'slope_riserun' 'profile_curvature' 'planform_curvature' 'curvature'
# info on atributes: https://richdem.readthedocs.io/en/latest/terrain_attributes.html
atributes = ['aspect']

# Flights
# pick which flight you want to be the earlier
first = 0
# pick the later flight to compare it to 
second = 1

# AOI - northing/easting coordinates for area of interest
# headwall
ymin = 4842400 
ymax = 4844000 
xmin = 571800  
xmax = 573800



# set the file name for later
file_name = [flights[first],
             flights[second]]



#%% Import Libraries

import fiona # shape file reader
import geopandas as gpd # organize datasets into matrices for analysis
import matplotlib as mpl
import matplotlib.pyplot as plt # all the maths
from matplotlib_scalebar.scalebar import ScaleBar # adds the scalebar
import numpy as np
import os
import pandas as pd
import rasterio as rio
from rasterio.mask import mask
import rasterio.plot
from rasterio import Affine # or from affine import Affine
import richdem as rd # pull this from https://github.com/r-barnes/richdem
import rioxarray
from shapely.geometry import Point
from sklearn.metrics import r2_score
from descartes import PolygonPatch



#%% Make Subfolder for Outputs

sub_dir = file_path + 'Output'

try:
    os.makedirs(sub_dir)
except OSError:
    print ("Did not create %s " % sub_dir)
else:
    print ("Successfully created the directory %s" % sub_dir)



#%% Print Raster Function: 
  
    # found on https://corteva.github.io/rioxarray/stable/examples/reproject_match.html

def print_raster(raster):
    print(
        f"shape: {raster.rio.shape}\n"
        f"resolution: {raster.rio.resolution()}\n"
        f"bounds: {raster.rio.bounds()}\n"
        f"sum: {raster.sum().item()}\n"
        f"CRS: {raster.rio.crs}\n"
    )
    
    
    
#%% Checker Function

def limit_checker(array):
    
    nd = array.flatten() # flatten eerything into an array
    nd = np.unique(nd[~np.isnan(nd)]) # remove the nan's 
    print("Max Value: "+ str(max(nd))) # what is the max value? 
    print("Min Value: "+ str(min(nd))) # what is the min value? 
    
    
    
#%% Load in Rasters for Arrays

print(  )
print('Loading the Rasters for Both Arrays')
i = 0 # indexing start
for name in file_name:
    print(name)
    # Load the fsr dem using rioxarray / rasterio commands. 
    DEM = rioxarray.open_rasterio(file_path + name)#fn, masked=True, default_name=name, long_name=name).squeeze(dim='band', drop=True)
    
    # print out metadaat
    print('Metadata for '+ name)
    print(DEM)
    print_raster(DEM)
    print( )
    
    # Change the fill values to nan    
    DEM = DEM.where(DEM != DEM._FillValue, np.nan)    
    
    # Reproject this data to a different CRS and save it
    # https://epsg.io/26912
    DEM = DEM.rio.write_crs(ref_sys) # writes a CRS for specified ref. system
    
    # assign labels based on iterations
    if i == 0:
        
        DEM1 = DEM # assign a new name to the DEM
        print("DEM 1 Info") 
        print_raster(DEM1) #print the raster info 
        
    if i == 1:
        DEM2 = DEM # assign a new name to the DEM
        print("DEM 2 Info")
        print_raster(DEM2)  #print the raster info 
    
    i = i+1 # increase the index
    
    
    

#%% Visualize both days

print( )
print('Visualizing both Arrays')
## uncomment these when ready to turn things into a for loop and do both DEMS
DEM_List = [DEM1, DEM2]

#start indexing
i = 0

for DEM in DEM_List:
    # NOTE: This reads just the metadata into memory, not the whole file
    path = file_path+file_name[i] # set the file path to be opened
        
    # open up the shape file using fiona
    with fiona.open(file_path + road, "r") as shapefile:
        features = [feature["geometry"] for feature in shapefile]
    
        # We can read this raster data as a numpy array to perform calculations
    with rasterio.open(path) as src:
        data = src.read(1) #read first band
        # print(type(data)) # check on the data type
        plt.rcParams['figure.dpi'] = dpi # increase the dpi resolution
        plt.imshow(data, cmap = 'winter') # show image of data
                        
        # color bar
        cbar = plt.colorbar() # implement colorbar      
        cbar.set_label('Elevation (m)', 
                       rotation=270, 
                       labelpad = 20) # colorbar labels
        plt.title(file_name[i] + '\n ') #set the plot title
        
        ax = plt.gca() # get the current axis
        ax.set_facecolor('k')
        
        # add the road
        # import the bounds of the road for graphing
        patches = [PolygonPatch(feature, # apply the geometry of the road
                                edgecolor=road_color,  # load the edge color
                                facecolor=road_color,# load the fill color for the road
                                linewidth=road_width # line width of borders
                                ) for feature in features]
        
        # add the road boundary o the graph
        ax.add_collection(mpl.collections.PatchCollection(patches, 
                                                          match_original=True))

        # scale bar
        scalebar = ScaleBar(1, # 1 pixel = 1 meter
                            location = 'upper left', # set location in upper left
                            border_pad=0.4) # add a bit of cushion to make it pretty
        plt.gca().add_artist(scalebar) # add the scalebar to the graph
        
        # north arrow
        x, y, arrow_length = 0.09, 0.87 , 0.1 # x position, y position, length of arrow
        ax.annotate('N', xy=(x, y), color = 'white',
        xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
                    ha='center', va='top', fontsize=10, xycoords=ax.transAxes)

        ax.set_xlabel('Easting')
        ax.set_ylabel('Northing')
        rio.plot.show(src, cmap = 'winter') # show the plot CONTRLES THE PLOT?!?!?!?!
    i = i + 1



#%% Reproject the DEMs

print("DEM 1 New Info") # label
print_raster(DEM1) # double check atributesin printout
print( )

# DEM2 = DEM2.rio.write_crs(ref_sys) # Set CRS for DEM 2
print("DEM 2 New Info") # label for what's happening
print_raster(DEM2) # double check attributes in printout

# allign both rasters
DEM2 = DEM2.rio.reproject_match(DEM1) # found at ~55 minutes in reference videos 



#%% Difference the DEM's

# subtract date 1 from date 2
Snow_Diff = DEM1-DEM2 

# set a path for the raster to be saved to
Snow_Diff_Path = sub_dir + '/Snow_Diff.tif'

# convert differencing to a raster# convert differencing to a raster
Snow_Diff.rio.to_raster(Snow_Diff_Path, # save to the snow diff path as Snow_Diff.tiff
                        driver= src.driver, # set the same driver as before
                        nodata = np.nan, # set a nodata value
                        mask_and_scale = DEM1.scale_factor) #https://www.cogeo.org



#%% Full Differencing


path = Snow_Diff_Path
with rasterio.open(path) as src:
    print(src.profile)
    
# open up the shape file using fiona
with fiona.open(file_path + road, "r") as shapefile:
    features = [feature["geometry"] for feature in shapefile]

    # We can read this raster data as a numpy array to perform calculations
with rasterio.open(path) as src:
    data = src.read(1) #read first band
    print(type(data)) # check on the data type
    
    # plotting
    plt.rcParams['figure.dpi'] = dpi # increase the dpi resolution
    plt.imshow(data, 
               cmap = color) # show image of data
    cbar = plt.colorbar(cmap = color) # implement colorbar      
    cbar.set_label('Elevation (m)', 
                   rotation=270, 
                   labelpad = 20) # colorbar labels
    
    plt.title('Unbound Differences Between Days for \n' + # set the title
                 str(file_name[0] + ' and ' + str(file_name[1]))) # label dates #set the plot title
    ax = plt.gca() # get the current axis
    
    ax.set_facecolor('k')

    # add the road
    # import the bounds of the road for graphing
    patches = [PolygonPatch(feature, 
                            edgecolor=road_color, 
                            facecolor=road_color, 
                            linewidth=road_width) for feature in features]
    # add the road boundary o the graph
    ax.add_collection(mpl.collections.PatchCollection(patches, 
                                                      match_original=True))
    
    
    # scale bar
    scalebar = ScaleBar(1, # 1 pixel = 1 meter
                        location = 'upper left', # set location in upper left
                        border_pad=0.4) # add a bit of cushion to make it pretty
    plt.gca().add_artist(scalebar) # add the scalebar to the graph
    
    # north arrow
    x, y, arrow_length = 0.09, 0.87 , 0.1 # x position, y position, length of arrow
    ax.annotate('N', xy=(x, y), color = 'white',
    xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
                ha='center', va='top', fontsize=10, xycoords=ax.transAxes)

    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    
    rio.plot.show(src, # plot the difference data
                  cmap = color) # choose the color of the map
    
    SD = data.copy() # assign snow difference data to this variable
                          # this appeases the gis gods....for what reason i know not




#%% Visualize the inliers

# filiter out of SD
SD[SD < mi] = np.nan # filter out numberrs less than -3
SD[SD > mx] = np.nan # filter out numbers greater than 3

# plot everything out
DEM_in = SD # This is where we set the DEM to be graphed
fig, ax = plt.subplots(figsize=(5, 5))
data = np.squeeze(DEM_in.data) # set the data
plt.rcParams['figure.dpi'] = dpi # set the figure dps

# use imshow so that we have something to map the colorbar to
image_hidden = ax.imshow(data, 
                         cmap=color)

 # add the road
 # import the bounds of the road for graphing
patches = [PolygonPatch(feature, 
                        edgecolor=road_color, 
                        facecolor=road_color, 
                        linewidth=road_width
                        ) for feature in features]
# add the road boundary o the graph
ax.add_collection(mpl.collections.PatchCollection(patches, 
                                                  match_original=True))


# plot on the same axis with rio.plot.show
image = rio.plot.show(data, # show the squeezed data
                       transform=src.transform, # SETS X AND Y VARIABLES 
                      ax=ax, # set the axis
                      cmap=color) # colormap

# add colorbar using the now hidden image
cbar = fig.colorbar(image_hidden, # set the image that it is attached to  
             ax=ax) # determine the axis for the colorbar
cbar.set_label('Snow Accumulation (m)', # set colorbar label
               rotation=270, # rotate it for vertical
               labelpad = 20) # colorbar labels
ax.set_xlabel('Easting') # x label
ax.set_ylabel('Northing') # y label
ax.set_facecolor('k') # set background color to blac
scalebar = ScaleBar(1, # 1 pixel = 1 meter
                    location = 'upper left', # set location in upper left
                    border_pad=0.4) # add a bit of cushion to make it pretty
plt.gca().add_artist(scalebar) # add the scalebar to the graph


# north arrow
x, y, arrow_length = 0.09, 0.87 , 0.1 # x position, y position, length of arrow
ax.annotate('N', xy=(x, y), color = 'white',
xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
            ha='center', va='top', fontsize=10, xycoords=ax.transAxes)


plt.title("Height of Snow Difference DEM Map for \n" + # set the title
             str(file_name[0] + ' and ' + str(file_name[1])), pad = 2) # label dates,pad=2) # Title Plt

# idea for saving figure
# plt.savefig('myfile.png', bbox_inches="tight") # IT FREAKING WORKS! maybe add output path?



#%% Plot Outlier Data
# SO = Snow_Diff.data.copy() # extract data from xarray use.copy() to create a different copy
SO = Snow_Diff.copy()

# filter out of SD
np.where(np.abs(SO)<=mx,np.nan,SO)

DEM_in = SO # This is where we set the DEM to be graphed
fig, ax = plt.subplots(figsize=(5, 5))
data = np.squeeze(DEM_in.data)

# use imshow so that we have something to map the colorbar to
image_hidden = ax.imshow(data, 
                         cmap = color)


# plot on the same axis with rio.plot.show
image = rio.plot.show(data, # show the squeezed data
                       transform=src.transform, # SETS X AND Y VARIABLES 
                      ax = ax, # set the axis
                      cmap = color) # colormap

 # add the road
 # import the bounds of the road for graphing
patches = [PolygonPatch(feature, 
                        edgecolor = road_color, 
                        facecolor = road_color, 
                        linewidth = road_width
                        ) for feature in features]
 # add the road boundary o the graph
ax.add_collection(mpl.collections.PatchCollection(patches,
                                                  match_original=True))
 
# add colorbar using the now hidden image
cbar = fig.colorbar(image_hidden, # set the image that it is attached to  
             ax=ax) # determine the axis for the colorbar
cbar.set_label('Snow Accumulation (m)', # set colorbar label
               rotation=270, # rotate it for vertical
               labelpad = 20) # colorbar labels
ax.set_xlabel('Easting')
ax.set_ylabel('Northing')
ax.set_facecolor('k')

scalebar = ScaleBar(1, # 1 pixel = 1 meter
                    location = 'upper left', # set location in upper left
                    border_pad=0.4) # add a bit of cushion to make it pretty
plt.gca().add_artist(scalebar) # add the scalebar to the graph

# north arrow
x, y, arrow_length = 0.09, 0.87 , 0.1 # x position, y position, length of arrow
ax.annotate('N', xy=(x, y), color = 'white',
xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
            ha='center', va='top', fontsize=10, xycoords=ax.transAxes)

plt.title('Height of Snow (Hs) Values outside of +/- ' + str(mx) + '\n for \n' + # set the title
             str(file_name[0] + ' and ' + str(file_name[1])),pad=2) # Title Plt


# idea for saving figure
# plt.savefig('myfile.png', bbox_inches="tight")

#%% Analysis of inlier and outlier data
print( )
print('Differenced Data Inlier and Outlier Printout: ')
print( )


#reshape data for analysis
data = Snow_Diff.data.flatten() # put all values into a 1d column for data analysis
data = np.unique(data[~np.isnan(data)]) # take out the nan values from the dataset

# total number of elements in the dataset
Nt = np.size(data) 

# tell me how many elements there are total in the dataset
print('There are ' + str(Nt) + ' elements in the entire dataset')
print( )


# inlier data
inliers = data[data >= mi ] # find the data that is greater than or equal to -2
inliers = inliers[inliers <= mx ]  # find the data that is less than or equal to 2
Nin = np.size(inliers) # Total number of inliers
Pin = np.around(Nin/Nt*100, decimals=2) # % of total dataset that is in inlier territory

print('There are ' + str(Nin) + ' elements between ' + str(mi) + ' and ' + str(mx) + ' in the dataset' )
print(str(Pin) + '% of the data is between ' + str(mi) + ' and ' + str(mx))
print( )

outliersd = data[data < mi ] # first gfrequent number
outliersu = data[data > mx] # second frequent number
outliers = np.concatenate((outliersd,outliersu))
Nout = np.size(outliers) # total number of outliers
Pout = np.around(Nout/Nt*100, decimals=2) # % of total dataset that is in inlier territory
print('There are ' + str(Nin) + ' elements outside the range of '+ str(mi) +' and '+ str(mx))
print(str(Pout) + '% of the data is outside the range of '+ str(mi) +' and '+ str(mx))

Tdata = [[inliers],[outliers]] # total data combined for visual processing


#%% Double check using richdem: pulled from https://richdem.readthedocs.io/en/latest/terrain_attributes.html

# Create the Richdem array for both DEMS

# DEM1
arr1 = rd.rdarray(DEM1, no_data = np.nan ) # create a richdem object
newarr1 = arr1.reshape(arr1.shape[1],arr1.shape[2]) # bring it down to 2D

# DEM2
arr2 = rd.rdarray(DEM2, no_data = np.nan ) # create a richdem object
newarr2 = arr2.reshape(arr2.shape[1],arr2.shape[2]) # bring it down to 2D


#%% Check for atribute offset

for atribute in atributes:
    # Creat an aspect calculation for each object
    at1 = rd.TerrainAttribute(newarr1, attrib = atribute) # aspect calculation
    at2 = rd.TerrainAttribute(newarr2, attrib = atribute) # aspect calculation
    
    # difference the aspects
    d_aspect = at1-at2
    
    # plot everything out
    DEM_in = d_aspect # This is where we set the DEM to be graphed
    title = ("Differences in " + atribute + " Between Raster for \n" + # set the title
                 str(file_name[0] + ' and ' + str(file_name[1]))) # label dates"
    c_label = 'Attribute Offset' # Set colorbar label
    
    # Plotting
    fig, ax = plt.subplots(figsize=(5, 5)) # set the figure 
    data = np.squeeze(DEM_in.data) # set the data
    plt.rcParams['figure.dpi'] = dpi # set the figure dps
    # use imshow so that we have something to map the colorbar to
    image_hidden = ax.imshow(data, 
                             cmap=color)


    # plot on the same axis with rio.plot.show
    image = rio.plot.show(data, # show the squeezed data
                           transform=src.transform, # SETS X AND Y VALUES 
                          ax=ax, # set the axis
                          cmap=color) # colormap
    # add the road
    # import the bounds of the road for graphing
    patches = [PolygonPatch(feature, 
                            edgecolor = road_color, 
                            facecolor = road_color, 
                            linewidth=road_width
                            ) for feature in features]
    # add the road boundary o the graph
    ax.add_collection(mpl.collections.PatchCollection(patches,
                                                      match_original=True))
    
    # add colorbar using the now hidden image
    cbar = fig.colorbar(image_hidden, # set the image that it is attached to  
                 ax=ax) # determine the axis for the colorbar
    cbar.set_label(c_label, # set colorbar label
                   rotation=270, # rotate it for vertical
                   labelpad = 20) # colorbar labels
    ax.set_xlabel('Easting')
    ax.set_ylabel('Northing')
    ax.set_facecolor('k')
    scalebar = ScaleBar(1, # 1 pixel = 1 meter
                        location = 'upper left', # set location in upper left
                        border_pad=0.4) # add a bit of cushion to make it pretty
    plt.gca().add_artist(scalebar) # add the scalebar to the graph
    
    # north arrow
    x, y, arrow_length = 0.09, 0.87 , 0.1 # x position, y position, length of arrow
    ax.annotate('N', xy=(x, y), color = 'white',
    xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
            ha='center', va='top', fontsize=10, xycoords=ax.transAxes)


    plt.title(title + '\n ',pad=2) # Title Plt
    
    # histogram and bar plot
    # set the new figure
    d_aspect = np.unique(d_aspect[~np.isnan(d_aspect)]) # remove the nan's 
    fig, axs = plt.subplots(2) # sets out the number of subplots
    
    
    # Plot the histogram.
    n, bins, patches = axs[0].hist(d_aspect, # attribute plotted 
                                   bins=300, # number of bins
                                   alpha=0.6)
    axs[0].set_title("Differences in " + atribute + " Between Rasters for \n" + # set the title
                 str(file_name[0] + ' and ' + str(file_name[1]))) # label dates") # set the image title
    axs[0].set_ylabel('Number of Occurances') # set the y axis label
    axs[0].set_xlabel('Range of Values') # set the x axis labe 
    
    #axs[0].set_title('Differenced Value between '+ str(mi) +' and '+ str(mx) +'  '  + str(Pin) + '% of data')
    
    # boxplot
    axs[1].boxplot(d_aspect) # plot everything in a boxplot
    axs[1].set_ylabel('Range of Values') # set the y label
    
    fig.tight_layout() # create a tight figure to fit both  figs
    fig.show() # command figure show or the next figure will be put within this one
    
    # idea for saving figure
    # plt.savefig(title, bbox_inches="tight") # IT FREAKING WORKS! maybe add output path?
    
#%% Plots
i = 0
for data in Tdata:
    fig, axs = plt.subplots(2) # sets out the number of subplots
    
    
    # Plot the histogram.
    
    axs[0].hist(data, bins=300, alpha=0.6) 
    
    axs[0].set_ylabel('Number of Occurances') # set the y axis label
    axs[0].set_xlabel('Range of Values (m)') # set the x axis labes
    if i == 0: # if it is the inliers
        axs[0].set_title('Differenced Value between '+ str(mi) +' and '+ str(mx) +'  '  + str(Pin) + '% of data for \n' + # set the title
                     str(file_name[0] + ' and ' + str(file_name[1]))) # label dates
        
    else: # otherwise it will be the outliers
 
        axs[0].set_title('Differenced Value Outliers: ' + str(Pout) + '% of data for \n' + # set the title
                     str(file_name[0] + ' and ' + str(file_name[1]))) # label dates
    
    
    # boxplot
    axs[1].boxplot(data)
    axs[1].set_ylabel('Range of Values') # set the x axis labes 
    # set titles for the boxplot
    if i == 0:
        axs[1].set_title('Differenced Values Between ' + str(mi) +' and '+ str(mx))
        
    else:
         axs[1].set_title('Differenced Value Outliers Boxplot for \n' + # set the title
                      str(file_name[0] + ' and ' + str(file_name[1]))) # label dates

    fig.tight_layout() # create a tight figure to fit both figs
    
    i = i + 1


#%% Building Dataframes 

## Saving as Masked Array

Road_File = [file_path+file_name[0], # First Raster
              file_path+file_name[1], # Second Raster
              Snow_Diff_Path] # Snow Difference Raster


i = 0

for path in Road_File:
    # open up the shape file using fiona
    with fiona.open(file_path + road, "r") as shapefile:
        features = [feature["geometry"] for feature in shapefile]
    
    
    with rasterio.open(path, "r+") as src:
        print(src.profile)
        print( )
        data = src.read(1) #read first band
        
        try: 
            src.nodata = data[0][0]
        except:
            print('nodata value already written')
        else: 
            print('nodata value added')
            
        out_img, out_transform = mask(src, # source data
                                      shapes = features, # polygon coordinates
                                      crop=True) # crop according to shape file
                                       
        out_mask = np.ma.masked_where(out_img<0,
                                      out_img)
        
        # if i == 2:
        #     fig, ax = plt.subplots(figsize=(8,12))
        #     outPlot = ax.imshow(np.squeeze(out_mask),
        #                         cmap='terrain')
            
        #     #copying metadata from original raster
        #     out_meta = src.meta.copy()
            
        #     #amending original metadata
        #     out_meta.update({'nodata': src.nodata,
        #                   'height' : out_img.shape[0],
        #                   'width' : out_img.shape[1],
        #                   'transform' : out_transform})
        # else:
            
            
        #writing and then re-reading the output data to see if it looks good 
        msk_data = out_mask.data
        
        #define no data value
        no_data=src.nodata
        # gdf = gpd.GeoDataFrame(data=msk_data)
    
        # extract the row, columns of the valid values
        elev = np.extract(msk_data != no_data, msk_data)  # extract the elevation data
        ind, row, col = np.where(out_mask != no_data)# extract the index data
            
        if i <= 1: # indexing and reducing elevation
            d = gpd.GeoDataFrame({ 'index':ind,
                              'col':col,
                              'row':row,
                              'Elevation (m)':elev})
            # show array size before and after to make sure we are eliminating NaN values
            print("Initial Size: " + str(d.size))
            d = d[~d['Elevation (m)'].isnull()]
            print("Reduced Size: " + str(d.size))

        else: # indexing and reducing height differences
            d = gpd.GeoDataFrame({ 'index':ind,
                              'col':col,
                              'row':row,
                              'Difference (m)':elev})  
            # show array size before and after to make sure we are eliminating NaN values
            print("Initial Size: " + str(d.size))
            d = d[~d['Difference (m)'].isnull()]
            print("Reduced Size: " + str(d.size))
            
        if i == 0:
          print( )
          print('First Road Data Loaded:')
          road1 = d.copy() # assign first road data
          print( )
          print(road1)
          print( )
        if i == 1:
          print( )
          print('Second Road Data Loaded:')
          road2 = d.copy() # assign second road data
          print( )
          print(road2)
          print( )
        if i == 2: 
          print( )
          print('Road Difference Data:')
          road_diff = d.copy() # assign difference data
          print( )
          print(road_diff)
          print( )
 
    i = i + 1 # INCREASE INDEX


## Merge GPD Arrays Based on row and column indices and Calculate Geometry

all_road = pd.merge(road1,road2, how = 'left')
print(all_road)
print( )


f = pd.merge(all_road,road_diff, how = 'left')
f = f[~f['Difference (m)'].isnull()]
print('Database Built')

print(f)
print( )

T1 = out_transform * Affine.translation(0.5, 0.5) # reference the pixel centre
rc2xy = lambda r, c: (c, r) * T1  
print( )

# coordinate transformation
f['x'] = f.apply(lambda row: rc2xy(row.row,row.col)[0], axis=1)
print('X coordinates Calculated')
print( )
f['y'] = f.apply(lambda row: rc2xy(row.row,row.col)[1], axis=1) 
print('Y coordinates Calculated')
print( )

# geometry
print( )
print('Calculating Geometry')
f['geometry'] = f.apply(lambda row: Point(row['x'], row['y']), axis=1)
print( )
print('Geometry Calculated')

print(f)


#%% Graphing Elevation vs. difference

f.plot.scatter(x = "Elevation (m)", # scatter plot with elevation as x axis
               y = 'Difference (m)', # set the y axis
               s = 0.3, # set the point size
               alpha = 0.075) # set the transparency

x = f['Elevation (m)'] # set x as the elevation
y = f['Difference (m)'] # set y as the differences in height

z = np.polyfit(x, y, 1) # fit z
y_hat = np.poly1d(z)(x) # fit the change in y
s =f.size # determine the size for the number of points

plt.plot(x, y_hat, "r--", lw=1.0) # plot the line over the scatter plot

# text outputs for linear equation, r2, and number of elements
text = f" $y={z[0]:0.3f}\;x{z[1]:+0.3f}$ \n $R^2 = {r2_score(y,y_hat):0.3f}$ \n $n = {s}$"

plt.gca().text(0.05, 0.95, 
               text,transform=plt.gca().transAxes, 
               fontsize=10,  # set font size
               verticalalignment='top') # set alignment

# set axis and titles
ax = plt.gca() # grab current axis
ax.set_title('Road Differences vs. Elevation for \n' + # set the title
             str(file_name[0] + ' and ' + str(file_name[1]))) # label dates
plt.tight_layout() # set the layout

plt.show() # show the plot


#%% Notes
# Raster Alignment Video: https://www.youtube.com/watch?v=ZXg70BH3wW4&
# Kernel Crash Fix: https://stackoverflow.com/questions/69786885/after-conda-update-python-kernel-crashes-when-matplotlib-is-used

# unsued but awesome tutoria overlaying road:
    # https://www.earthdatascience.org/courses/scientists-guide-to-plotting-data-in-python/plot-spatial-data/customize-raster-plots/plotting-extents/

# try/except/else statement for getting around errors
#    https://stackoverflow.com/questions/29065179/python-an-if-statement-if-an-exception-occurs


#%% To-Do List

# incorporate GCP's:
# https://rasterio.groups.io/g/dev/topic/use_of_gcps_in_rasterio/32990267?p=


# Necessary switches
# Saving or not saving figures 
