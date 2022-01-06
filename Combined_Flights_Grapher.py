# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 05:29:23 2021

@author: Mitch Creelman
"""
#%% Notes on Usage

# This script is optimized to compare 6 different DEM's side-by-side in order
# to give a quick idea of which ones to analyze together for more 
# in-depth analysis



# All files should be in .tif format
# file path is input on line 31
# flight names are input starting at line 38

#%% Levers and Adjustments

## Image options
# set dpi for images
dpi = 1500

color = 'winter' # colormap for difference map. make divergent



#%% File Names and Directories

# set the file path for where the rasters are located
file_path = "E:/Local Thesis Testing Ground/DEMs/"

# # Shape file for road extraction
# road = "Bogus_Road_20m Buffer_planar.shp"

# determine the flights (be sure to include the '.tif' ending)
flights = ["M20Jan.tif",
           "M04Feb.tif",
           "M17Feb.tif",
           "M24Feb.tif",
           "M03Mar.tif",
           "M10Mar.tif"] 



#%% Import Libraries

import fiona
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import rasterio as rio
import os
import rasterio.plot
from descartes import PolygonPatch

#%% Make Subfolder for Outputs

sub_dir = file_path + 'Output'

try:
    os.makedirs(sub_dir)
except OSError:
    print ("Did not create %s " % sub_dir)
else:
    print ("Successfully created the directory %s" % sub_dir)



#%% Visualize both days Bogus Basin Resort Zoom


fig, axs = plt.subplots(2 ,3 ,sharex = True, sharey = True)
plt.subplots_adjust(wspace=0)
fig.suptitle('Printouts of all flights')


#start indexing
i = 0
j = 0
p = 0


for DEM in flights:
    # NOTE: This reads just the metadata into memory, not the whole file
    if p == 6:
        break
    print('Times up.... Let''s do this!')
        
    # open up the shape file using fiona

    for i in range(0,2):
        
        for j in range (0,3):
            # We can read this raster data as a numpy array to perform calculations
            path = file_path+flights[p] # set the file path to be opened

            with rio.open(path) as src:
                data = src.read(1) #read first band
                
                plt.rcParams['figure.dpi'] = dpi # increase the dpi resolution
                axs[i, j].imshow(data,
                                 cmap = color) # show image of data

         
                scalebar = ScaleBar(1, # 1 pixel = 1 meter
                                    location = 'lower right', # set location in upper left
                                    border_pad=0.4) # add a bit of cushion to make it pretty
                axs[i, j].add_artist(scalebar) # add the scalebar to the graph
        
                axs[i, j].set_title(flights[p][:-4]) #set the plot title
                
                axs[i, j].get_xaxis().set_visible(False)
                axs[i, j].get_yaxis().set_visible(False)
                axs[i, j].set_facecolor('k')
                p = p + 1
                
                # north arrow
                x, y, arrow_length = 0.1, 0.9 , 0.1 # x position, y position, length of arrow
                axs[i, j].annotate('N', xy=(x, y), color = 'white',
                xytext=(x, y-arrow_length), arrowprops=dict(facecolor='white', width=5, headwidth=15),
                            ha='center', va='top', fontsize=10, xycoords=axs[i, j].transAxes)
                

    # add colorbar and other overall figure settings 
    axs = plt.gca()
    cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=color ),
                        cax=cbar_ax)   
    cbar.set_label('Elevation (m)', 
                    rotation=270, 
                    labelpad = 20) # colorbar label   


    plt.show()
    
#%% To-Do

# automate colorbar tick marks to max/min values of overall graphs
# add easting/northing labes/ranges somewhere (unsure how to make this look good)
# Add road overlay onto graphs
# make the top font bigger
# figure out how to reduce value between graphs on each row
# change scale bar to black background and white text
# move north arrow farther up

