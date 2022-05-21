#! /usr/bin/env python3
import numpy as np
import pandas as pd
# import os
import matplotlib.pyplot as plt
import fiona
import matplotlib.colors as colors
# from matplotlib import cm
from matplotlib_scalebar.scalebar import ScaleBar
import geopandas as gpd
import rasterio
import earthpy.spatial as es
import richdem as rd


def getdata(f):
    with rasterio.open(f) as src:
        el = src.read(1)
        BBox = (src.bounds[0],  src.bounds[2], src.bounds[1],  src.bounds[3])
        # Set masked values to np.nan
        el[el < 0] = np.nan
        hs = es.hillshade(el, altitude=10)

        gt = src.transform
        #PIXELSIZE
        piX = gt[0]
        piY = gt[4]
        #print(piX, piY)

    return(el, hs, BBox)


def plotHS(demsHEK):
    # set bounding box
    BBox_s = (652750, 653250, 5188750, 5189400)
    # set some stuff for plot
    al = 0.9
    fig, ax = plt.subplots(3, 4, figsize=(8, 8), sharex=True, sharey=True)
    ax = ax.flatten()
    # set DEMs to be plotted
    interesting_keys = [1953, 1969, 1971, 1977, 1990, 1997, 2006, 2010, 2018, 2019, 2020, 2021]
    subdict = {x: demsHEK[x] for x in interesting_keys if x in demsHEK}
    # Open the DEM with Rasterio, make and plot hillshades with earthpy
    for i, d in enumerate(subdict):
        print(d)
        el, hs, BBox = getdata(demsHEK[d])
        ax[i].imshow(hs, cmap="Greys", alpha=al, extent=BBox)
        ax[i].set_title(str(d))

    ax[0].set_xlim([652750, 653250])
    ax[0].set_ylim([5188750, 5189400])

    for i, a in enumerate(ax):
        ax[i].annotate('A', xy=(652900, 5188900),  xycoords='data',
                       fontsize=15, horizontalalignment='right', verticalalignment='top',)
        ax[i].annotate('B', xy=(652850, 5189100),  xycoords='data',
                       fontsize=15, horizontalalignment='right', verticalalignment='top',)

    scalebar = ScaleBar(1) # 1 pixel = 1 meter
    ax[11].add_artist(scalebar)
    fig.savefig('figs/hillshades.png', dpi=150)


def overview(Lines):
    #set crs of lines to match DEMs
    for L in Lines:
        Lines[L]['gdf'].to_crs(epsg=32632, inplace=True)

    # load DEMs:
    el, hs, BBox = getdata('/Users/leahartl/Desktop/HEK/aligned-resam_2017.tif')
    dem = rd.LoadGDAL('/Users/leahartl/Desktop/HEK/aligned-resam_2017.tif')
    slope = rd.TerrainAttribute(dem, attrib='slope_degrees')

    # load centerline:
    ln = gpd.read_file('data/centerline_1.shp')  
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))#, sharex=True, sharey=True)

    bounds = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    s = ax.imshow(slope, cmap="viridis", norm=norm, alpha=0.9, extent=BBox)
    ax.imshow(hs, cmap="Greys", alpha=0.2, extent=BBox)
    cbar = fig.colorbar(s, label='slope') #Adding colorbar with label
    scalebar = ScaleBar(1)  # 1 pixel = 1 meter
    ax.add_artist(scalebar)

    Lines['L0']['gdf'][Lines['L0']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=10, label='L0', color='k')
    Lines['L1']['gdf'][Lines['L1']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=10, label='L1', color='k')
    Lines['L2']['gdf'][Lines['L2']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=10, label='L2', color='k')
    Lines['L3']['gdf'][Lines['L3']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=10, label='L3', color='k')
    
    # comment out this line to remove longditudinal profile
    Lines['L4']['gdf'][Lines['L4']['gdf']['Date']>'2015-01-01'].plot(ax=ax, markersize=5, label='L4', color='k')

    ln.plot(ax=ax, color='k', linewidth=0.5)
    ax.annotate('P0', xy=(652816, 5189154),  xycoords='data',
                fontsize=15,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P1', xy=(652839, 5188913),  xycoords='data',
                fontsize=15,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P2', xy=(652985, 5188664),  xycoords='data',
                fontsize=15,
                horizontalalignment='right', verticalalignment='top',)
    ax.annotate('P3', xy=(653090, 5188504),  xycoords='data',
                fontsize=15,
                horizontalalignment='right', verticalalignment='top',)
    # ax.legend()
    fig.savefig('figs/overview.png', dpi=150)
    # plt.show()

