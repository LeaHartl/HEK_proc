#! /usr/bin/env python3
import numpy as np
import pandas as pd
# import xarray as xr
# import os
import matplotlib.pyplot as plt
from datetime import datetime
import fiona
import geopandas as gpd
from shapely.geometry import Point
# import files with code to make figures:
import figures as fg
import rasterPlots as rp

# turn off setting with copy warning - CAREFUL!
pd.options.mode.chained_assignment = None  # default='warn'

# ---some helper functions----
# read csv files exported from martin's excel file. need to replace , with . in csv. decimal is inconsistent! renamed Stones
# in L1 after line was moved in 2008 to retain unique Stone IDs
def readFiles(fname, L):
    data = pd.read_csv(fname,  delimiter=';', header=1)
    data = data[['Date', 'Stone', 'North [m]', 'East [m]', 'Elevation [m]']]
    data.columns = ['Date', 'Stone', 'x', 'y', 'z']
    data['dxyz'] = np.nan
    data['dxyz_a'] = np.nan
    data['dxy'] = np.nan
    data['dz'] = np.nan
    data['Line'] = L
    data['Date'] = pd.to_datetime(data['Date'])

    for s in data.Stone.unique():
        # print(s)
        # print(data.loc[data['Stone'] == s, 'z'])
        d = data.loc[data['Stone'] == s]
        d['dz'] = d['z'].diff()
        data.loc[data['Stone'] == s, 'dz'] = d['dz'].values

        d['dxy'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2))
        data.loc[data['Stone'] == s, 'dxy'] = d['dxy'].values

        d['dxyz'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2 + d['z'].diff()**2))
        data.loc[data['Stone'] == s, 'dxyz'] = d['dxyz'].values

        # print(d['Date'].diff().astype(int))
        d['dxyz_a'] = (d['dxyz'] / d['Date'].diff().dt.days) * 365
        data.loc[data['Stone'] == s, 'dxyz_a'] = d['dxyz_a'].values

    geometry = [Point(xy) for xy in zip(data.y, data.x)]
    gdf = gpd.GeoDataFrame(data, crs="EPSG:31254", geometry=geometry)
    return(data, gdf)


# compute annual displacement
def annualavg(val, date1, date2):
    delta = datetime.strptime(date2, "%Y-%m-%d")-datetime.strptime(date1, "%Y-%m-%d")
    val2 = (val/abs(delta.days))*365
    return(val2)


# convert to decimal degrees
def convertDMS(DMS):
    deg = [i.split('°')[0] for i in DMS]
    minutes = [i.partition('°')[2].partition('′')[0] for i in DMS]
    seconds = [i.partition('′')[2].partition('′′')[0] for i in DMS]
    deg = np.array(deg, dtype=np.float32)
    minutes = np.array(minutes, dtype=np.float32)
    seconds = np.array(seconds, dtype=np.float32)

    dec = (deg.astype(float) + minutes.astype(float)/60 + seconds.astype(float)/(60*60)) #* (-1 if direction in ['W', 'S'] else 1)
    return(dec)


# make shapefile with all stone coordinates
def make_shp(Lines):
    #make one shapefile of stone coordinates
    L01 = Lines['L0']['gdf'].append(Lines['L1']['gdf'])
    L012 = L01.append(Lines['L2']['gdf'])
    L0123 = L012.append(Lines['L3']['gdf'])
    L01234 = L0123.append(Lines['L4']['gdf'])
    L01234.dropna(axis=0, subset=['Date'], inplace=True)
    L01234.Date = L01234.Date.astype(str)

    L01234.to_crs(epsg=32632, inplace=True)
    L01234.to_file('data/Steine_32632.shp')
    L01234.to_file('data/Steine_32632.geojson', driver='GeoJSON')

# ----------------------------

# load shapefile of RG outline
bg = gpd.read_file('data/HEK2018.shp')  
# set crs to match stone coordinates
bg.to_crs(epsg=31254, inplace=True)

# initiatlize dictionary of line data
Lines = {
        'L0': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 0-Table 1.csv',
        },
        'L1': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 1-Table 1.csv',
                    },
        'L2': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 2-Table 1.csv',
                    },
        'L3': {
            'fname': 'data/HEK_Bewegung_seit1996/Linie 3-Table 1.csv',
        },
        'L4': {
            'fname': 'data/HEK_Bewegung_seit1996/Längs-Table 1.csv',
        }
        }


# initiatlize DF for mean velocities (only for cross profiles)
meanv = pd.DataFrame(index=np.arange(1997, 2022), columns=['L0', 'L1', 'L2', 'L3'])


demsHEK = {1953: 'data/DEMs/transf_32632_1953.tif',
           1969: 'data/DEMs/transf_32632_1969.tif',
           1971: 'data/DEMs/transf_32632_1971.tif',
           1977: 'data/DEMs/transf_32632_1977.tif',
           1990: 'data/DEMs/transf_32632_1990_geo.tif',
           1997: 'data/DEMs/transf_32632_1997_geo.tif',
           2006: 'data/DEMs/aligned-resam_2006.tif',
           2009: 'data/DEMs/aligned-2009_test_geo.tif',
           2010: 'data/DEMs/aligned-2010_test_geo.tif',
           2017: 'data/DEMs/aligned-resam_2017.tif',
           2018: 'data/DEMs/aligned-resam2018.tif',
           2019: 'data/DEMs/aligned-resam2019.tif',
           2020: 'data/DEMs/aligned-resam2020.tif',
           2021: 'data/DEMs/aligned-resam2021.tif'}



# add gdf to Lines dictionary and populate DF of mean vel.
for L in Lines:
    # print(L)
    dat, Lines[L]['gdf'] = readFiles(Lines[L]['fname'], L)
    dat['year'] = dat['Date'].dt.year
    Lines[L]['dat'] = dat
    # Lines[L]['Line'] = L

    for y in meanv.index:
        meanv.loc[y, L] = dat.loc[dat['year'] == y, 'dxyz_a'].mean()


# deal with reference points
dat = pd.read_csv('data/fixpunkte.csv')
dat['decx'] = convertDMS(dat['startx'].astype(str).values)
dat['decy'] = convertDMS(dat['starty'].astype(str).values)

dat = dat[['p', 'decx', 'decy']]
gdf = gpd.GeoDataFrame(
    dat, geometry=gpd.points_from_xy(dat['decx'], dat['decy']))
gdf.set_crs(epsg=4326, inplace=True)
gdf.to_crs(epsg=31254, inplace=True)
gdf['x'] = gdf.geometry.y
gdf['y'] = gdf.geometry.x

# ---------make shapefile------------
# make_shp(Lines)

# ---------plot some stuff-----------

# flowline plots (surface elevation along flowline)
fg.FlowlineAll(demsHEK)
fg.FlowlineSubplots(demsHEK)

# subplots of single blocks (2016-2021)
fg.SingleBlocks(Lines, gdf)
# simple map to show where the lines are
fg.figMap(Lines, bg)
# basic 3d plot of longditudinal profile (L4) and L0
fg.plot3d(Lines['L4']['gdf'], Lines['L0']['gdf'])
# basic time series of mean profil velocity
fg.figTS(meanv)
# test plot to show tmeseries of single blocks
fg.figTS_Stones(Lines)
# nicer time series plot, this contains a lot of hardcoded stuff, should fix
fg.timeseriesBoth()

# flowline plots:
fg.FlowlineAll(demsHEK)

# subplots of hillshades (selection of years is hard coded..)
rp.plotHS(demsHEK)
# overview plot
rp.overview(Lines)

plt.show()



