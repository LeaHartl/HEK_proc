#! /usr/bin/env python3
import numpy as np
import pandas as pd
import fiona
import geopandas as gpd
import rasterio
import matplotlib.pyplot as plt
from matplotlib import cm


# plot shapefile of outline and stone positions
def figMap(Lines, bg):
    fig, ax = plt.subplots(figsize=(7, 7))
    bg.geometry.boundary.plot(color=None, edgecolor='k', linewidth=2, ax=ax, label='2018 outline')
    Lines['L0']['gdf'].plot(ax=ax, markersize=1, label='L0')
    Lines['L1']['gdf'].plot(ax=ax, markersize=1, label='L1')
    Lines['L2']['gdf'].plot(ax=ax, markersize=1, label='L2')
    Lines['L3']['gdf'].plot(ax=ax, markersize=1, label='L3')
    Lines['L4']['gdf'].plot(ax=ax, markersize=1, label='L_L')
    ax.legend()
    fig.savefig('figs/figMap.png', dpi=150)


# time series plot of mean vel. per line
def figTS(meanv):
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(meanv.index, meanv.L0, label='L0')
    ax.plot(meanv.index, meanv.L1, label='L1')
    ax.plot(meanv.index, meanv.L2, label='L2')
    ax.plot(meanv.index, meanv.L3, label='L3')
    ax.grid()
    ax.set_ylabel('velocity (xyz) in m/a')
    ax.legend()
    fig.savefig('figs/meanv.png', dpi=150)


# time series of single blocks
def figTS_Stones(Lines):
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    ax = ax.flatten()

    Lines['L1']['dat']['St'] = Lines['L1']['dat']['Stone'].str.replace(r'\D', '')
    Lines['L1']['dat'] = Lines['L1']['dat'].dropna()

    Lines['L0']['dat']['St'] = Lines['L0']['dat']['Stone']
    Lines['L0']['dat'] = Lines['L0']['dat'].dropna()

    Lines['L2']['dat']['St'] = Lines['L2']['dat']['Stone']
    Lines['L2']['dat'] = Lines['L2']['dat'].dropna()

    Lines['L3']['dat']['St'] = Lines['L3']['dat']['Stone']
    Lines['L3']['dat'] = Lines['L3']['dat'].dropna()

    interesting_keys = ('L0', 'L1', 'L2', 'L3')
    subdict = {x: Lines[x] for x in interesting_keys if x in Lines}
    for i, L in enumerate(subdict):
        for s in Lines[L]['dat']['St'].unique():
            temp = Lines[L]['dat'].loc[Lines[L]['dat']['St'] == s]
            ax[i].scatter(temp['Date'], temp['dxyz_a'], label=str(s))
        ax[i].grid(axis='both')
        ax[i].set_title(L)
        ax[i].set_ylabel('xyz displacement, m/a')

    ax[1].legend(loc=(1.04, 0))
    ax[3].legend(loc=(1.04, 0))
    ax[0].legend(loc=(-0.3, 0))
    ax[2].legend(loc=(-0.3, 0))
    # ax[0].grid(axis='both')
    fig.savefig('figs/TS_SingleBlocks.png', dpi=150)


# subplots of single blocks
def SingleBlocks(Lines, pts):
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=False, sharey=True)
    ax = ax.flatten()

    Lines['L0']['dat']['St'] = Lines['L0']['dat']['Stone']
    Lines['L0']['dat'] = Lines['L0']['dat'].dropna()

    Lines['L1']['dat']['St'] = Lines['L1']['dat']['Stone'].str.replace(r'\D', '')
    Lines['L1']['dat'] = Lines['L1']['dat'].dropna()

    Lines['L2']['dat']['St'] = Lines['L2']['dat']['Stone']
    Lines['L2']['dat'] = Lines['L2']['dat'].dropna()

    Lines['L3']['dat']['St'] = Lines['L3']['dat']['Stone']
    Lines['L3']['dat'] = Lines['L3']['dat'].dropna()

    yr = Lines['L1']['dat']['Date'].dt.year.unique()

    interesting_keys = ('L0', 'L1', 'L2', 'L3')
    subdict = {x: Lines[x] for x in interesting_keys if x in Lines}

    for i, L in enumerate(subdict):
        # print(i, L)
        for y in yr[yr > 2015]:
            # temp = temp[(temp.Date <'2009-01-01') & (temp.Date >'2001-01-01')]
            temp = Lines[L]['dat'].loc[Lines[L]['dat']['Date'].dt.year == y]
            pts2 = pts[pts.p == 'P'+str(i)+'S']
            dis = np.sqrt(((pts2['x'].values-temp['x'])**2 + (pts2['y'].values-temp['y'])**2))
            # ax[i].scatter(temp['x'], temp['dxyz_a'], label=str(y))
            ax[i].scatter(dis, temp['dxyz_a'], label=str(y))
            # print(temp[['Date','dxyz_a']])
        ax[i].grid(axis='both')
        ax[i].set_title('P'+str(i))
        ax[i].set_ylabel('xyz displacement (m/a)')
        ax[i].set_xlabel('distance from reference point (m)')
    ax[0].legend()
    plt.tight_layout()
    fig.savefig('figs/SingleBlocks.png', dpi=150)


# 3d plot of longditudinal profile and P0
def plot3d(gdf, gdf2):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    temp = gdf2.loc[gdf2['Date'].dt.year > 2015]
    # for s in gdf.Stone.dropna().unique():
    for s in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]: # only use stones available in all years 2016-2021
        ax.scatter(gdf['x'][gdf['Stone'] == s], gdf['y'][gdf['Stone'] == s], gdf['z'][gdf['Stone'] == s], label=str(s))
    # for s in temp.Stone.dropna().unique():
    for s in [1.0, 2.0, 3.0, 5.0]: # only use stones available in all years 2016-2021
        ax.scatter(temp['x'][temp['Stone'] == s], temp['y'][temp['Stone'] == s], temp['z'][temp['Stone'] == s], color='k')

    ax.legend()
    ax.view_init(elev=15., azim=-31)
    fig.savefig('figs/L0_L4_3d.png', dpi=150)


def timeseriesBoth():
    data = pd.read_csv('data/2022_mittelProfile.csv', delimiter=';')
    dat2 = pd.read_csv('data/daten_christoph.csv', delimiter=',')

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(data.year, data.P0, color='k', linestyle='--', marker='s', markersize=4, label='P0')
    ax.plot(data.year, data.P1, color='k', marker='o', markersize=4, label='P1')
    ax.plot(data.year, data.P2, color='k', marker='^', markersize=4, linestyle='-.', label='P2')
    ax.plot(data.year, data.P3, color='k', marker='*', markersize=4, linestyle=':', label='P3')

    ax.hlines(y=dat2.meanv.values, xmin=dat2.y1.values, xmax=dat2.y2.values, linewidth=3, color='grey', label='mean horiz. displacement (Klug et al., 2012)')
    ax.hlines(y=dat2.maxv.values, xmin=dat2.y1.values, xmax=dat2.y2.values, linewidth=3, color='grey', linestyles='--', label='max. horiz. displacement (Klug et al., 2012)')

    ax.grid('both')
    ax.legend()
    ax.set_ylabel('displacement, m/a')
    ax.set_ylim([0, 13])
    fig.savefig('figs/longterm_velocity.png', dpi=150)


def extractFlowline(fname, xs, ys):
    with rasterio.open(fname) as src:
        out = np.asarray([smpl[0] for smpl in src.sample(zip(xs, ys))])
        out[out < 0] = np.nan
    return (out)


def collectFlowline(demsHEK):
    # exctract points from flowline
    ln = gpd.read_file('data/centerline_1.shp')  
    line =ln.geometry[0]
    distance_delta = 1
    distances = np.arange(0, line.length, distance_delta)
    points = [line.interpolate(distance) for distance in distances] + [line.boundary[1]]
    xs = [point.x for point in points]
    ys = [point.y for point in points]
    d = pd.DataFrame(columns=['x', 'y'])
    d['x'] = xs
    d['y'] = ys
    d['dxy'] = np.sqrt((d['x'].diff()**2 + d['y'].diff()**2))
    # extract data from DEMs for each point along flowline
    flowdata = pd.DataFrame(index=d['dxy'].values)
    for i, d in enumerate(demsHEK):
        flowdata[d] = extractFlowline(demsHEK[d], xs, ys)

    ug = extractFlowline('data/DEMs/test_aligned-untergrund1_32632_GEO.tif', xs, ys)

    # hardcoded removal of artefacts in dems.
    ug[1600:] = np.nan

    flowdata[1953][flowdata[1953] < 2400] = np.nan
    flowdata[1953].iloc[0:153] = np.nan
    flowdata[1953].iloc[214:220] = np.nan

    flowdata[1969][flowdata[1969] < 2400] = np.nan
    flowdata[1969][flowdata[1969] > 2800] = np.nan

    temp = flowdata[1969].iloc[1300:]
    temp[temp < 2800] = np.nan
    flowdata[1969].iloc[1300:] = temp
    flowdata[1969].iloc[199:247] = np.nan
    flowdata[1969].iloc[463:483] = np.nan

    flowdata[1971].iloc[1730:] = np.nan

    flowdata[1990][flowdata[1990] < 2380] = np.nan

    flowdata[2019][flowdata[2019] > 1000] = np.nan

    ug[0:266] = flowdata[1953].iloc[0:266]
    # ug[0:152] = flowdata[1971].iloc[0:152]
    ug[0:200] = flowdata[1971].iloc[0:200]


    return(flowdata, ug)


def FlowlineAll(demsHEK):
    flowdata, ug = collectFlowline(demsHEK)
    flowdata.index.values[0] = 0

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(flowdata.index.values.cumsum(), ug, label='Bedrock', color='k', linestyle='--', linewidth=0.5)

    interesting_keys = [1953, 1969, 1971, 1977, 1990, 1997, 2006, 2010, 2017, 2018, 2019, 2020, 2021]
    clr = cm.cividis(np.linspace(0, 1, len(interesting_keys)))

    for i, c in enumerate(flowdata[interesting_keys].columns):
        ax.plot(flowdata.index.values.cumsum(), flowdata[c].values, label=c, color=clr[i])

    ax.grid('both')
    ax.legend()
    ax.set_ylabel('m.a.s.l.')
    ax.set_xlabel('distance along central flowline (m)')
    ax.set_ylim([2300, 2880])
    ax.set_xlim([0, 1780])
    fig.savefig('figs/flowline_all.png')


def FlowlineSubplots(demsHEK):
    flowdata, ug = collectFlowline(demsHEK)
    flowdata.index.values[0] = 0

    interesting_keys = [1953, 1969, 1971, 1977, 1990, 1997, 2006, 2009, 2010, 2017, 2018, 2019, 2020, 2021]
    # clr = cm.cividis(np.linspace(0, 1, len(interesting_keys)))
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    ax = ax.flatten()

    for i, c in enumerate(flowdata[interesting_keys[0:4]].columns):
        clr = cm.cividis(np.linspace(0, 1, len(interesting_keys[0:4])))
        ax[0].plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=0.5, color=clr[i])
    ax[0].plot(flowdata.index.values.cumsum(), ug, label='BR', color='k', linestyle='--', linewidth=0.5)

    for i, c in enumerate(flowdata[interesting_keys[3:6]].columns):
        clr = cm.cividis(np.linspace(0, 1, len(interesting_keys[3:6])))
        ax[1].plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=0.5, color=clr[i])
    ax[1].plot(flowdata.index.values.cumsum(), ug, label='BR', color='k', linestyle='--', linewidth=0.5)

    for i, c in enumerate(flowdata[interesting_keys[5:9]].columns):
        clr = cm.cividis(np.linspace(0, 1, len(interesting_keys[5:9])))
        ax[2].plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=0.5, color=clr[i])
    ax[2].plot(flowdata.index.values.cumsum(), ug, label='BR', color='k', linestyle='--', linewidth=0.5)

    for i, c in enumerate(flowdata[interesting_keys[8:]].columns):
        clr = cm.cividis(np.linspace(0, 1, len(interesting_keys[8:])))
        ax[3].plot(flowdata.index.values.cumsum(), flowdata[c], label=c, linewidth=0.5, color=clr[i])
    ax[3].plot(flowdata.index.values.cumsum(), ug, label='BR', color='k', linestyle='--', linewidth=0.5)

    ax[3].legend()
    ax[3].set_xlim([100,600])
    ax[3].set_ylim([2300,2700])
# ax[3].legend()

    ax[0].set_ylabel('m.a.s.l.')
    ax[2].set_xlabel('distance along central flowline (m)')
    for a in ax:
        a.legend()
        a.grid('both')
    fig.savefig('figs/flowline_sub.png')
