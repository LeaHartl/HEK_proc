# HEK_proc

proc_data.py is the main file and calls figures.py and rasterplots.py

## proc_data: 
- reads csv files generated from Martin's excel and turns them into a dictionary of geodataframes for each line
- computes xy and xyz displacement and annual movement
- contains function to write point coordinates with velocities to shapefile
- contains funtion to convert degrees/min/sec to decimal degrees
- sets up dictionary with file names of DEMs (key = year, fname = filename)

## raster_plots:
contains code to plot 
- figure with subplots of hillshades
- overview figure

## figures
- various plots of mean profile velocities and single blocks
- flow line polts and funtcions to extract elevation from flowline
