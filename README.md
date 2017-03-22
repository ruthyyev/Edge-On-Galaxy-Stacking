# Edge-On2

This repository contains SPIRE 250, 350 and 500 micron data for a sample of edge-on galaxies in the DustPedia database. The stacking_pipeline.txt script is used to background subtract, rotate, regrid, crop and stack the maps to produce one final stacked map. 

NB This pipeline is intended for the use of edge-on galaxies, as we use a kpc-per-pixel regime to regrid the maps. This is because d25 does not provide an accurate representation of galaxy size for edge-on galaxies. 

The pipeline is available provided you have SPIRE 250, 350 and/or 500 micron maps of edge-on galaxies, and a .csv file containing source names and co-ordinates (edit table indexes in script according to correct columns). Also edit the index of your galaxy with the largest linear scale (i.e. worst resolution) when re-gridding the maps. Happy stacking! 
