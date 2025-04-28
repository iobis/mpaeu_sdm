# Rasterized version of EUSeaMap

This rasterized version was produced based on the EUSeaMap 2023 (https://emodnet.ec.europa.eu/geonetwork/srv/eng/catalog.search#/metadata/0a1cb988-22de-48b2-8cda-d90947ef77d1).

We focused on two columns: EUNIScomb and EUNIS2019C.

For each EUNIScomb, we got all features available. We then grouped those by the EUNIS2019C (code) and aggregated the shapefile, producing one single layer (POLYGON or MULTIPOLYGON) per code. Then, each feature was rasterized, assigning to each cell the percentage covered by the polygon (using a base raster at 0.05 degrees resolution). Finally, we aggregated the multiple rasters to produce one single raster per EUNIScomb. For that, we assigned, for each cell, the code relative to the EUNIS2019C which had the highest cover of that particular cell.

So, for example, for EUNIScomb A6.11 we have the codes "ME11" "ME12" "ME15" "MF11" "MF12" "MF15" "MG11" "MG12" "MG15". Each one became a raster layer, showing in each cell the percentage covered by the polygon of that EUNIS feature. If a cell was 50% covered by ME11, and only 2% covered by any of the other codes, then the final code for that cell was ME11.

Each EUNIScomb becomes an individual file in the folder `proc-layers` (identified by the number). You can find which is the EUNIScomb code for that layer by looking at the table on file `features_list_EUSeaMap.csv`:

|EUNIS2019C |EUNIScomb | proc_layer| inside_id|EUNIScombD           |EUNIS2019D                        |
|:----------|:---------|----------:|---------:|:--------------------|:---------------------------------|
|ME11       |A6.11     |          1|         1|A6.11: Deep-sea rock |ME11: Arctic upper bathyal rock   |
|ME12       |A6.11     |          1|         2|A6.11: Deep-sea rock |ME12: Atlantic upper bathyal rock |

So, for the layer called `proc_layer_EUComb1.tif`, you can filter in the column `proc_layer` which is equal 1. All the EUNIScomb for that layer will be A6.11 (A6.11: Deep-sea rock). What will change is the `inside_id`, the value you can extract from the layer. If the value you extracted is 2, that means that the EUNIS2019C code is ME12 (ME12: Atlantic upper bathyal rock).

The code to produce this version is at `rasterize-euseamap.R`