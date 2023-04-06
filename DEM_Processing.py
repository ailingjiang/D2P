# -*- coding: utf-8 -*-
"""
@author: ailing

"""

## Workflow example
import os
import richdem as rd
import numpy as np
from GRASS_GIS_Processing import GRASS_WatershedAnalysis
from Depression_Delineation import DelineateDepressions_Mod, np2rdarray, polygonize
from CRS_Carving import CRS_Carving
from Depression_Characterization import Depression_Characterization
from osgeo import gdal
from osgeo import ogr
import matplotlib.pyplot as plt


data_dir = './Example_Input_Data/'
rawdem_name = 'GCEW_10m_DEM_MedianF.tif'
ws_name = 'GCEW_10m_Watershed.tif'
obsDepName = 'NWI_wetlands.shp'


## minimum number of pixels to be considered as a depression
min_size = 9

## minimum number of pixels to be considered as a nested depression
min_size_sub = 1

## minimum average depth to be considered as a depression
min_depth = 0.0036

## minimum flow accumulation values to be considered as a river
threshold = 0.001*216656.13529464 # 0.001*max(flowacc)

## minimum slope for non-depression cells
minslope = 1e-3

## maximum slope for entire domain
maxslope = 1.5


out_dir = os.path.join(os.getcwd(), "Output")

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

#%% Import Raw DEM

dem_raw = rd.LoadGDAL(os.path.join(data_dir, rawdem_name))

## Get infomration from the dem
no_data = dem_raw.no_data
projection = dem_raw.projection
geotransform = dem_raw.geotransform
cell_size = np.round(geotransform[1], decimals=2)
ny, nx = dem_raw.shape

## Set value outside the watershed boundary as No Data
if ws_name != None:
    ws_mask =  rd.LoadGDAL(os.path.join(data_dir, ws_name))
    dem_raw = np.ma.masked_array(dem_raw,ws_mask==0, fill_value=no_data)
    dem_raw = dem_raw.filled()

rd.rdShow(dem_raw, ignore_colours=[no_data],figsize=(8,5.5))
#%% Delineate Depressions

####-----------------------------------Fill DEM using Priority Flow with D4 topology------------------------------------------------------------------
print("Priority flow transverse ...")
dem_filled = rd.FillDepressions(dem_raw,topology='D4')
dem_diff =  dem_filled-dem_raw

rd.SaveGDAL(os.path.join(out_dir,'dem_diff.tif'), dem_diff)
rd.SaveGDAL(os.path.join(out_dir,'dem_filled.tif'), dem_filled)


#####----------------------------------Delineate Depressions and embed them into dem------------------------------------------------------------------
fdep_id, fdep_level, fsmoothElev = DelineateDepressions_Mod(os.path.join(out_dir,'dem_filled.tif'),os.path.join(data_dir, rawdem_name), min_size, min_size_sub, min_depth, out_dir, 0)
    
dep_id = rd.LoadGDAL(fdep_id)
smoothElev = rd.LoadGDAL(fsmoothElev)

## Keep real depression in the filled DEM and also smoothed elevation
dem_refined = dem_filled.copy()
dem_refined[dep_id > 0] = dem_raw[dep_id > 0] # reverse to elevations forming depressions
dem_refined[smoothElev > 0] = smoothElev[smoothElev > 0] # use smoothened elevations for depression cells
rd.SaveGDAL(os.path.join(out_dir,'dem_refined.tif'), dem_refined)


#%% Filter out Unwanted Depressions

#####-------------------------------------------Perform Watershed Analysis-------------------------------------------------------------
## r.watershed analysis to obtain flow direction using A* LCP
flwdir_Path, subbasin_Path, flwacc_Path, strseg_Path, strvector_Path = GRASS_WatershedAnalysis(os.path.join(out_dir,'dem_refined.tif'),threshold,out_dir)

if ws_name == None:
    LPC_flrdir = rd.LoadGDAL(flwdir_Path,no_data=255)
else:
    LPC_flrdir = rd.LoadGDAL(flwdir_Path)
    
LPC_subbasin = rd.LoadGDAL(subbasin_Path)
LPC_strseg = rd.LoadGDAL(strseg_Path)

## Load stream seg raster but set no data and change the data type to allow visualization
LPC_strseg_plot = rd.LoadGDAL(strseg_Path,no_data = no_data)
LPC_strseg_plot= LPC_strseg_plot.astype('float32')
if ws_name != None:
    LPC_strseg_plot[ws_mask == 0] = 0
rd.rdShow(LPC_strseg_plot,ignore_colours=[0],figsize=(8,5.5),cmap='jet') 


#####-------------------------------------------Filter Depressions using NWI data------------------------------------------------------------
## Obtain the parent and leaf depressions from depression id and level information
dep_level = rd.LoadGDAL(fdep_level)
dep_id = rd.LoadGDAL(fdep_id)
depTop, depLeaf = Depression_Characterization(dep_id, dep_level, out_dir)

## Find if top depressions and obervations overlaps
srcDS1 = gdal.OpenEx(os.path.join(out_dir,'top_depressionIDs.shp'))
srcDS2 = gdal.OpenEx(os.path.join(data_dir,obsDepName))
lyrDS1 = srcDS1.GetLayer(0)
lyrDS2 = srcDS2.GetLayer(0)
tempDB = os.path.join(out_dir,'temp.splite')
outShp = os.path.join(out_dir,'Dep_Match_NWI.shp')

## Select only depressions larger than the minium size threshold from observation datasets
lyrDS2.SetAttributeFilter("AreaM2 > "+str(min_size*cell_size**2))

## If using other observation data, the following attribute names have to be changed
dtemp = gdal.VectorTranslate(tempDB, srcDS1, format = 'SQLite',layerName = 'DEM_Dep', datasetCreationOptions=['SPATIALITE=YES'])
dtemp = gdal.VectorTranslate(tempDB, srcDS2, format = 'SQLite', accessMode ='update',layerName = 'NWI',SQLDialect='sqlite')
sqlclip = 'SELECT D.id, D.geometry, W.ATTRIBUTE, W.AreaM2 from DEM_Dep D, NWI W '\
    'WHERE ST_INTERSECTS(D.geometry, W.geometry) AND '\
        'D.rowid  IN (SELECT rowid  FROM SpatialIndex WHERE f_table_name = "DEM_Dep" AND search_frame = D.geometry)'
dclip = gdal.VectorTranslate(outShp, gdal.OpenEx(tempDB),SQLDialect='sqlite', SQLStatement=sqlclip ,accessMode = 'overwrite')
del dclip,dtemp
os.remove(tempDB)

## Load matched depressions
shp_match = gdal.OpenEx(outShp,gdal.OF_UPDATE)
lyrmatch = shp_match.GetLayer(0)


## Calculate RtD ratio for the matched depressions
new_field = ogr.FieldDefn('ratio', ogr.OFTReal)
lyrmatch.CreateField(new_field)
obsRatio = []

for feat in lyrmatch:
    did = float(feat.GetFieldAsString('id'))
    idLevel = dep_level[depTop==did]
    idRiverLevel = dep_level[(depTop==did) & (LPC_strseg!=LPC_strseg.no_data)]
    idRiverLeaf = np.unique(depLeaf[(depTop==did) & (LPC_strseg!=LPC_strseg.no_data)])
    idLeaf = np.unique(depLeaf[depTop==did])

    ratio = np.count_nonzero(idRiverLevel==1)/np.count_nonzero(idLevel==1)*idRiverLeaf.size/idLeaf.size
    obsRatio.append(ratio)

    lyrmatch.SetFeature(feat)

## RtD threshold can then be determines as maximum RtD values excluding 1
ratioThes = max(list(filter(lambda a: a != 1, obsRatio)))


## Filter out depressions that are larger than RtD threshold
depTopRefined = depTop.copy()
depList = np.unique(depTop)[1:]
depList = np.sort(depList)

for did in depList:
    idLevel = dep_level[depTop==did]
    idRiverLevel = dep_level[(depTop==did) & (LPC_strseg!=LPC_strseg.no_data)]
    idRiverLeaf = np.unique(depLeaf[(depTop==did) & (LPC_strseg!=LPC_strseg.no_data)])
    idLeaf = np.unique(depLeaf[depTop==did])
    ratio = np.count_nonzero(idRiverLevel==1)/np.count_nonzero(idLevel==1)*idRiverLeaf.size/idLeaf.size
    if (ratio > ratioThes):
        depTopRefined[depTopRefined==did] = no_data

rd.rdShow(depTopRefined,figsize=(8,5.5),ignore_colours=[0],cmap='jet')

rd.SaveGDAL(os.path.join(out_dir,'top_depressionIDs_Refined.tif'), depTopRefined)
polygonize(os.path.join(out_dir,'top_depressionIDs_Refined.tif'),os.path.join(out_dir,'top_depressionIDs_Refined.shp'))

#%% Smoothen River Segments
####-----------------------------------------------CRS Carving (Can also see Schwanghart et al. 2017 for theories behind)--------------------------------------------
# 50th quantile regressions
tau = 0.5
## Degree of smoothening for river segemnt outside depressions
K = 2 
## Degree of smoothening for river segemnt withen depressions
Kdep = 0.08 #
## Minimun number of river pixels to be smoothened together
minStrLen = 100 

dem_RivSmooth =  CRS_Carving(dem_refined,depTopRefined, LPC_strseg,LPC_flrdir,tau,K,Kdep, minStrLen,minslope,maxslope)
rd.SaveGDAL(os.path.join(out_dir,'dem_RivSmooth.tif'), dem_RivSmooth)

####------------------------ Fill the depression after river smoothing ----------------------------------------------
dem_excStream = dem_RivSmooth.copy()
dem_excStream[LPC_strseg!=LPC_strseg.no_data] = dem_excStream.no_data
dem_excStream[depTopRefined>0] =  dem_excStream.no_data
dem_final = rd.FillDepressions(dem_RivSmooth,topology='D4')
dem_final[LPC_strseg!=LPC_strseg.no_data] = dem_RivSmooth[LPC_strseg!=LPC_strseg.no_data]
dem_final[depTopRefined>0] = dem_RivSmooth[depTopRefined>0]
rd.rdShow(dem_final,figsize=(8,5.5),cmap='jet')

rd.SaveGDAL(os.path.join(out_dir,'dem_final.tif'), dem_final) 
rd.SaveGDAL(os.path.join(out_dir,'dem_final_minus_rivSmooth.tif'), dem_final-dem_RivSmooth)     

#%% Calculate Slopes

####-------------------------- Slope calculations ------------------------------------------------------------------------------
slopex = np2rdarray(np.ones((ny,nx))*no_data, no_data, projection, geotransform)
slopey = slopex.copy()
slopex[:,:-1] = (dem_final[:,1:]-dem_final[:,:-1])/cell_size
slopey[1:,:] = (dem_final[:-1,:]-dem_final[1:,:])/cell_size

if ws_name ==  None:
    ws_mask = np.ones((ny,nx))
    slopex[:,-1] = slopex.no_data
    slopey[0,:] = slopey.no_data
else:
    slopex[ws_mask==0] = slopex.no_data
    slopey[ws_mask==0] = slopey.no_data
     
## Loop through all the cells within the provided mask
## Adjust the slopes to fullfill downwinding approach to be consistent with ParFlow's OverlandFlow boundary 
for j in range(0,ny):
    xlist = np.where(ws_mask[j,:] == 1)[0]
    for i in xlist:
        if (i == 0) | (slopex[j,max(i-1,0)] == no_data):
            xlf = i
        else:
            xlf = i-1
            
        if (j == (ny-1)) | (slopey[min(j+1,ny-1),i] == no_data):
            ybf = j
        else:
            ybf = j+1
            
        if slopex[j,min(i+1,nx-1)] == no_data:
            slopex[j,i] = slopex[j,i-1]
            
        if slopey[max(j-1,0),i] == no_data:
            slopey[j,i] = slopey[j+1,i]
            
        # North
        if abs(LPC_flrdir[j,i]) == 2:
            if (slopey[j,i] == 0) & (minslope > 0):
                slopey[j,i] = -minslope
            if ((depTopRefined[j,i]) == depTopRefined.no_data) & (slopey[j,i] > 0):
                slopey[j,i] = -abs(slopey[j,i])
            if minslope > 0:
                slopey[j,i] = np.sign(slopey[j,i])*max(abs(slopey[j,i]),minslope)
                
        # West
        elif abs(LPC_flrdir[j,i]) == 4:
            if (slopex[j,xlf]  == 0) & (minslope > 0):
                slopex[j,xlf] = minslope
            if ((depTopRefined[j,i]) == depTopRefined.no_data) & (slopex[j,xlf] < 0):
                slopex[j,xlf] = abs(slopex[j,xlf])
            if minslope > 0:
                slopex[j,xlf] = np.sign(slopex[j,xlf])*max(abs(slopex[j,xlf]),minslope)
        # South
        elif abs(LPC_flrdir[j,i]) == 6:
            if (slopey[ybf,i] == 0) & (minslope > 0):
                slopey[ybf,i] = minslope
            if ((depTopRefined[j,i]) == depTopRefined.no_data) & (slopey[ybf,i] < 0):
                slopey[ybf,i] = abs(slopey[ybf,i])
            if minslope > 0:
                slopey[ybf,i] = np.sign(slopey[ybf,i])*max(abs(slopey[ybf,i]),minslope)

        # East
        elif abs(LPC_flrdir[j,i]) == 8:
            if (slopex[j,i] == 0) & (minslope > 0):
                slopex[j,i] == -minslope
            if ((depTopRefined[j,i]) == depTopRefined.no_data) & (slopex[j,i] > 0):
                slopex[j,i]  = -abs(slopex[j,i])
            if minslope > 0:
                slopex[j,i] = np.sign(slopex[j,i])*max(abs(slopex[j,i]),minslope)

if maxslope > 0:
    boolx = (abs(slopex) > maxslope) & (slopex != no_data)
    booly = (abs(slopey) > maxslope) & (slopey != no_data)
    slopex[boolx] = np.sign(slopex[boolx])*maxslope
    slopey[booly] = np.sign(slopey[booly])*maxslope

## Process the no neighboring cell in one direction data
boolx = (ws_mask==1) & (slopex==no_data)
slopex[boolx] = 0
booly = (ws_mask==1) & (slopey==no_data)
slopey[booly] = 0

rd.rdShow(slopex,figsize=(8,5.5),cmap='jet')
plt.ylabel('slope x')
rd.rdShow(slopey,figsize=(8,5.5),cmap='jet') 
plt.ylabel('slope y')
rd.SaveGDAL(os.path.join(out_dir,'slopex.tif'), slopex) 
rd.SaveGDAL(os.path.join(out_dir,'slopey.tif'), slopey) 


####--------------------------------------Convert geoTiff to ascii files--------------------------------------------------------------
ds = gdal.Translate(os.path.join(out_dir,'dem_final.asc'),os.path.join(out_dir,'dem_final.tif'),format='AAIGrid')
ds = gdal.Translate(os.path.join(out_dir,'slopex.asc'),os.path.join(out_dir,'slopex.tif'),format='AAIGrid')
ds = gdal.Translate(os.path.join(out_dir,'slopey.asc'),os.path.join(out_dir,'slopey.tif'),format='AAIGrid')