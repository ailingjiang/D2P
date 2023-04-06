# -*- coding: utf-8 -*-
"""
@author: ailing

"""
import os
import numpy as np
import richdem as rd
from Depression_Delineation import np2rdarray, regionGroup, polygonize
import pandas as pd

def Depression_Characterization(dep_id, dep_level, dep_dir):
    
    # Generate emplty array for top deression array and leaf depression array
    no_data = dep_level.no_data
    projection = dep_level.projection
    geotransform = dep_level.geotransform
    ny, nx = dep_level.shape

    depTop = np2rdarray(np.ones((ny,nx))*no_data, no_data, projection, geotransform)
    depLeaf = np2rdarray(np.ones((ny,nx))*no_data, no_data, projection, geotransform)

    depLeaf[dep_level==1] = dep_id[dep_level==1]
    lbl_objects, n_labels = regionGroup(dep_id, 0, no_data)
    depData = pd.read_csv(os.path.join(dep_dir,'depressions_info.csv'))
    
    for label in range(1,n_labels+1):
        depIds = np.unique(np.ma.masked_array(dep_id,mask = lbl_objects!=label).filled(0))
        depData_sub = depData[depData.id.isin(depIds)]
        maxLevelid = depData_sub.id[depData_sub.level.idxmax()]
        depTop[lbl_objects==label] = maxLevelid
    
    rd.SaveGDAL(os.path.join(dep_dir,'top_depressionIDs.tif'), depTop) 
    rd.SaveGDAL(os.path.join(dep_dir,'leaf_depressionIDs.tif'), depLeaf) 
    #rd.rdShow(lbl_objects,figsize=(8,5.5),cmap='jet')
    polygonize(os.path.join(dep_dir,'top_depressionIDs.tif'),os.path.join(dep_dir,'top_depressionIDs.shp'))
    polygonize(os.path.join(dep_dir,'leaf_depressionIDs.tif'),os.path.join(dep_dir,'leaf_depressionIDs.shp'))
    
    return depTop, depLeaf